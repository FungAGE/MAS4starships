#!/usr/bin/env python
"""
Annotation Dataset Preparation Script

This script helps prepare annotation data by:
1. Validating existing annotations
2. Identifying missing or incomplete annotations
3. Ensuring consistency across similar sequences
4. Filling in missing information where possible

Uses the mas dev conda environment.
"""

import os
import sys
from dotenv import load_dotenv
import django

# Load environment variables from .env file
load_dotenv()

# Set up Django environment
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'MAS.settings')

# Override database settings for script running outside Docker
if not os.getenv('DOCKER_COMPOSE'):
    os.environ['DB_HOST'] = '127.0.0.1'  # Use localhost when running outside Docker
    os.environ['DB_PORT'] = '3307'  # Use the exposed port from docker-compose

django.setup()

# Now we can import Django models
from starship.models import Annotation, Feature
from starship.views import flag_options_reverse
from MAS.celery import app

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import mysql.connector
from tqdm import tqdm
import argparse
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqRecord
from io import StringIO
import subprocess
import tempfile
from result_viewer.api.tasks import run_single_search
from result_viewer.models import Blastp_Result

class AnnotationPreparation:
    def __init__(self, db_connection, dry_run=False):
        self.connection = db_connection
        self.cursor = self.connection.cursor()
        self.dry_run = dry_run
        self.flag_options = {
            "GREEN": 0,
            "YELLOW": 1,
            "RED": 2,
            "REVIEW NAME": 3,
            "N/A": 4,
            "ORANGE": 5,
            "UNANNOTATED": 7
        }
        # Create temporary BLAST database if it doesn't exist
        self._setup_blast_db()
        
    def _setup_blast_db(self):
        """Create a temporary BLAST database from current annotations"""
        # Create a temporary file for sequences
        self.temp_fasta = tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False)
        
        # Write all sequences to the temporary file
        for ann in Annotation.objects.all():
            if ann.sequence:
                record = SeqRecord.SeqRecord(
                    Seq(ann.sequence),
                    id=str(ann.id),
                    description=f"{ann.accession}|{ann.annotation}"
                )
                SeqIO.write(record, self.temp_fasta, "fasta")
        
        self.temp_fasta.close()
        
        # Create BLAST database
        makeblastdb_cmd = f"makeblastdb -in {self.temp_fasta.name} -dbtype prot"
        subprocess.run(makeblastdb_cmd, shell=True, check=True)

    def validate_annotation(self, annotation):
        """Validate a single annotation entry"""
        issues = []
        
        # Check required fields
        if not annotation.sequence:
            issues.append("Missing sequence")
        if not annotation.annotation:
            issues.append("Missing annotation name")
        if annotation.flag not in self.flag_options.values():
            issues.append(f"Invalid flag value: {annotation.flag}")
            
        # Validate sequence is valid protein sequence and starts with a start codon
        if annotation.sequence:
            if not all(c in "ACDEFGHIKLMNPQRSTVWY*" for c in annotation.sequence.upper()):
                issues.append("Invalid protein sequence characters")
            if not annotation.sequence.upper().startswith("M"):
                issues.append("Protein sequence does not start with a start codon")
                
        return issues

    def find_similar_sequences(self, annotation, threshold=0.9):
        """
        Find annotations with similar sequences using existing BLAST infrastructure
        
        Args:
            annotation (Annotation): Annotation object to find similarities for
            threshold (float): Minimum identity threshold (0-1)
            
        Returns:
            list: List of Annotation objects with similar sequences
        """
        # Check if we already have recent results
        existing_results = Blastp_Result.objects.filter(
            annotation=annotation,
            database='internal',
            status=0  # Completed successfully
        ).order_by('-run_date').first()
        
        if not existing_results:
            # Queue the BLAST search using existing infrastructure
            run_single_search.delay(
                accession=annotation.accession,
                tool='blastp',
                database='internal',
                site='localhost'  # Or appropriate site name
            )
            return []  # Results will be available later
        
        # Parse existing BLAST results
        similar_annotations = []
        blast_results = existing_results.result
        
        for hit in blast_results:
            try:
                hit_identity = float(hit.get('pident', 0)) / 100.0
                if hit_identity >= threshold:
                    hit_annotation = Annotation.objects.get(accession=hit.get('sseqid'))
                    similar_annotations.append(hit_annotation)
            except (Annotation.DoesNotExist, ValueError):
                continue
            
        return similar_annotations

    def standardize_annotation(self, annotation, similar_annotations):
        """Standardize annotation based on similar sequences"""
        if not similar_annotations:
            return annotation
            
        # Get most common annotation name among similar sequences
        annotation_counts = defaultdict(int)
        for similar in similar_annotations:
            annotation_counts[similar.annotation] += 1
            
        if annotation_counts:
            most_common = max(annotation_counts.items(), key=lambda x: x[1])
            if most_common[1] >= len(similar_annotations) * 0.7:  # 70% consensus
                annotation.annotation = most_common[0]
                
        return annotation

    def process_all_annotations(self, annotations=None):
        """Process and validate all annotations in the database"""
        print("Processing annotations..." + (" (DRY RUN)" if self.dry_run else ""))
        
        if annotations is None:
            annotations = Annotation.objects.all()
        
        stats = {
            "total": 0,
            "issues": 0,
            "would_update": 0 if self.dry_run else "updated",
            "standardized": 0
        }
        
        for annotation in tqdm(annotations):
            stats["total"] += 1
            
            # Validate
            issues = self.validate_annotation(annotation)
            if issues:
                stats["issues"] += 1
                tqdm.write(f"Issues with annotation {annotation.accession}:")
                for issue in issues:
                    tqdm.write(f"  - {issue}")
                    
            # Find similar sequences using existing infrastructure
            similar = self.find_similar_sequences(annotation)
            
            # Standardize annotation
            if similar:  # Only standardize if we have results
                updated_annotation = self.standardize_annotation(annotation, similar)
                if updated_annotation != annotation:
                    stats["standardized"] += 1
                    if self.dry_run:
                        stats["would_update"] += 1
                        print(f"\nWould update annotation {annotation.accession}:")
                        print(f"  Old annotation: {annotation.annotation}")
                        print(f"  New annotation: {updated_annotation.annotation}")
                    else:
                        annotation = updated_annotation
                        annotation.save()
                
        return stats

    def __del__(self):
        """Cleanup temporary files"""
        try:
            # Remove BLAST database files
            for ext in ['.phr', '.pin', '.psq']:
                os.remove(self.temp_fasta.name + ext)
            # Remove temporary FASTA file
            os.remove(self.temp_fasta.name)
        except:
            pass

def main():
    # Add argument parsing
    parser = argparse.ArgumentParser(description='Prepare and validate annotation dataset.')
    parser.add_argument('--dry-run', action='store_true', 
                       help='Show what changes would be made without actually making them')
    parser.add_argument('--limit', type=int, 
                       help='Limit the number of annotations to process (for testing)')
    args = parser.parse_args()

    # Database connection using environment variables
    connection = mysql.connector.connect(
        host=os.getenv('DB_HOST', '127.0.0.1'),
        port=int(os.getenv('DB_PORT', '3307')),
        user=os.getenv('DB_USER', 'root'),
        password=os.getenv('DB_PASSWORD', 'changeme'),
        database=os.getenv('DB_NAME', 'mas')
    )
    
    prep = AnnotationPreparation(connection, dry_run=args.dry_run)
    
    # Modify the query to limit annotations if specified
    annotations = Annotation.objects.all()
    if args.limit:
        print(f"Limiting to {args.limit} annotations")
        annotations = annotations[:args.limit]
        
    stats = prep.process_all_annotations(annotations)  # Pass the annotations queryset
    
    print("\nProcessing complete!")
    print(f"Total annotations processed: {stats['total']}")
    print(f"Annotations with issues: {stats['issues']}")
    print(f"Annotations that would be standardized: {stats['standardized']}")
    if args.dry_run:
        print(f"Updates that would be made: {stats['would_update']}")
    else:
        print(f"Database updates made: {stats['updated']}")

if __name__ == "__main__":
    main()