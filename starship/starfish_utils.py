"""
Utilities for managing starfish-nextflow input and output data
"""

import os
import shutil
import pandas as pd
import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from django.conf import settings
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
import logging

from .models import StarfishRun, StarfishRunGenome, StarfishElement

logger = logging.getLogger(__name__)


class StarfishDataManager:
    """Manager class for starfish input/output data operations"""
    
    def __init__(self, run: StarfishRun):
        self.run = run
        self.base_dir = os.path.join(settings.MEDIA_ROOT, 'starfish_runs', run.run_name)
        self.results_dir = os.path.join(self.base_dir, 'results')
        self.work_dir = os.path.join(self.results_dir, 'work')
        
    def setup_directories(self):
        """Create necessary directories for the run"""
        directories = [
            self.base_dir,
            self.results_dir,
            self.work_dir,
            os.path.join(self.results_dir, 'elements'),
            os.path.join(self.results_dir, 'annotations'),
            os.path.join(self.results_dir, 'pairViz'),
            os.path.join(self.results_dir, 'locusViz'),
            os.path.join(self.results_dir, 'summary')
        ]
        
        for directory in directories:
            os.makedirs(directory, exist_ok=True)
            
        logger.info(f"Created directories for run: {self.run.run_name}")
    
    def validate_input_files(self) -> Tuple[bool, List[str]]:
        """Validate that all input files exist and are accessible"""
        errors = []
        
        for genome in self.run.genomes.all():
            # Check required files
            required_files = [
                ('fna_path', genome.fna_path, 'Genome assembly'),
                ('gff3_path', genome.gff3_path, 'GFF3 annotation')
            ]
            
            for field_name, file_path, file_type in required_files:
                if not os.path.exists(file_path):
                    errors.append(f"{file_type} file not found for {genome.genome_id}: {file_path}")
                elif not os.access(file_path, os.R_OK):
                    errors.append(f"Cannot read {file_type} file for {genome.genome_id}: {file_path}")
            
            # Check optional files if provided
            optional_files = [
                ('emapper_path', genome.emapper_path, 'eggNOG mapper annotations'),
                ('cds_path', genome.cds_path, 'CDS sequences'),
                ('faa_path', genome.faa_path, 'Protein sequences')
            ]
            
            for field_name, file_path, file_type in optional_files:
                if file_path and file_path.strip():
                    if not os.path.exists(file_path):
                        errors.append(f"{file_type} file not found for {genome.genome_id}: {file_path}")
                    elif not os.access(file_path, os.R_OK):
                        errors.append(f"Cannot read {file_type} file for {genome.genome_id}: {file_path}")
        
        return len(errors) == 0, errors
    
    def create_samplesheet(self) -> str:
        """Create samplesheet CSV file for the run"""
        samplesheet_path = os.path.join(self.base_dir, 'samplesheet.csv')
        
        data = []
        for genome in self.run.genomes.all():
            row = {
                'genomeID': genome.genome_id,
                'taxID': genome.tax_id or '',
                'fna': genome.fna_path,
                'gff3': genome.gff3_path,
                'emapper': genome.emapper_path or '',
                'cds': genome.cds_path or '',
                'faa': genome.faa_path or ''
            }
            data.append(row)
        
        df = pd.DataFrame(data)
        df.to_csv(samplesheet_path, index=False)
        
        logger.info(f"Created samplesheet for {len(data)} genomes: {samplesheet_path}")
        return samplesheet_path
    
    def copy_input_files(self) -> bool:
        """Copy input files to run directory for processing"""
        try:
            input_dir = os.path.join(self.base_dir, 'input')
            os.makedirs(input_dir, exist_ok=True)
            
            for genome in self.run.genomes.all():
                genome_dir = os.path.join(input_dir, genome.genome_id)
                os.makedirs(genome_dir, exist_ok=True)
                
                # Copy required files
                required_files = [
                    (genome.fna_path, 'genome.fasta'),
                    (genome.gff3_path, 'annotation.gff3')
                ]
                
                for src_path, dst_name in required_files:
                    if os.path.exists(src_path):
                        dst_path = os.path.join(genome_dir, dst_name)
                        shutil.copy2(src_path, dst_path)
                        logger.info(f"Copied {src_path} to {dst_path}")
                
                # Copy optional files
                optional_files = [
                    (genome.emapper_path, 'emapper.annotations'),
                    (genome.cds_path, 'cds.fasta'),
                    (genome.faa_path, 'proteins.faa')
                ]
                
                for src_path, dst_name in optional_files:
                    if src_path and os.path.exists(src_path):
                        dst_path = os.path.join(genome_dir, dst_name)
                        shutil.copy2(src_path, dst_path)
                        logger.info(f"Copied {src_path} to {dst_path}")
            
            return True
            
        except Exception as e:
            logger.error(f"Error copying input files: {str(e)}")
            return False
    
    def parse_results(self) -> Dict[str, int]:
        """Parse starfish results and return statistics"""
        stats = {
            'total_elements': 0,
            'genomes_processed': 0,
            'elements_by_genome': {}
        }
        
        try:
            results_dir = os.path.join(self.results_dir, 'results', self.run.run_name)
            
            if not os.path.exists(results_dir):
                logger.warning(f"Results directory not found: {results_dir}")
                return stats
            
            # Parse elements from regionFinder BED files
            region_finder_dir = os.path.join(results_dir, 'regionFinder')
            if os.path.exists(region_finder_dir):
                for bed_file in os.listdir(region_finder_dir):
                    if bed_file.endswith('.elements.bed'):
                        genome_id = bed_file.replace('.elements.bed', '')
                        elements = self._parse_bed_file(os.path.join(region_finder_dir, bed_file))
                        stats['elements_by_genome'][genome_id] = len(elements)
                        stats['total_elements'] += len(elements)
                        
                        # Create StarfishElement objects
                        self._create_element_objects(elements, genome_id)
            
            # Update genome status
            for genome in self.run.genomes.all():
                genome.num_elements = stats['elements_by_genome'].get(genome.genome_id, 0)
                genome.status = 'completed'
                genome.save()
            
            stats['genomes_processed'] = len(self.run.genomes.all())
            
            # Update run statistics
            self.run.num_elements_found = stats['total_elements']
            self.run.save()
            
            logger.info(f"Parsed {stats['total_elements']} elements from {stats['genomes_processed']} genomes")
            
        except Exception as e:
            logger.error(f"Error parsing results: {str(e)}")
        
        return stats
    
    def _parse_bed_file(self, bed_file_path: str) -> List[Dict]:
        """Parse a BED file and return element data"""
        elements = []
        
        try:
            with open(bed_file_path, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 6:
                            element = {
                                'contig_id': parts[0],
                                'start': int(parts[1]),
                                'end': int(parts[2]),
                                'element_id': parts[3],
                                'score': parts[4] if len(parts) > 4 else '0',
                                'strand': parts[5],
                                'line_number': line_num
                            }
                            elements.append(element)
        except Exception as e:
            logger.error(f"Error parsing BED file {bed_file_path}: {str(e)}")
        
        return elements
    
    def _create_element_objects(self, elements: List[Dict], genome_id: str):
        """Create StarfishElement objects from parsed data"""
        try:
            genome = self.run.genomes.get(genome_id=genome_id)
            
            for element in elements:
                StarfishElement.objects.create(
                    element_id=element['element_id'],
                    run=self.run,
                    genome=genome,
                    contig_id=element['contig_id'],
                    start=element['start'],
                    end=element['end'],
                    strand=element['strand'],
                    sequence='',  # Will be filled from FASTA if available
                    notes=f"Parsed from BED file, line {element['line_number']}"
                )
        except Exception as e:
            logger.error(f"Error creating element objects: {str(e)}")
    
    def get_results_summary(self) -> Dict:
        """Get a summary of results for the run"""
        summary = {
            'run_name': self.run.run_name,
            'status': self.run.status,
            'num_genomes': self.run.num_genomes or 0,
            'num_elements': self.run.num_elements_found or 0,
            'created_at': self.run.created_at,
            'started_at': self.run.started_at,
            'completed_at': self.run.completed_at,
            'duration': str(self.run.duration) if self.run.duration else None,
            'error_message': self.run.error_message
        }
        
        # Add per-genome statistics
        summary['genomes'] = []
        for genome in self.run.genomes.all():
            genome_summary = {
                'genome_id': genome.genome_id,
                'tax_id': genome.tax_id,
                'num_elements': genome.num_elements or 0,
                'status': genome.status
            }
            summary['genomes'].append(genome_summary)
        
        return summary
    
    def export_results(self, format: str = 'json') -> str:
        """Export results in specified format"""
        if format == 'json':
            return self._export_json()
        elif format == 'csv':
            return self._export_csv()
        else:
            raise ValueError(f"Unsupported export format: {format}")
    
    def _export_json(self) -> str:
        """Export results as JSON"""
        results = {
            'run_summary': self.get_results_summary(),
            'elements': []
        }
        
        for element in self.run.elements.all():
            element_data = {
                'element_id': element.element_id,
                'genome_id': element.genome.genome_id,
                'contig_id': element.contig_id,
                'start': element.start,
                'end': element.end,
                'strand': element.strand,
                'length': element.length,
                'family': element.family,
                'navis': element.navis,
                'haplotype': element.haplotype,
                'quality_score': element.quality_score,
                'confidence': element.confidence,
                'created_at': element.created_at.isoformat()
            }
            results['elements'].append(element_data)
        
        output_file = os.path.join(self.results_dir, 'results.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        return output_file
    
    def _export_csv(self) -> str:
        """Export results as CSV"""
        data = []
        for element in self.run.elements.all():
            row = {
                'element_id': element.element_id,
                'genome_id': element.genome.genome_id,
                'contig_id': element.contig_id,
                'start': element.start,
                'end': element.end,
                'strand': element.strand,
                'length': element.length,
                'family': element.family or '',
                'navis': element.navis or '',
                'haplotype': element.haplotype or '',
                'quality_score': element.quality_score or '',
                'confidence': element.confidence or '',
                'created_at': element.created_at.isoformat()
            }
            data.append(row)
        
        df = pd.DataFrame(data)
        output_file = os.path.join(self.results_dir, 'results.csv')
        df.to_csv(output_file, index=False)
        
        return output_file
    
    def cleanup_work_directory(self):
        """Clean up temporary work files"""
        try:
            if os.path.exists(self.work_dir):
                shutil.rmtree(self.work_dir)
                logger.info(f"Cleaned up work directory: {self.work_dir}")
        except Exception as e:
            logger.error(f"Error cleaning up work directory: {str(e)}")
    
    def archive_results(self) -> str:
        """Create a compressed archive of results"""
        try:
            archive_path = os.path.join(self.base_dir, f"{self.run.run_name}_results.tar.gz")
            
            with tarfile.open(archive_path, "w:gz") as tar:
                tar.add(self.results_dir, arcname="results")
            
            logger.info(f"Created results archive: {archive_path}")
            return archive_path
            
        except Exception as e:
            logger.error(f"Error creating archive: {str(e)}")
            return None


class StarfishDataValidator:
    """Validator for starfish input data"""
    
    @staticmethod
    def validate_genome_files(genome_data: Dict) -> Tuple[bool, List[str]]:
        """Validate genome file paths and formats"""
        errors = []
        
        # Check required files
        required_files = ['fna_path', 'gff3_path']
        for field in required_files:
            file_path = genome_data.get(field)
            if not file_path:
                errors.append(f"Missing required field: {field}")
            elif not os.path.exists(file_path):
                errors.append(f"File not found: {file_path}")
            elif not os.access(file_path, os.R_OK):
                errors.append(f"Cannot read file: {file_path}")
        
        # Validate FASTA file format
        fna_path = genome_data.get('fna_path')
        if fna_path and os.path.exists(fna_path):
            if not StarfishDataValidator._is_valid_fasta(fna_path):
                errors.append(f"Invalid FASTA format: {fna_path}")
        
        # Validate GFF3 file format
        gff3_path = genome_data.get('gff3_path')
        if gff3_path and os.path.exists(gff3_path):
            if not StarfishDataValidator._is_valid_gff3(gff3_path):
                errors.append(f"Invalid GFF3 format: {gff3_path}")
        
        return len(errors) == 0, errors
    
    @staticmethod
    def _is_valid_fasta(file_path: str) -> bool:
        """Check if file is valid FASTA format"""
        try:
            from Bio import SeqIO
            with open(file_path, 'r') as f:
                records = list(SeqIO.parse(f, 'fasta'))
                return len(records) > 0
        except:
            return False
    
    @staticmethod
    def _is_valid_gff3(file_path: str) -> bool:
        """Check if file is valid GFF3 format"""
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 9:
                            return True
                return False
        except:
            return False


class StarfishResultsAnalyzer:
    """Analyzer for starfish results"""
    
    def __init__(self, run: StarfishRun):
        self.run = run
    
    def get_element_statistics(self) -> Dict:
        """Get statistics about found elements"""
        elements = self.run.elements.all()
        
        stats = {
            'total_elements': elements.count(),
            'elements_by_genome': {},
            'elements_by_family': {},
            'elements_by_contig': {},
            'length_distribution': {
                'min': 0,
                'max': 0,
                'mean': 0,
                'median': 0
            }
        }
        
        if elements.exists():
            lengths = [elem.length for elem in elements]
            stats['length_distribution'] = {
                'min': min(lengths),
                'max': max(lengths),
                'mean': sum(lengths) / len(lengths),
                'median': sorted(lengths)[len(lengths) // 2]
            }
            
            # Group by genome
            for element in elements:
                genome_id = element.genome.genome_id
                stats['elements_by_genome'][genome_id] = stats['elements_by_genome'].get(genome_id, 0) + 1
                
                # Group by family
                if element.family:
                    stats['elements_by_family'][element.family] = stats['elements_by_family'].get(element.family, 0) + 1
                
                # Group by contig
                stats['elements_by_contig'][element.contig_id] = stats['elements_by_contig'].get(element.contig_id, 0) + 1
        
        return stats
    
    def get_quality_metrics(self) -> Dict:
        """Get quality metrics for the run"""
        elements = self.run.elements.all()
        
        metrics = {
            'total_elements': elements.count(),
            'elements_with_family': elements.filter(family__isnull=False).count(),
            'elements_with_navis': elements.filter(navis__isnull=False).count(),
            'elements_with_haplotype': elements.filter(haplotype__isnull=False).count(),
            'elements_with_quality_score': elements.filter(quality_score__isnull=False).count(),
            'average_quality_score': 0
        }
        
        quality_scores = elements.filter(quality_score__isnull=False).values_list('quality_score', flat=True)
        if quality_scores:
            metrics['average_quality_score'] = sum(quality_scores) / len(quality_scores)
        
        return metrics
