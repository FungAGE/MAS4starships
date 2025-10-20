from tempfile import TemporaryDirectory
import os
import subprocess
from configparser import ConfigParser
import json
import pandas as pd
from datetime import datetime
import shutil
import logging

from celery import shared_task

# from celery.app.control import Inspect
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.db import transaction
from django.conf import settings
from django.db.models.signals import post_save

from starship.genomic_loci_conversions import *
from starship import gene_calling
from starship import models as starship_models
from starship.forms import parse_prots_from_coords
from starship.models import JoinedShips, StarfishRun, StarfishRunGenome, StarfishElement
from MAS.celery import app

logger = logging.getLogger(__name__)


def on_starship_uploaded_or_removed(
    sender, **kwargs
):  # TODO: ADD SIGNAL WHICH CALLS THIS UPON GENOME DELETION
    if settings.NUCLEOTIDE_DB_PATH:
        i = app.control.inspect(settings.CELERY_WORKERS)
        queue = i.scheduled()
        if queue:
            for worker in queue.keys():
                for task in queue[worker]:
                    if (
                        task["request"]["name"]
                        == "starship.tasks.create_internal_nucleotide_blastdb"
                    ):
                        app.control.revoke(task["request"]["id"])

            create_internal_nucleotide_blastdb.apply_async(countdown=5 * 60)


def on_annotation_changed(sender, **kwargs):
    i = app.control.inspect(settings.CELERY_WORKERS)
    queue = i.scheduled()
    if queue:
        for worker in queue.keys():
            for task in queue[worker]:
                if (
                    task["request"]["name"]
                    == "starship.tasks.create_internal_protein_blastdb"
                ):
                    app.control.revoke(task["request"]["id"])

        create_internal_protein_blastdb.apply_async(countdown=5 * 60)


starship_models.starship_upload_complete.connect(on_starship_uploaded_or_removed)
post_save.connect(on_annotation_changed, sender=starship_models.Annotation)
starship_models.starship_upload_complete.connect(on_annotation_changed)


@shared_task
def create_internal_nucleotide_blastdb():
    print("Creating internal nucleotide BLAST database.")
    db_path = settings.NUCLEOTIDE_DB_PATH

    if settings.NUCLEOTIDE_FASTA_PATH:
        fasta_file_path = settings.NUCLEOTIDE_FASTA_PATH
        starships = starship_models.JoinedShips.objects.all()
        starship_list = []
        for starship in starships:
            sequence = SeqRecord(
                Seq(starship.starship_sequence),
                id=starship.starship_name,
                description=starship.starship_name,
            )
            starship_list.append(sequence)

        SeqIO.write(starship_list, fasta_file_path, "fasta")
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                fasta_file_path,
                "-input_type",
                "fasta",
                "-dbtype",
                "nucl",
                "-title",
                "MAS Starships",
                "-out",
                db_path,
            ],
            check=True,
        )


@shared_task
def create_internal_protein_blastdb():
    print("Creating internal protein BLAST database.")

    fasta_file_path = settings.PROTEIN_FASTA_PATH
    db_path = settings.PROTEIN_DB_PATH
    annotations = starship_models.Annotation.objects.all()
    annotations = annotations.prefetch_related("feature_set")
    annotations = annotations.prefetch_related("feature_set__starship")

    annotation_list = []
    for annotation in annotations:
        anno = annotation.annotation
        aa_sequence = annotation.sequence
        public_note = annotation.public_notes
        private_note = annotation.private_notes
        flag = annotation.get_flag_display()
        genomes = ", ".join(
            {feature.starship.starship_name for feature in annotation.feature_set.all()}
        )
        record = SeqRecord(
            Seq(aa_sequence),
            id=annotation.accession + " |",
            description="%s | %s | %s | %s | %s"
            % (anno, public_note, private_note, flag, genomes),
        )
        annotation_list.append(record)
    print("writing fasta")
    SeqIO.write(annotation_list, fasta_file_path, "fasta")
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            '"%s"' % fasta_file_path,
            "-input_type",
            "fasta",
            "-dbtype",
            "prot",
            "-parse_seqids",
            "-title",
            "MAS Internal Protein Database",
            "-out",
            db_path,
        ],
        check=True,
    )


# Create new annotation objects and associated features for CDS
def create_CDS_annotations(
    glimmer_results_file_path, genome, assign_to, new_annotations, new_features
):
    for cds in gene_calling.parse_glimmer_results(glimmer_results_file_path):
        sequence = Seq(JoinedShips.objects.get(starship_name=genome).starship_sequence)
        protein = get_protein_sequence(cds.start, cds.stop, cds.strand, sequence)

        if starship_models.Annotation.objects.filter(sequence=protein).count() > 0:
            new_annotation = starship_models.Annotation.objects.get(sequence=protein)
        elif protein in new_annotations:
            new_annotation = new_annotations[protein]
        else:
            new_annotation = starship_models.Annotation()
            new_annotation.sequence = protein
            new_annotation.assigned_to = assign_to
            new_annotations[protein] = new_annotation

        new_feature = starship_models.Feature()
        new_feature.genome = genome
        new_feature.start = cds.start
        new_feature.stop = cds.stop
        new_feature.type = "CDS"
        new_feature.strand = cds.strand
        new_feature.annotation = new_annotation
        new_features.append(new_feature)


# Create an annotation and feature objects for given coordinate file
def create_custom_CDS_annotations(
    coordinate_file, translation_table, genome, assign_to, new_annotations, new_features
):
    starship_sequence = Seq(JoinedShips.objects.get(starship_name=genome).starship_sequence)

    for protein_sequence, cds in parse_prots_from_coords(
        coordinate_file, starship_sequence, translation_table
    ):

        if (
            starship_models.Annotation.objects.filter(sequence=protein_sequence).count()
            > 0
        ):
            new_annotation = starship_models.Annotation.objects.get(
                sequence=protein_sequence
            )

        elif protein_sequence in new_annotations:
            new_annotation = new_annotations[protein_sequence]

        else:
            new_annotation = starship_models.Annotation()
            new_annotation.sequence = protein_sequence
            new_annotation.assigned_to = assign_to
            new_annotations[protein_sequence] = new_annotation

        new_feature = starship_models.Feature()
        new_feature.genome = genome
        new_feature.start = cds.start
        new_feature.stop = cds.stop
        new_feature.type = "CDS"
        new_feature.strand = cds.strand
        new_feature.annotation = new_annotation
        new_features.append(new_feature)


# Create new annotation objects and associated features for tRNAs
def create_trna_annotations(
    trnascan_results_file_path, genome, assign_to, new_annotations, new_features
):
    for tRNA in gene_calling.parse_trnascan_results(trnascan_results_file_path):
        sequence = Seq(JoinedShips.objects.get(starship_name=genome).starship_sequence)
        rna = get_rna_sequence(tRNA.start, tRNA.stop, tRNA.strand, sequence)

        if starship_models.Annotation.objects.filter(sequence=rna).count() > 0:
            new_annotation = starship_models.Annotation.objects.get(sequence=rna)
        elif rna in new_annotations:
            new_annotation = new_annotations[rna]
        else:
            new_annotation = starship_models.Annotation()
            new_annotation.annotation = "tRNA-%s-%s" % (tRNA.amino, tRNA.codon)
            new_annotation.public_notes = (
                "Called with tRNAscan-SE, Inf score = %f" % tRNA.score
            )
            new_annotation.private_notes = (
                "This annotation was automatically generated."
            )
            new_annotation.sequence = rna
            new_annotation.flag = 8
            new_annotation.assigned_to = assign_to
            new_annotations[rna] = new_annotation

        new_feature = starship_models.Feature()
        new_feature.genome = genome
        new_feature.start = tRNA.start
        new_feature.stop = tRNA.stop
        new_feature.type = "tRNA"
        new_feature.strand = tRNA.strand
        new_feature.annotation = new_annotation
        new_features.append(new_feature)


def add_annotations_and_features_to_db(new_annotations, new_features):
    # Adding objects to db gives them an id
    starship_models.Annotation.objects.bulk_create(list(new_annotations.values()))

    # We must assign newly created id in feature object before creation
    for feat in new_features:
        feat.annotation_id = feat.annotation.id

    starship_models.Feature.objects.bulk_create(new_features)


# Starfish-specific tasks
@shared_task(bind=True)
def run_starfish_pipeline(self, run_id, resume=False):
    """
    Celery task to run the starfish-nextflow pipeline
    
    Args:
        run_id: ID of the StarfishRun
        resume: If True, use Nextflow's -resume flag to continue from previous run
    """
    run = None
    try:
        # Get the run object
        run = StarfishRun.objects.get(id=run_id)
        
        # Check if run is already running (unless we're resuming)
        if run.status == 'running' and not resume:
            logger.warning(f"Run {run.run_name} is already running, skipping")
            return
        
        # For resume, verify the run was previously failed or cancelled
        if resume and run.status not in ['failed', 'cancelled', 'running']:
            logger.error(f"Cannot resume run {run.run_name} with status {run.status}")
            return
        
        # Update status to running
        run.status = 'running'
        run.started_at = datetime.now()
        run.celery_task_id = self.request.id
        run.save()
        
        if resume:
            logger.info(f"Resuming starfish pipeline for run: {run.run_name} (using -resume flag)")
        else:
            logger.info(f"Starting starfish pipeline for run: {run.run_name}")
        
        # Set up paths using run ID for uniqueness
        run_dir_name = f"{run.id}_{run.run_name}"
        base_dir = os.path.join(settings.MEDIA_ROOT, 'starfish_runs', run_dir_name)
        os.makedirs(base_dir, exist_ok=True)
        
        # Update paths in the run object
        run.samplesheet_path = os.path.join(base_dir, 'samplesheet.csv')
        run.output_dir = os.path.join(base_dir, 'results')
        run.log_file = os.path.join(base_dir, 'starfish.log')
        run.save()
        
        # Create samplesheet from genomes
        create_samplesheet(run)
        
        # Build nextflow command
        nextflow_cmd = build_nextflow_command(run, resume=resume)
                
        # Run the pipeline
        # When resuming, append to log file instead of overwriting
        log_mode = 'a' if resume else 'w'
        with open(run.log_file, log_mode) as log_f:
            if resume:
                log_f.write(f"\n\n{'='*80}\n")
                log_f.write(f"RESUMING RUN AT {datetime.now()}\n")
                log_f.write(f"{'='*80}\n\n")
            
            process = subprocess.run(
                nextflow_cmd,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                cwd=base_dir,
                check=False
            )
        
        # Check if pipeline succeeded
        if process.returncode == 0:
            # Verify that all expected output files and directories exist
            if verify_starfish_outputs(run):
                run.status = 'completed'
                run.completed_at = datetime.now()
                
                # Parse results
                parse_starfish_results(run)
                
                logger.info(f"Starfish pipeline completed successfully for run: {run.run_name}")
            else:
                run.status = 'failed'
                run.error_message = "Pipeline completed but expected output files are missing or empty"
                logger.error(f"Starfish pipeline outputs incomplete for run: {run.run_name}")
        else:
            run.status = 'failed'
            run.error_message = f"Pipeline failed with return code {process.returncode}"
            logger.error(f"Starfish pipeline failed for run: {run.run_name}")
            
    except Exception as e:
        logger.error(f"Error in starfish pipeline for run {run_id}: {str(e)}")
        if run is not None:
            run.status = 'failed'
            run.error_message = str(e)
        
    finally:
        if run is not None:
            run.save()

def verify_starfish_outputs(run):
    """Verify that all expected starfish output files and directories exist and are not empty"""
    # TODO: detect status from nextflow output logs, in addition to manual checks?
    results_dir = os.path.join(run.output_dir, 'results', run.run_name)
    
    if not os.path.exists(results_dir):
        logger.error(f"Results directory not found: {results_dir}")
        return False
    
    # Check required files
    required_files = [
        # elementFinder outputs
        '*.insert.bed',
        '*.insert.stats',
        '*.flank.bed',
        '*.flank.singleDR.stats',
        # regionFinder outputs
        '*.elements.bed',
        '*.elements.feat',
        # '*.elements.fna',
        '*.elements.named.stats',
    ]
    
    # Check required directories
    required_dirs = [
        'pairViz'
        # 'locusViz'
    ]
    
    # Verify files exist and are not empty
    for file in required_files:
        
        # Check for files matching pattern
        import glob
        matching_files = glob.glob(file)
        if not matching_files:
            logger.error(f"No files found matching pattern {file}")
            return False
        
        # Check that files are not empty
        for file_path in matching_files:
            if os.path.getsize(file_path) == 0:
                logger.error(f"File is empty: {file_path}")
                return False
    
    # Verify directories exist and are not empty
    for dir_name in required_dirs:
        dir_path = os.path.join(results_dir, dir_name)
        if not os.path.exists(dir_path):
            logger.error(f"Required directory missing: {dir_path}")
            return False
        
        # Check if directory is empty
        if not os.listdir(dir_path):
            logger.error(f"Required directory is empty: {dir_path}")
            return False
    
    logger.info(f"All required starfish outputs verified for run: {run.run_name}")
    return True


def create_samplesheet(run):
    """Create samplesheet CSV file from run genomes"""
    genomes = run.genomes.all()
    
    samplesheet_data = []
    for genome in genomes:
        row = {
            'genomeID': genome.genome_id,
            'taxID': genome.tax_id or '',
            'fna': genome.fna_path,
            'gff3': genome.gff3_path,
            'emapper': genome.emapper_path or '',
            'cds': genome.cds_path or '',
            'faa': genome.faa_path or ''
        }
        samplesheet_data.append(row)
    
    # Create DataFrame and save as CSV
    df = pd.DataFrame(samplesheet_data)
    df.to_csv(run.samplesheet_path, index=False)
    
    # Update run with number of genomes
    run.num_genomes = len(genomes)
    run.save()

def build_nextflow_command(run, resume=False):
    """Build the nextflow command for running starfish
    
    Args:
        run: StarfishRun instance
        resume: If True, add -resume flag to continue from previous run
    """
    pipeline_path = os.getenv('STARFISH_NEXTFLOW_PATH',
                              '/mnt/sda/johannesson_lab/adrian/starfish_pipeline/starfish-nextflow'
    )
    main_nf = os.path.join(pipeline_path, 'main.nf')
    
    # Add -resume flag if requested
    resume_flag = '-resume ' if resume else ''
        
    # Build command that:
    # 1. Saves the current JAVA_CMD and JAVA_HOME
    # 2. Activates starfish environment
    # 3. Restores Java paths so nextflow can run
    # 4. Runs nextflow (which will use starfish env for its processes)
    cmd_string = (
        f'SAVED_JAVA_CMD=$JAVA_CMD && '
        f'SAVED_JAVA_HOME=$JAVA_HOME && '
        f'nextflow run {main_nf} '
        f'{resume_flag}'
        f'-profile local '
        f'--samplesheet {run.samplesheet_path} '
        f'--run_name {run.run_name} '
        f'--model {run.model} '
        f'--threads {run.threads} '
        f'--missing {run.missing} '
        f'--maxcopy {run.maxcopy} '
        f'--pid {run.pid} '
        f'--hsp {run.hsp} '
        f'--flank {run.flank} '
        f'--neighbourhood {run.neighbourhood} '
        f'-w {os.path.join(run.output_dir, "work")} '
        f'--outdir {run.output_dir}'
    )
    
    cmd = ['bash', '-c', cmd_string]
    
    return cmd


def parse_starfish_results(run):
    """Parse starfish results and populate database"""
    try:
        results_dir = os.path.join(run.output_dir, 'results', run.run_name)
        
        # Count total elements found
        total_elements = 0
        
        # Process each genome's results from regionFinder outputs
        for genome in run.genomes.all():
            genome_elements = 0
            
            # Look for genome-specific result files in regionFinder directory
            region_finder_dir = os.path.join(results_dir, 'regionFinder')
            if os.path.exists(region_finder_dir):
                for bed_file in os.listdir(region_finder_dir):
                    if bed_file.endswith('.elements.bed') and genome.genome_id in bed_file:
                        elements = parse_bed_file(os.path.join(region_finder_dir, bed_file), run, genome)
                        genome_elements += len(elements)
                        total_elements += len(elements)
            
            # Update genome with element count
            genome.num_elements = genome_elements
            genome.status = 'completed'
            genome.save()
        
        # Update run with total elements found
        run.num_elements_found = total_elements
        run.save()
        
        logger.info(f"Parsed {total_elements} elements for run {run.run_name}")
        
    except Exception as e:
        logger.error(f"Error parsing starfish results for run {run.run_name}: {str(e)}")
        run.error_message = f"Error parsing results: {str(e)}"
        run.save()


def parse_bed_file(bed_file_path, run, genome):
    """Parse BED file and create StarfishElement objects"""
    elements = []
    
    try:
        with open(bed_file_path, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        contig_id = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        element_id = parts[3]
                        strand = parts[5]
                        
                        # Create StarfishElement
                        element = StarfishElement.objects.create(
                            element_id=element_id,
                            run=run,
                            genome=genome,
                            contig_id=contig_id,
                            start=start,
                            end=end,
                            strand=strand,
                            sequence='',  # Will be filled from FASTA file if available
                            notes=f'Parsed from {os.path.basename(bed_file_path)}'
                        )
                        elements.append(element)
    
    except Exception as e:
        logger.error(f"Error parsing BED file {bed_file_path}: {str(e)}")
    
    return elements


@shared_task
def cleanup_starfish_run(run_id):
    """Clean up temporary files for a completed starfish run"""
    try:
        run = StarfishRun.objects.get(id=run_id)
        
        # Only clean up if run is completed or failed
        if run.status in ['completed', 'failed']:
            # Clean up work directory
            work_dir = os.path.join(run.output_dir, 'work')
            if os.path.exists(work_dir):
                shutil.rmtree(work_dir)
                logger.info(f"Cleaned up work directory for run: {run.run_name}")
        
    except Exception as e:
        logger.error(f"Error cleaning up run {run_id}: {str(e)}")


@shared_task
def import_starfish_elements_to_mas(run_id):
    """Import starfish elements into the main MAS database"""
    try:
        run = StarfishRun.objects.get(id=run_id)
        
        if run.status != 'completed':
            logger.warning(f"Cannot import elements from incomplete run: {run.run_name}")
            return
        
        imported_count = 0
        
        for element in run.elements.all():
            try:
                # Create or get accession
                accession, created = starship_models.Accessions.objects.get_or_create(
                    accession_tag=f"{run.run_name}_{element.element_id}",
                    defaults={
                        'ship_name': element.element_id,
                        'version_tag': run.run_name
                    }
                )
                
                # Create ship if sequence is available
                if element.sequence:
                    ship, created = starship_models.Ships.objects.get_or_create(
                        accession=accession,
                        defaults={
                            'sequence': element.sequence
                        }
                    )
                
                # Create joined ship record
                joined_ship, created = starship_models.JoinedShips.objects.get_or_create(
                    starshipID=element.element_id,
                    defaults={
                        'contigID': element.contig_id,
                        'elementBegin': element.start,
                        'elementEnd': element.end,
                        'ship_family': None,  # Will be filled later
                        'genome': None,  # Will be linked to existing genome if found
                        'source': 'starfish-nextflow',
                        'evidence': 'computational'
                    }
                )
                
                if created:
                    imported_count += 1
                    logger.info(f"Imported element: {element.element_id}")
                
            except Exception as e:
                logger.error(f"Error importing element {element.element_id}: {str(e)}")
        
        logger.info(f"Imported {imported_count} elements from run {run.run_name}")
        
    except Exception as e:
        logger.error(f"Error importing elements for run {run_id}: {str(e)}")
