from django.db import models
import re
from django.core.exceptions import ValidationError
from django.shortcuts import redirect
from django.contrib.auth.models import User
from simple_history.models import HistoricalRecords
from django.conf import settings
import django.dispatch
from django.db.models import Count, Q, F
from django.utils import timezone

from starship.starbase_models import JoinedShips, Accessions, Ships


# Validators
def validate_starship_name(value):
    if settings.STARSHIP_NAME_FORMAT:
        pattern = re.compile(settings.STARSHIP_NAME_FORMAT)
        if pattern.match(value) is None:
            raise ValidationError(
                "%s is not a valid Starship name. Would You like to continue with this name?"
                % value,
                code=1,
            )


def validate_starship_error(option):
    if option != "Yes":
        raise ValidationError("Starship was not accepted", code=2)
    else:
        redirect("starship:starship_list")


def validate_duplicate_name(name):
    starships = JoinedShips.objects.filter(starshipID=name)
    if starships.exists():
        raise ValidationError("%s is already a Starship name." % name, code=3)

starship_upload_complete = django.dispatch.Signal()


# TODO: add features of interest here?
class Feature(models.Model):
    class Meta:
        db_table = 'starship_feature'

    feature_options = (
        ("gene", "Gene Annotation"),
        ("CDS", "Coding Sequence"),
        ("Repeat Region", "Repeat Region"),
        ("tRNA", "tRNA"),
    )

    """
    on_delete=models.CASCADE() - deletes features of a Starship when a Starship is deleted
    """
    starship = models.ForeignKey(
        JoinedShips,
        on_delete=models.CASCADE,
        db_column='starship_id'
    )

    start = models.IntegerField(default=0)
    stop = models.IntegerField(default=0)

    type = models.CharField(max_length=50, choices=feature_options)

    strand = models.CharField(max_length=1)

    # obtain information from the Starship annotations
    annotation = models.ForeignKey(
        "Annotation", blank=True, null=True, on_delete=models.PROTECT
    )

    def __str__(self):
        return "%s: %s %s..%s %s" % (
            self.type,
            self.starship,
            self.start,
            self.stop,
            self.strand,
        )


class Annotation(models.Model):
    class Meta:
        db_table = 'starship_annotation'

    # if you update this, you need to update the flag dict in confirm_upload_annotations in starship views
    flag_options = (
        (0, "GREEN"),
        (1, "YELLOW"),
        (2, "RED"),
        (3, "REVIEW NAME"),
        (4, "N/A"),
        (5, "ORANGE"),
        (7, "UNANNOTATED"),
    )
    annotation = models.CharField(
        max_length=255, blank=True, null=True, default="No Annotation"
    )
    # Amino Acid Sequence
    sequence = models.TextField(max_length=10000, unique=True)
    public_notes = models.TextField(max_length=30000, blank=True, null=True, default="")
    private_notes = models.TextField(
        max_length=30000, blank=True, null=True, default=""
    )
    # default 7 for unannotated Starship
    flag = models.IntegerField(default=7, choices=flag_options)
    assigned_to = models.ForeignKey(User, on_delete=models.SET_NULL, null=True)
    history = HistoricalRecords()

    @property
    def accession(self):
        chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        integer = abs(self.id)
        result = ""
        padding = []

        while integer > 0:
            integer, remainder = divmod(integer, 36)
            padding.insert(0, chars[remainder])
        if len(padding) < 5:
            while len(padding) < 5:
                padding.insert(0, "0")
        for i in padding:
            result = result + i
        return result

    def __str__(self):
        return "%s | %s" % (self.annotation, self.get_flag_display())


class StarfishRun(models.Model):
    """Model to track starfish-nextflow pipeline runs"""
    class Meta:
        db_table = 'starfish_run'
        ordering = ['-created_at']
    
    # Run status choices
    STATUS_CHOICES = (
        ('pending', 'Pending'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
        ('cancelled', 'Cancelled'),
    )
    
    # Run name and description
    run_name = models.CharField(max_length=255, unique=True)
    description = models.TextField(blank=True, null=True)
    
    # Status and timing
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    created_at = models.DateTimeField(auto_now_add=True)
    started_at = models.DateTimeField(null=True, blank=True)
    completed_at = models.DateTimeField(null=True, blank=True)
    
    # User who created the run
    created_by = models.ForeignKey(User, on_delete=models.CASCADE, related_name='starfish_runs')
    
    # Pipeline parameters
    model = models.CharField(max_length=50, default='tyr', help_text='Gene model for de novo annotations')
    threads = models.IntegerField(default=20, help_text='Number of CPU threads')
    missing = models.IntegerField(default=1, help_text='Maximum missing genes in orthogroup filtering')
    maxcopy = models.IntegerField(default=5, help_text='Maximum copy number in orthogroup filtering')
    pid = models.IntegerField(default=90, help_text='Minimum percent identity for BLAST')
    hsp = models.IntegerField(default=1000, help_text='Minimum HSP length for BLAST')
    flank = models.IntegerField(default=6, help_text='Flank size for repeat detection')
    neighbourhood = models.IntegerField(default=10000, help_text='Neighborhood size for sourmash sketch')
    
    # File paths
    samplesheet_path = models.CharField(max_length=500, help_text='Path to samplesheet CSV file')
    output_dir = models.CharField(max_length=500, help_text='Output directory for results')
    log_file = models.CharField(max_length=500, blank=True, null=True, help_text='Path to log file')
    
    # Celery task tracking
    celery_task_id = models.CharField(max_length=255, blank=True, null=True)
    
    # Error information
    error_message = models.TextField(blank=True, null=True)
    
    # Results summary
    num_genomes = models.IntegerField(null=True, blank=True, help_text='Number of genomes processed')
    num_elements_found = models.IntegerField(null=True, blank=True, help_text='Number of starship elements found')
    
    def __str__(self):
        return f"StarfishRun: {self.run_name} ({self.status})"
    
    @property
    def duration(self):
        """Calculate run duration if completed"""
        if self.started_at and self.completed_at:
            return self.completed_at - self.started_at
        return None
    
    def get_absolute_url(self):
        from django.urls import reverse
        return reverse('starship:starfish_run_detail', kwargs={'pk': self.pk})


class StarfishRunGenome(models.Model):
    """Model to track individual genomes in a starfish run"""
    class Meta:
        db_table = 'starfish_run_genome'
        unique_together = ['run', 'genome_id']
    
    run = models.ForeignKey(StarfishRun, on_delete=models.CASCADE, related_name='genomes')
    genome_id = models.CharField(max_length=255, help_text='Unique identifier for the genome')
    tax_id = models.CharField(max_length=50, blank=True, null=True, help_text='Taxonomy ID')
    fna_path = models.CharField(max_length=500, help_text='Path to genome assembly file')
    gff3_path = models.CharField(max_length=500, help_text='Path to GFF3 annotation file')
    emapper_path = models.CharField(max_length=500, blank=True, null=True, help_text='Path to emapper annotations')
    cds_path = models.CharField(max_length=500, blank=True, null=True, help_text='Path to CDS sequences')
    faa_path = models.CharField(max_length=500, blank=True, null=True, help_text='Path to protein sequences')
    
    # Results
    num_elements = models.IntegerField(null=True, blank=True, help_text='Number of elements found in this genome')
    status = models.CharField(max_length=20, default='pending', choices=StarfishRun.STATUS_CHOICES)
    error_message = models.TextField(blank=True, null=True)
    
    def __str__(self):
        return f"{self.run.run_name} - {self.genome_id}"


class StarfishElement(models.Model):
    """Model to store starship elements found by starfish-nextflow"""
    class Meta:
        db_table = 'starfish_element'
        ordering = ['-created_at']
    
    # Element identification
    element_id = models.CharField(max_length=255, unique=True)
    run = models.ForeignKey(StarfishRun, on_delete=models.CASCADE, related_name='elements')
    genome = models.ForeignKey(StarfishRunGenome, on_delete=models.CASCADE, related_name='elements')
    
    # Element coordinates and sequence
    contig_id = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.CharField(max_length=1)
    sequence = models.TextField(help_text='Element sequence')
    
    # Element classification
    family = models.CharField(max_length=255, blank=True, null=True)
    navis = models.CharField(max_length=255, blank=True, null=True)
    haplotype = models.CharField(max_length=255, blank=True, null=True)
    
    # Quality metrics
    quality_score = models.FloatField(null=True, blank=True)
    confidence = models.CharField(max_length=50, blank=True, null=True)
    
    # Metadata
    created_at = models.DateTimeField(auto_now_add=True)
    notes = models.TextField(blank=True, null=True)
    
    def __str__(self):
        return f"Element {self.element_id} from {self.genome.genome_id}"
    
    @property
    def length(self):
        return self.end - self.start + 1


class StagingStarship(models.Model):
    """Staging model for starship submissions before approval"""
    class Meta:
        db_table = 'staging_starship'
        ordering = ['-submitted_at']
    
    # Submission status choices
    STATUS_CHOICES = (
        ('pending', 'Pending Review'),
        ('approved', 'Approved'),
        ('rejected', 'Rejected'),
        ('needs_revision', 'Needs Revision'),
    )
    
    # Basic submission data
    starshipID = models.CharField(max_length=255, help_text='Starship name/ID')
    sequence = models.TextField(help_text='Starship sequence')
    species = models.CharField(max_length=255, help_text='Species name')
    
    # Optional additional data
    evidence = models.CharField(max_length=255, blank=True, null=True, help_text='Evidence for this starship')
    source = models.CharField(max_length=255, blank=True, null=True, help_text='Source of the data')
    notes = models.TextField(blank=True, null=True, help_text='Additional notes')
    
    # File uploads (stored as text for simplicity)
    annotation_file_content = models.TextField(blank=True, null=True, help_text='Annotation file content')
    annotation_file_name = models.CharField(max_length=255, blank=True, null=True, help_text='Original annotation file name')
    
    # Terminal repeat information
    terminal_repeat_length = models.IntegerField(blank=True, null=True, help_text='Length of terminal repeat')
    terminal_repeat_sequence = models.TextField(blank=True, null=True, help_text='Terminal repeat sequence')
    
    # Submission tracking
    submitted_by = models.ForeignKey(User, on_delete=models.CASCADE, related_name='staging_submissions')
    submitted_at = models.DateTimeField(auto_now_add=True)
    reviewed_by = models.ForeignKey(User, on_delete=models.SET_NULL, null=True, blank=True, related_name='reviewed_submissions')
    reviewed_at = models.DateTimeField(null=True, blank=True)
    status = models.CharField(max_length=20, choices=STATUS_CHOICES, default='pending')
    review_notes = models.TextField(blank=True, null=True, help_text='Reviewer notes')
    
    # Validation results
    validation_passed = models.BooleanField(default=False)
    validation_errors = models.TextField(blank=True, null=True, help_text='Validation error messages')
    
    def __str__(self):
        return f"Staging: {self.starshipID} ({self.status})"
    
    def approve(self, reviewer, notes=""):
        """Approve this submission"""
        self.status = 'approved'
        self.reviewed_by = reviewer
        self.reviewed_at = timezone.now()
        self.review_notes = notes
        self.save()
    
    def reject(self, reviewer, notes=""):
        """Reject this submission"""
        self.status = 'rejected'
        self.reviewed_by = reviewer
        self.reviewed_at = timezone.now()
        self.review_notes = notes
        self.save()
    
    def request_revision(self, reviewer, notes=""):
        """Request revision of this submission"""
        self.status = 'needs_revision'
        self.reviewed_by = reviewer
        self.reviewed_at = timezone.now()
        self.review_notes = notes
        self.save()
    
    def migrate_to_main_database(self):
        """Migrate this staging entry to the main JoinedShips database"""
        try:
            # Check if this starshipID already exists in JoinedShips
            if JoinedShips.objects.filter(starshipID=self.starshipID).exists():
                raise ValidationError("This starshipID already exists in the main database.")

            # TODO: implement logic for detecting duplicate sequences here, before assigning accessions
            print(f"Creating Accessions for {self.starshipID}")
            # Create Accessions object first
            accession = Accessions(
                ship_name=self.starshipID,
                accession_tag=f"{self.starshipID}_{timezone.now().strftime('%Y%m%d_%H%M%S')}",
                version_tag="1.0"
            )
            accession.save()
            print(f"Accession created with ID: {accession.id}")
            
            print(f"Creating Ships for {self.starshipID}")
            # Create Ships object with the sequence
            ship = Ships(
                sequence=self.sequence,
                accession=accession
            )
            ship.save()
            print(f"Ship created with ID: {ship.id}")
            
            print(f"Creating JoinedShips for {self.starshipID}")
            # Create the main JoinedShips entry
            main_starship = JoinedShips(
                starshipID=self.starshipID,
                evidence=self.evidence,
                source=self.source,
                curated_status='staged_import',
                ship=ship
            )
            main_starship.save()
            print(f"JoinedShips created with ID: {main_starship.id}")
            
            return main_starship
            
        except Exception as e:
            print(f"Error in migrate_to_main_database: {str(e)}")
            import traceback
            print(traceback.format_exc())
            raise e
