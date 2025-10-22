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
    starships = JoinedShips.objects.all()
    for starship in starships:
        if starship.starshipID == name:
            raise ValidationError("%s is already a Starship name." % name, code=3)


starship_upload_complete = django.dispatch.Signal()


# Models to match schema structure of starbase SQLite
# having the schema/tables matching will be useful when we develop mariadb -> sqlite migration
# TODO: the data displayed in tables should be created from these updated models
class Accessions(models.Model):
    """Accessions table from SQLite schema"""
    class Meta:
        db_table = 'starbase_accessions'
        app_label = 'starship'
    
    ship_name = models.CharField(max_length=255, null=True, blank=True)
    accession_tag = models.CharField(max_length=255, unique=True)
    version_tag = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"{self.accession_tag} - {self.ship_name}"


class Ships(models.Model):
    """Ships table from SQLite schema"""
    class Meta:
        db_table = 'starbase_ships'
    
    sequence = models.TextField()
    md5 = models.CharField(max_length=32, null=True, blank=True)
    rev_comp_md5 = models.CharField(max_length=32, null=True, blank=True)
    accession = models.ForeignKey(Accessions, on_delete=models.CASCADE, related_name='ships')

    def __str__(self):
        return f"Ship {self.id} - {self.accession.accession_tag}"


class Captains(models.Model):
    """Captains table from SQLite schema"""
    class Meta:
        db_table = 'starbase_captains'
    
    captainID = models.CharField(max_length=255, unique=True)
    sequence = models.TextField()
    ship = models.ForeignKey(Ships, on_delete=models.CASCADE, related_name='captains')
    reviewed = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return self.captainID


class Taxonomy(models.Model):
    """Taxonomy table from SQLite schema"""
    class Meta:
        db_table = 'starbase_taxonomy'
    
    name = models.CharField(max_length=255, null=True, blank=True)
    taxID = models.CharField(max_length=255, null=True, blank=True)
    superkingdom = models.CharField(max_length=255, null=True, blank=True)
    clade = models.CharField(max_length=255, null=True, blank=True)
    kingdom = models.CharField(max_length=255, null=True, blank=True)
    subkingdom = models.CharField(max_length=255, null=True, blank=True)
    phylum = models.CharField(max_length=255, null=True, blank=True)
    subphylum = models.CharField(max_length=255, null=True, blank=True)
    class_field = models.CharField(max_length=255, null=True, blank=True, db_column='class')  # class is reserved
    subclass = models.CharField(max_length=255, null=True, blank=True)
    order = models.CharField(max_length=255, null=True, blank=True)
    suborder = models.CharField(max_length=255, null=True, blank=True)
    family = models.CharField(max_length=255, null=True, blank=True)
    genus = models.CharField(max_length=255, null=True, blank=True)
    species = models.CharField(max_length=255, null=True, blank=True)
    section = models.CharField(max_length=255, null=True, blank=True)
    species_group = models.CharField(max_length=255, null=True, blank=True)
    subgenus = models.CharField(max_length=255, null=True, blank=True)
    strain = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"{self.name} ({self.taxID})"


class Genome(models.Model):
    """Genome table from SQLite schema"""
    class Meta:
        db_table = 'starbase_genomes'
    
    ome = models.CharField(max_length=50, null=True, blank=True)
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='genomes', null=True)
    version = models.CharField(max_length=50, null=True, blank=True)
    genomeSource = models.CharField(max_length=50, null=True, blank=True)
    citation = models.CharField(max_length=50, null=True, blank=True)
    biosample = models.CharField(max_length=50, null=True, blank=True)
    acquisition_date = models.IntegerField(null=True, blank=True)
    assembly_accession = models.CharField(max_length=50, null=True, blank=True)

    def __str__(self):
        return f"{self.ome} v{self.version}"


class Papers(models.Model):
    """Papers table from SQLite schema"""
    class Meta:
        db_table = 'starbase_papers'
    
    Key = models.CharField(max_length=255, null=True, blank=True)
    ItemType = models.CharField(max_length=255, null=True, blank=True)
    PublicationYear = models.IntegerField(null=True, blank=True)
    Author = models.CharField(max_length=255, null=True, blank=True)
    Title = models.CharField(max_length=255, null=True, blank=True)
    PublicationTitle = models.CharField(max_length=255, null=True, blank=True)
    DOI = models.CharField(max_length=255, null=True, blank=True)
    Url = models.CharField(max_length=255, null=True, blank=True)
    AbstractNote = models.TextField(null=True, blank=True)
    Date = models.CharField(max_length=255, null=True, blank=True)
    starshipMentioned = models.CharField(max_length=255, null=True, blank=True)
    typePaper = models.CharField(max_length=255, null=True, blank=True)
    shortCitation = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"{self.Author} ({self.PublicationYear}): {self.Title}"


class FamilyNames(models.Model):
    """FamilyNames table from SQLite schema"""
    class Meta:
        db_table = 'starbase_family_names'
    
    longFamilyID = models.CharField(max_length=255, null=True, blank=True)
    oldFamilyID = models.CharField(max_length=255, null=True, blank=True)
    clade = models.IntegerField(null=True, blank=True)
    newFamilyID = models.IntegerField(null=True, blank=True)
    familyName = models.CharField(max_length=255, null=True, blank=True)
    type_element_reference = models.CharField(max_length=255, null=True, blank=True)
    notes = models.CharField(max_length=255, null=True, blank=True)
    otherFamilyID = models.CharField(max_length=255, null=True, blank=True)
    paper = models.ForeignKey(Papers, on_delete=models.CASCADE, null=True, blank=True)
    
    # Many-to-many relationship with Papers
    papers = models.ManyToManyField(Papers, related_name='family_names', blank=True)

    def __str__(self):
        return self.familyName or f"Family {self.id}"

class StarshipFeatures(models.Model):
    """StarshipFeatures table from SQLite schema - more comprehensive than current Feature model"""
    class Meta:
        db_table = 'starbase_starship_features'

    contigID = models.CharField(max_length=255, null=True, blank=True)
    starshipID = models.CharField(max_length=255, null=True, blank=True)
    captainID = models.CharField(max_length=255, null=True, blank=True)
    elementBegin = models.CharField(max_length=255, null=True, blank=True)
    elementEnd = models.CharField(max_length=255, null=True, blank=True)
    elementLength = models.CharField(max_length=255, null=True, blank=True)
    strand = models.CharField(max_length=255, null=True, blank=True)
    boundaryType = models.CharField(max_length=255, null=True, blank=True)
    emptySiteID = models.CharField(max_length=255, null=True, blank=True)
    emptyContig = models.CharField(max_length=255, null=True, blank=True)
    emptyBegin = models.CharField(max_length=255, null=True, blank=True)
    emptyEnd = models.CharField(max_length=255, null=True, blank=True)
    emptySeq = models.CharField(max_length=255, null=True, blank=True)
    upDR = models.CharField(max_length=255, null=True, blank=True)
    downDR = models.CharField(max_length=255, null=True, blank=True)
    DRedit = models.CharField(max_length=255, null=True, blank=True)
    upTIR = models.CharField(max_length=255, null=True, blank=True)
    downTIR = models.CharField(max_length=255, null=True, blank=True)
    TIRedit = models.CharField(max_length=255, null=True, blank=True)
    nestedInside = models.CharField(max_length=255, null=True, blank=True)
    containNested = models.CharField(max_length=255, null=True, blank=True)
    ship = models.ForeignKey(Ships, on_delete=models.CASCADE, null=True, blank=True)
    captain = models.ForeignKey(Captains, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return f"Feature {self.starshipID} - {self.captainID}"


class Navis(models.Model):
    """Navis table from SQLite schema"""
    class Meta:
        db_table = 'starbase_navis_names'
    
    navis_name = models.CharField(max_length=255, null=True, blank=True)
    previous_navis_name = models.CharField(max_length=255, null=True, blank=True)
    ship_family = models.ForeignKey(FamilyNames, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return self.navis_name or f"Navis {self.id}"


class Haplotype(models.Model):
    """Haplotype table from SQLite schema"""
    class Meta:
        db_table = 'starbase_haplotype_names'
    
    haplotype_name = models.CharField(max_length=255, null=True, blank=True)
    previous_haplotype_name = models.CharField(max_length=255, null=True, blank=True)
    navis = models.ForeignKey(Navis, on_delete=models.CASCADE, null=True, blank=True)
    ship_family = models.ForeignKey(FamilyNames, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return self.haplotype_name or f"Haplotype {self.id}"


class Gff(models.Model):
    """Gff table from SQLite schema"""
    class Meta:
        db_table = 'starbase_gff'
    
    accession = models.ForeignKey(Accessions, on_delete=models.CASCADE, null=True, blank=True)
    source = models.CharField(max_length=255, null=True, blank=True)
    type = models.CharField(max_length=255, null=True, blank=True)
    start = models.IntegerField(null=True, blank=True)
    end = models.IntegerField(null=True, blank=True)
    phase = models.IntegerField(null=True, blank=True)
    strand = models.CharField(max_length=255, null=True, blank=True)
    score = models.CharField(max_length=255, null=True, blank=True)
    attributes = models.CharField(max_length=255, null=True, blank=True)
    ship = models.ForeignKey(Ships, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return f"GFF {self.contigID}:{self.start}-{self.end}"

class JoinedShips(models.Model):
    """JoinedShips table from SQLite schema - comprehensive starship data"""
    class Meta:
        db_table = 'starbase_joined_ships'
    
    # Starship quality flag options
    QUALITY_FLAG_CHOICES = (
        (0, "COMPLETE"),           # All critical information present
        (1, "MISSING_BOUNDARIES"), # Missing DR/TIR boundary information
        (2, "MISSING_CAPTAIN"),    # Missing captain/transposase information  
        (3, "MISSING_CLASSIFICATION"), # Missing family/navis classification
        (4, "INCOMPLETE"),         # Multiple critical pieces missing
        (5, "MISSING_CARGO_ANNOTATIONS"), # Missing cargo annotations
    )

    id = models.IntegerField(primary_key=True)
    starshipID = models.CharField(max_length=255, null=True, blank=True)
    evidence = models.CharField(max_length=255, null=True, blank=True)
    source = models.CharField(max_length=255, null=True, blank=True)
    curated_status = models.CharField(max_length=255, null=True, blank=True)
    
    # Direct fields for simplified upload process
    sequence = models.TextField(null=True, blank=True, help_text='Starship sequence')
    species = models.CharField(max_length=255, null=True, blank=True, help_text='Species name')
    
    ship_family = models.ForeignKey('FamilyNames', on_delete=models.SET_NULL, null=True, blank=True, db_column='ship_family_id')
    taxonomy = models.ForeignKey('Taxonomy', on_delete=models.SET_NULL, null=True, blank=True, db_column='tax_id')
    ship = models.ForeignKey('Ships', on_delete=models.SET_NULL, null=True, blank=True, db_column='ship_id')
    genome = models.ForeignKey('Genome', on_delete=models.SET_NULL, null=True, blank=True, db_column='genome_id')
    captain = models.ForeignKey('Captains', on_delete=models.SET_NULL, null=True, blank=True, db_column='captain_id')
    ship_navis = models.ForeignKey('Navis', on_delete=models.SET_NULL, null=True, blank=True, db_column='ship_navis_id')
    ship_haplotype = models.ForeignKey('Haplotype', on_delete=models.SET_NULL, null=True, blank=True, db_column='ship_haplotype_id')
    
    created_at = models.CharField(max_length=255, null=True, blank=True)  # text field in starbase
    updated_at = models.CharField(max_length=255, null=True, blank=True)  # text field in starbase

    # Property to access the sequence - now uses direct field
    @property
    def starship_sequence(self):
        return self.sequence or ""

    # TODO: integrate `quality_flag` into existing tables

    def calculate_quality_flag(self):
        """
        Calculate the appropriate quality flag based on available data.
        Returns the flag value (integer).
        """
        missing_data = []
        
        # Check for captain/transposase information
        has_captain = bool(self.captain)
        if not has_captain:
            missing_data.append('captain')
        
        # Check for classification information
        has_classification = bool(
            self.ship_family or 
            (self.navis and self.haplotype)
        )
        if not has_classification:
            missing_data.append('classification')
        
        # Check for basic sequence information
        has_basic_info = bool(
            self.starshipID and 
            self.accession
        )
        if not has_basic_info:
            missing_data.append('basic_info')
        
        # Check for cargo annotations (if we have a ship with annotations)
        has_cargo_annotations = False
        if self.accession and hasattr(self.accession, 'ships'):
            for ship in self.accession.ships.all():
                if hasattr(ship, 'annotations') and ship.annotations.count() > 0:
                    has_cargo_annotations = True
                    break
        
        if not has_cargo_annotations:
            missing_data.append('cargo_annotations')

        # Determine flag based on missing data
        if not missing_data:
            return 0  # COMPLETE
        elif len(missing_data) >= 3:
            return 4  # INCOMPLETE
        elif 'captain' in missing_data:
            return 2  # MISSING_CAPTAIN
        elif 'classification' in missing_data:
            return 3  # MISSING_CLASSIFICATION
        elif 'cargo_annotations' in missing_data:
            return 5  # MISSING_CARGO_ANNOTATIONS
        else:
            return 4  # INCOMPLETE
    
    def update_quality_flag(self):
        """Update the quality flag based on current data and save."""
        self.quality_flag = self.calculate_quality_flag()
        self.save(update_fields=['quality_flag'])
    
    def get_missing_data_summary(self):
        """
        Return a human-readable summary of what data is missing.
        """
        missing_items = []
        
        # Check captain
        if not self.captain:
            missing_items.append("Captain/transposase identification")
        
        # Check classification
        if not (self.ship_family or (self.navis and self.haplotype)):
            missing_items.append("Family/classification information")
        
        # Check basic info
        if not (self.starshipID and self.accession):
            missing_items.append("Basic sequence information (ID, accession)")
        
        # Check cargo annotations
        has_cargo_annotations = False
        if self.accession and hasattr(self.accession, 'ships'):
            for ship in self.accession.ships.all():
                if hasattr(ship, 'annotations') and ship.annotations.count() > 0:
                    has_cargo_annotations = True
                    break
        
        if not has_cargo_annotations:
            missing_items.append("Cargo annotations")
        
        return missing_items

    def __str__(self):
        return f"{self.starshipID}"


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
        from django.utils import timezone
        
        try:
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
