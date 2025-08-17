from django.db import models
import re
from django.core.exceptions import ValidationError
from django.shortcuts import redirect
from django.contrib.auth.models import User
from simple_history.models import HistoricalRecords
from django.conf import settings
import django.dispatch


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
    # starships = Features.objects.filter(name=name)
    starships = Starship.objects.all()
    for starship in starships:
        if starship.starship_name == name:
            raise ValidationError("%s is already a Starship name." % name, code=3)


starship_upload_complete = django.dispatch.Signal()


# Create your models here.
class Starship(models.Model):
    class Meta:
        db_table = 'starship_starship'
    
    starship_name = models.CharField(max_length=100, unique=True)
    starship_sequence = models.TextField(max_length=15000000)
    species = models.CharField(max_length=100, unique=False)
    elementBegin = models.IntegerField(null=True, blank=True)
    elementEnd = models.IntegerField(null=True, blank=True)
    contigID = models.CharField(max_length=100, unique=False)
    starship_family = models.CharField(max_length=100, unique=False)
    starship_navis = models.CharField(max_length=100, unique=False)
    starship_haplotype = models.CharField(max_length=100, unique=False)

    def __str__(self):
        return self.starship_name


class ShipFeatures(models.Model):
    contigID = models.CharField(max_length=100)
    starshipID = models.CharField(max_length=100, unique=True)
    elementBegin = models.IntegerField(default=0)
    elementEnd = models.IntegerField(default=0)


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
        Starship,
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


# Models to match schema structure of starbase SQLite
# having the schema/tables matching will be useful when we develop mariadb -> sqlite migration
# TODO: the data displayed in tables should be created from these updated models
class Accessions(models.Model):
    """Accessions table from SQLite schema"""
    class Meta:
        db_table = 'accessions'
    
    ship_name = models.CharField(max_length=255, null=True, blank=True)
    accession_tag = models.CharField(max_length=255, unique=True)
    version_tag = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"{self.accession_tag} - {self.ship_name}"


class Ships(models.Model):
    """Ships table from SQLite schema"""
    class Meta:
        db_table = 'ships'
    
    sequence = models.TextField()
    md5 = models.CharField(max_length=32, null=True, blank=True)
    accession = models.ForeignKey(Accessions, on_delete=models.CASCADE, related_name='ships')

    def __str__(self):
        return f"Ship {self.id} - {self.accession.accession_tag}"


class Captains(models.Model):
    """Captains table from SQLite schema"""
    class Meta:
        db_table = 'captains'
    
    captainID = models.CharField(max_length=255, unique=True)
    sequence = models.TextField()
    ship = models.ForeignKey(Accessions, on_delete=models.CASCADE, related_name='captains')
    reviewed = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return self.captainID


class Taxonomy(models.Model):
    """Taxonomy table from SQLite schema"""
    class Meta:
        db_table = 'taxonomy'
    
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

    def __str__(self):
        return f"{self.name} ({self.taxID})"


class Genome(models.Model):
    """Genome table from SQLite schema"""
    class Meta:
        db_table = 'genomes'
    
    ome = models.CharField(max_length=50, null=True, blank=True)
    taxonomy = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='genomes', null=True)
    version = models.CharField(max_length=50, null=True, blank=True)
    genomeSource = models.CharField(max_length=50, null=True, blank=True)
    citation = models.CharField(max_length=50, null=True, blank=True)
    biosample = models.CharField(max_length=50, null=True, blank=True)
    acquisition_date = models.IntegerField(null=True, blank=True)

    def __str__(self):
        return f"{self.ome} v{self.version}"


class Papers(models.Model):
    """Papers table from SQLite schema"""
    class Meta:
        db_table = 'papers'
    
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
        db_table = 'family_names'
    
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
        db_table = 'starship_features'
    
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
    ship = models.ForeignKey(Accessions, on_delete=models.CASCADE, null=True, blank=True)
    captain = models.ForeignKey(Captains, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return f"Feature {self.starshipID} - {self.captainID}"


class NavisHaplotype(models.Model):
    """NavisHaplotype table from SQLite schema"""
    class Meta:
        db_table = 'navis_haplotype'
    
    navis_name = models.CharField(max_length=255, null=True, blank=True)
    haplotype_name = models.CharField(max_length=255, null=True, blank=True)
    ship_family = models.ForeignKey(FamilyNames, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return f"{self.navis_name} - {self.haplotype_name}"


class Gff(models.Model):
    """Gff table from SQLite schema"""
    class Meta:
        db_table = 'gff'
    
    contigID = models.CharField(max_length=255, null=True, blank=True)
    accession = models.CharField(max_length=255, null=True, blank=True)
    source = models.CharField(max_length=255, null=True, blank=True)
    type = models.CharField(max_length=255, null=True, blank=True)
    start = models.IntegerField(null=True, blank=True)
    end = models.IntegerField(null=True, blank=True)
    phase = models.IntegerField(null=True, blank=True)
    strand = models.CharField(max_length=255, null=True, blank=True)
    score = models.CharField(max_length=255, null=True, blank=True)
    attributes = models.CharField(max_length=255, null=True, blank=True)
    ship = models.ForeignKey(Accessions, on_delete=models.CASCADE, null=True, blank=True)

    def __str__(self):
        return f"GFF {self.contigID}:{self.start}-{self.end}"


class JoinedShips(models.Model):
    """JoinedShips table from SQLite schema - comprehensive starship data"""
    class Meta:
        db_table = 'joined_ships'
    
    starshipID = models.CharField(max_length=255, null=True, blank=True)
    genus = models.CharField(max_length=255, null=True, blank=True)
    species = models.CharField(max_length=255, null=True, blank=True)
    strain = models.CharField(max_length=255, null=True, blank=True)
    evidence = models.CharField(max_length=255, null=True, blank=True)
    source = models.CharField(max_length=255, null=True, blank=True)
    contigID = models.CharField(max_length=255, null=True, blank=True)
    captainID = models.CharField(max_length=255, null=True, blank=True)
    elementBegin = models.IntegerField(null=True, blank=True)
    elementEnd = models.IntegerField(null=True, blank=True)
    size = models.IntegerField(null=True, blank=True)
    strand = models.CharField(max_length=255, null=True, blank=True)
    boundaryType = models.CharField(max_length=255, null=True, blank=True)
    emptySiteID = models.CharField(max_length=255, null=True, blank=True)
    emptyContig = models.CharField(max_length=255, null=True, blank=True)
    emptyBegin = models.IntegerField(null=True, blank=True)
    emptyEnd = models.IntegerField(null=True, blank=True)
    emptySeq = models.CharField(max_length=255, null=True, blank=True)
    upDR = models.CharField(max_length=255, null=True, blank=True)
    downDR = models.CharField(max_length=255, null=True, blank=True)
    DRedit = models.CharField(max_length=255, null=True, blank=True)
    upTIR = models.CharField(max_length=255, null=True, blank=True)
    downTIR = models.CharField(max_length=255, null=True, blank=True)
    TIRedit = models.CharField(max_length=255, null=True, blank=True)
    nestedInside = models.CharField(max_length=255, null=True, blank=True)
    containNested = models.CharField(max_length=255, null=True, blank=True)
    dr = models.CharField(max_length=255, null=True, blank=True)
    tir = models.CharField(max_length=255, null=True, blank=True)
    navis_name = models.CharField(max_length=255, null=True, blank=True)
    haplotype_name = models.CharField(max_length=255, null=True, blank=True)
    target = models.CharField(max_length=255, null=True, blank=True)
    spok = models.CharField(max_length=255, null=True, blank=True)
    ars = models.CharField(max_length=255, null=True, blank=True)
    other = models.CharField(max_length=255, null=True, blank=True)
    hgt = models.CharField(max_length=255, null=True, blank=True)
    ship_family = models.ForeignKey(FamilyNames, on_delete=models.CASCADE, null=True, blank=True)
    curated_status = models.CharField(max_length=255, null=True, blank=True)
    taxid = models.IntegerField(null=True, blank=True)
    ship = models.ForeignKey(Accessions, on_delete=models.CASCADE, null=True, blank=True)
    genome_id = models.CharField(max_length=255, null=True, blank=True)
    ome = models.CharField(max_length=255, null=True, blank=True)
    orphan = models.CharField(max_length=255, null=True, blank=True)
    captainID_new = models.IntegerField(null=True, blank=True)

    def __str__(self):
        return f"{self.starshipID} - {self.genus} {self.species}"
