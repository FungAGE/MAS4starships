"""
Starbase models for connecting to the SQLite starbase database.

This file contains the models that connect directly to the SQLite starbase database.
"""

from django.db import models

class Accessions(models.Model):
    class Meta:
        db_table = 'accessions'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    ship_name = models.CharField(max_length=255, null=True, blank=True)
    accession_tag = models.CharField(max_length=255, unique=True)
    version_tag = models.IntegerField(null=True, blank=True)

    def __str__(self):
        return f"{self.accession_tag} - {self.ship_name}"



class Ships(models.Model):
    class Meta:
        db_table = 'ships'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    sequence = models.TextField(null=True, blank=True)
    md5 = models.CharField(max_length=32, null=True, blank=True)
    rev_comp_md5 = models.CharField(max_length=32, null=True, blank=True)
    accession_id = models.ForeignKey(Accessions, on_delete=models.CASCADE, related_name='ships', db_column='accession_id')

    def __str__(self):
        return f"Ship {self.id} - {self.accession_id.accession_tag}"



class Captains(models.Model):
    class Meta:
        db_table = 'captains'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    captainID = models.CharField(max_length=255, unique=True)
    sequence = models.TextField(null=True, blank=True)
    ship_id = models.ForeignKey(Ships, on_delete=models.CASCADE, related_name='captains', db_column='ship_id')
    reviewed = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return self.captainID


class Taxonomy(models.Model):
    class Meta:
        db_table = 'taxonomy'
        managed = False

    id = models.IntegerField(primary_key=True)    
    name = models.CharField(max_length=255, null=True, blank=True)
    taxID = models.CharField(max_length=255, null=True, blank=True)
    superkingdom = models.CharField(max_length=255, null=True, blank=True)
    clade = models.CharField(max_length=255, null=True, blank=True)
    kingdom = models.CharField(max_length=255, null=True, blank=True)
    subkingdom = models.CharField(max_length=255, null=True, blank=True)
    phylum = models.CharField(max_length=255, null=True, blank=True)
    subphylum = models.CharField(max_length=255, null=True, blank=True)
    class_name = models.CharField(max_length=255, null=True, blank=True, db_column='class')
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
    class Meta:
        db_table = 'genomes'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    ome = models.CharField(max_length=50, null=True, blank=True)
    taxonomy_id = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='genomes', null=True, blank=True, db_column='taxonomy_id')
    version = models.CharField(max_length=50, null=True, blank=True)
    genomeSource = models.CharField(max_length=50, null=True, blank=True)
    citation = models.CharField(max_length=50, null=True, blank=True)
    biosample = models.CharField(max_length=50, null=True, blank=True)
    acquisition_date = models.CharField(max_length=255, null=True, blank=True)
    assembly_accession = models.CharField(max_length=50, null=True, blank=True)

    def __str__(self):
        return f"{self.ome} v{self.version}"


class Papers(models.Model):
    class Meta:
        db_table = 'papers'
        managed = False
    
    id = models.IntegerField(primary_key=True)
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

class FamilyNames(models.Model):
    class Meta:
        db_table = 'family_names'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    longFamilyID = models.CharField(max_length=255, null=True, blank=True)
    oldFamilyID = models.CharField(max_length=255, null=True, blank=True)
    clade = models.IntegerField(null=True, blank=True)
    newFamilyID = models.IntegerField(null=True, blank=True)
    familyName = models.CharField(max_length=255, null=True, blank=True)
    type_element_reference = models.CharField(max_length=255, null=True, blank=True)
    notes = models.CharField(max_length=255, null=True, blank=True)
    otherFamilyID = models.CharField(max_length=255, null=True, blank=True)
    paper_id = models.ForeignKey(Papers, on_delete=models.CASCADE, related_name='family_names', db_column='paper_id')
    
class StarshipFeatures(models.Model):
    class Meta:
        db_table = 'starship_features'
        managed = False

    id = models.IntegerField(primary_key=True)
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
    ship_id = models.ForeignKey(Ships, on_delete=models.CASCADE, null=True, blank=True, db_column='ship_id')
    captain_id = models.ForeignKey(Captains, on_delete=models.CASCADE, null=True, blank=True, db_column='captain_id')

    def __str__(self):
        return f"Feature {self.starshipID} - {self.captainID}"


class Navis(models.Model):
    class Meta:
        db_table = 'navis_names'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    navis_name = models.CharField(max_length=255, null=True, blank=True)
    previous_navis_name = models.CharField(max_length=255, null=True, blank=True)
    ship_family_id = models.ForeignKey(FamilyNames, on_delete=models.CASCADE, related_name='navis', null=True, blank=True, db_column='ship_family_id')

    def __str__(self):
        return self.navis_name or f"Navis {self.id}"


class Haplotype(models.Model):
    class Meta:
        db_table = 'haplotype_names'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    haplotype_name = models.CharField(max_length=255, null=True, blank=True)
    previous_haplotype_name = models.CharField(max_length=255, null=True, blank=True)
    navis_id = models.ForeignKey(Navis, on_delete=models.CASCADE, related_name='haplotype', null=True, blank=True, db_column='navis_id')
    ship_family_id = models.ForeignKey(FamilyNames, on_delete=models.CASCADE, related_name='haplotype', null=True, blank=True, db_column='ship_family_id')

    def __str__(self):
        return self.haplotype_name or f"Haplotype {self.id}"


class Gff(models.Model):
    class Meta:
        db_table = 'gff'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    accession_id = models.ForeignKey(Accessions, on_delete=models.CASCADE, related_name='gff', null=True, blank=True, db_column='accession_id')
    source = models.CharField(max_length=255, null=True, blank=True)
    type = models.CharField(max_length=255, null=True, blank=True)
    start = models.IntegerField(null=True, blank=True)
    end = models.IntegerField(null=True, blank=True)
    phase = models.IntegerField(null=True, blank=True)
    strand = models.CharField(max_length=255, null=True, blank=True)
    score = models.CharField(max_length=255, null=True, blank=True)
    attributes = models.CharField(max_length=255, null=True, blank=True)
    ship_id = models.ForeignKey(Ships, on_delete=models.CASCADE, related_name='gff', null=True, blank=True, db_column='ship_id')

    def __str__(self):
        return f"GFF {self.start}-{self.end}"


class JoinedShips(models.Model):
    class Meta:
        db_table = 'joined_ships'
        managed = False
    
    id = models.IntegerField(primary_key=True)
    starshipID = models.CharField(max_length=255, null=True, blank=True)
    evidence = models.CharField(max_length=255, null=True, blank=True)
    source = models.CharField(max_length=255, null=True, blank=True)
    curated_status = models.CharField(max_length=255, null=True, blank=True)
    
    accession_id = models.ForeignKey(Accessions, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='accession_id')
    ship_id = models.ForeignKey(Ships, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='ship_id')
    ship_family_id = models.ForeignKey(FamilyNames, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='ship_family_id')
    tax_id = models.ForeignKey(Taxonomy, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='tax_id')
    genome_id = models.ForeignKey(Genome, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='genome_id')
    captain_id = models.ForeignKey(Captains, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='captain_id')
    ship_navis_id = models.ForeignKey(Navis, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='ship_navis_id')
    ship_haplotype_id = models.ForeignKey(Haplotype, on_delete=models.CASCADE, related_name='joined_ships', null=True, blank=True, db_column='ship_haplotype_id')
    
    created_at = models.CharField(max_length=255, null=True, blank=True)
    updated_at = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"{self.starshipID}"