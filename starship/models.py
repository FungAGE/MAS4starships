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
    elementBegin = models.CharField(max_length=100, unique=False)
    elementEnd = models.CharField(max_length=100, unique=False)
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
