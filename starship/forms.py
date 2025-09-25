from django import forms
from django.contrib.auth.models import User
from django.forms import CheckboxInput
from django.db.models.fields.files import FieldFile
from django.core.files.uploadedfile import TemporaryUploadedFile, InMemoryUploadedFile, File
from django.core.exceptions import ValidationError
from django.utils.safestring import mark_safe
from django.urls import reverse_lazy
from django.utils.functional import lazy
from django.utils.text import format_lazy
from django.templatetags.static import static

import re
from copy import deepcopy
import io
from argparse import Namespace

from crispy_forms.helper import FormHelper
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
import pandas as pd

from starship import models as starship_models
from starship import views as starship_views
from starship.genomic_loci_conversions import *

# TODO: add/update validation methods that already exist in starbase

def validate_fasta_file(instance):

    try:
        file_text = get_file_handle(instance, mode='r')
        #, alphabet=SingleLetterAlphabet
        #accepts fasta files and tbl files <- not good
        fasta_file = SeqIO.read(file_text, "fasta")
        if len(fasta_file) < 1:
            raise ValidationError("%s has no sequence." % instance.name)

        # Check for non-IUPAC chars
        match = re.search('[^GACTRYSWKMBDHVN]', str(fasta_file.seq).upper())
        if match:
            raise ValidationError('FASTA file sequence contains unaccepted character: {}'.format(match.group()))

    except ValueError as e:
        raise ValidationError('The following error occurred when attempting to validate fasta file: {}'.format(str(e)))


def validate_excel_file(instance):
    try:
        excel_file = get_file_handle(instance)
        file = pd.read_excel(excel_file)
    except:
        raise ValidationError("%s is not an Excel file." % instance.name)


def validate_coordinate_file(instance):
    for line in get_file_handle(instance, mode='r'):
        line_list = line.split()

        if len(line_list) != 3:
            raise ValidationError("%s incorrectly formatted. Each line should contain 3 whitespace separated values: "
                                  "strand (+ or -), start coordinate, and stop coordinate" % instance.name)

        if line_list[0] not in ['+', '-']:
            raise ValidationError("{} incorrectly formatted. First column should be '+' or '-'. Value detected: {}"
                                  "".format(instance.name, line_list[0]))

        try:
            int(line_list[1])
        except ValueError as e:
            raise ValidationError("{} incorrectly formatted. Second column should be an integer. Value detected: {}"
                                  "".format(instance.name, line_list[1]))

        try:
            int(line_list[2])
        except ValueError as e:
            raise ValidationError("{} incorrectly formatted. Third column should be an integer. Value detected: {}"
                                  "".format(instance.name, line_list[2]))


def get_file_handle(instance, mode='rb'):
    re_b = re.compile('b')
    if isinstance(instance, FieldFile):
        instance = instance.file
    if isinstance(instance, TemporaryUploadedFile):
        return open(instance.file.name, mode)
    elif isinstance(instance, InMemoryUploadedFile):
        if re_b.search(mode):
            return deepcopy(instance)
        else:
            return io.TextIOWrapper(deepcopy(instance))
    elif isinstance(instance, File):
        return open(instance.file.name, mode)
    else:
        raise ValidationError('%s may not be a file.' % instance.name)

# TODO: write/edit method to accept general feature annotations (so it does not fail if perfect CDS regions are not given)
# ! I've temporarily disabled the translation check because it was causing issues with the CDS regions
# TODO: handle upload of GFF files?
def parse_prots_from_coords(cds_fh, starship_rec, selected_table):
    for cds_line in cds_fh:
        # Get coordinates
        strand_, start_, stop_ = cds_line.split()
        start, stop, strand = coordinate_file_to_db_standard(int(start_), int(stop_), strand_)

        # produce protein sequence
        prot = get_protein_sequence(start, stop, strand, starship_rec, table=selected_table)

        # # Ensure entire sequence translated
        # if len(prot) + 1.0 != (1 + int(stop_) - int(start_)) / 3:
        #     raise TranslationError('"{line}" is not a valid CDS'.format(line=cds_line.strip()))

        yield prot, Namespace(start=start, stop=stop, strand=strand)


def get_set_of_used_speciess():
    return tuple((x, x) for x in starship_models.JoinedShips.objects.all().values_list('species', flat=True).distinct())

### START classes moved from home.forms in LIMS ###
class CrispyModelForm(forms.ModelForm):
    class Meta:
        abstract = True

    def __init__(self, *args, **kwargs):
        super(CrispyModelForm, self).__init__(*args, **kwargs)
        self.helper = CrispyHorizontalFormHelper()
        self.helper.form_tag = False
        self.helper.disable_csrf = True


class CrispyHorizontalFormHelper(FormHelper):
    def __init__(self, *args, **kwargs):
        super(CrispyHorizontalFormHelper, self).__init__(*args, **kwargs)
        self.label_class = 'col-lg-1'
        self.field_class = 'col-lg-5'
        self.form_class = 'form-horizontal'


class CrispyHorizontalInternalFormHelper(CrispyHorizontalFormHelper):
    def __init__(self, *args, **kwargs):
        super(CrispyHorizontalInternalFormHelper, self).__init__(*args, **kwargs)
        self.form_tag = False
        self.disable_csrf = True
### START classes moved from home.forms in LIMS ###

class BaseAnnotationFormset(forms.BaseInlineFormSet):

    def __int__(self, *args, **kwargs):
        super(BaseAnnotationFormset, self).__init__(*args, **kwargs)
        self.form.helper = CrispyHorizontalInternalFormHelper()
        for form in self.forms:
            form.empty_permitted = False


class Annotation_Form(CrispyModelForm):
    class Meta:
        model = starship_models.Annotation
        fields = (
            'annotation',
            'public_notes',
            'private_notes',
            'flag'
        )

    def save(self, commit=True):
        annotation = super(Annotation_Form, self).save(commit=False)
        if commit:
            annotation.save()
        return annotation


class Starship_Upload_Form(forms.Form):
    user_choices = User.objects.none()

    name = forms.CharField(
        validators=[
            starship_models.validate_starship_name,
            starship_models.validate_duplicate_name
        ],
        max_length=100,
        required=True,
    )

    upload = forms.FileField(
        validators=[validate_fasta_file],
        help_text='Must be a single fasta file!',
    )

    assign_to = forms.ModelChoiceField(
        queryset=user_choices,
        required=False,
    )

    def __init__(self, *args, **kwargs):
        super(Starship_Upload_Form, self).__init__(*args, **kwargs)
        # User choices should be limited to those with permissions to change annotations
        self.fields['assign_to'].queryset = starship_views.get_annotation_editors()


# Does not perform validation to ensure value is one of the selected options
class DynamicChoiceField(forms.ChoiceField):
    def validate(self, value):
        super(forms.ChoiceField, self).validate(value)

# displayed in upload_starship.html
class StarshipUploadForm(forms.Form):
    name = forms.CharField(
        validators=[
            starship_models.validate_starship_name,
            starship_models.validate_duplicate_name
        ],
        max_length=100,
        required=True,
    )

    upload = forms.FileField(
        validators=[validate_fasta_file],
        help_text='Must be a single FASTA file!',
    )

    annotation_file = forms.FileField(
        required=False,
        help_text='Optional: BED or GFF file containing gene annotations'
    )

    terminal_repeat = forms.IntegerField(
        required=False,
        initial=0,
        help_text='Length of terminal repeat if known'
    )

    species = forms.CharField(
        max_length=100,
        required=True,
        help_text='Species name'
    )

    assign_to = forms.ModelChoiceField(
        queryset=User.objects.none(),
        required=False,
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['assign_to'].queryset = starship_views.get_annotation_editors()

class Starship_Delete(forms.Form):
    starship = forms.ModelChoiceField(queryset=starship_models.JoinedShips.objects.all())

# A single form for manual entry of starship data across multiple models
class ComprehensiveDataForm(CrispyModelForm):
    """Single form for manual entry of starship data across multiple models"""
    
    # Accession fields
    accession_tag = forms.CharField(max_length=255, required=True)
    ship_name = forms.CharField(max_length=255, required=False)
    version_tag = forms.CharField(max_length=255, required=False)
    
    # Basic starship info (JoinedShips key fields)
    starshipID = forms.CharField(max_length=255, required=True)
    genus = forms.CharField(max_length=255, required=False)
    species = forms.CharField(max_length=255, required=False)
    strain = forms.CharField(max_length=255, required=False)
    
    # Captain info
    captainID = forms.CharField(max_length=255, required=False)
    captain_sequence = forms.CharField(widget=forms.Textarea, required=False)
    
    # Location info
    contigID = forms.CharField(max_length=255, required=False)
    elementBegin = forms.IntegerField(required=False)
    elementEnd = forms.IntegerField(required=False)
    
    # Family/Classification
    family_name = forms.CharField(max_length=255, required=False)
    navis_name = forms.CharField(max_length=255, required=False)
    haplotype_name = forms.CharField(max_length=255, required=False)
    
    # Sequence data
    ship_sequence = forms.CharField(widget=forms.Textarea, required=False)
    
    def save(self):
        """Custom save method to create records across multiple models"""
        # Create Accession
        accession = starship_models.Accessions.objects.create(
            accession_tag=self.cleaned_data['accession_tag'],
            ship_name=self.cleaned_data['ship_name'],
            version_tag=self.cleaned_data['version_tag']
        )
        
        # Create Ship if sequence provided
        if self.cleaned_data['ship_sequence']:
            ship = Ships.objects.create(
                sequence=self.cleaned_data['ship_sequence'],
                accession=accession
            )
        
        # Create Captain if provided
        if self.cleaned_data['captainID']:
            captain = Captains.objects.create(
                captainID=self.cleaned_data['captainID'],
                sequence=self.cleaned_data['captain_sequence'] or '',
                ship=accession
            )
        
        # Create JoinedShips record
        joined_ship = JoinedShips.objects.create(
            starshipID=self.cleaned_data['starshipID'],
            genus=self.cleaned_data['genus'],
            species=self.cleaned_data['species'],
            strain=self.cleaned_data['strain'],
            captainID=self.cleaned_data['captainID'],
            contigID=self.cleaned_data['contigID'],
            elementBegin=self.cleaned_data['elementBegin'],
            elementEnd=self.cleaned_data['elementEnd'],
            navis_name=self.cleaned_data['navis_name'],
            haplotype_name=self.cleaned_data['haplotype_name'],
            ship=accession
        )
        
        return joined_ship

# Bulk upload forms for SQLite data import
class BulkDataUploadForm(forms.Form):
    """Form for uploading bulk data from SQLite database"""
    sqlite_file = forms.FileField(
        help_text='SQLite database file to import data from',
        required=False
    )
    
    csv_file = forms.FileField(
        help_text='CSV file containing data to import',
        required=False
    )
    
    data_type = forms.ChoiceField(
        choices=[
            ('accessions', 'Accessions'),
            ('ships', 'Ships'),
            ('captains', 'Captains'),
            ('taxonomy', 'Taxonomy'),
            ('genomes', 'Genomes'),
            ('papers', 'Papers'),
            ('family_names', 'Family Names'),
            ('starship_features', 'Starship Features'),
            ('navis_haplotype', 'Navis Haplotype'),
            ('gff', 'GFF'),
            ('joined_ships', 'Joined Ships'),
        ],
        required=True,
        help_text='Type of data to import'
    )

    def __init__(self, *args, **kwargs):
        super(Starship_Delete, self).__init__(*args, **kwargs)
        self.helper = CrispyHorizontalFormHelper()
        self.helper.form_tag = False
        self.helper.disable_csrf = True


class Annotations_Delete(forms.Form):
    annotations = forms.ModelMultipleChoiceField(
        queryset=starship_models.Annotation.objects.all(),
        required=False
    )

    def __init__(self, *args, **kwargs):
        annotations_to_delete = kwargs.pop('annotations_to_delete', None)
        super(Annotations_Delete, self).__init__(*args, **kwargs)
        if annotations_to_delete is not None:
            self.fields['annotations'].queryset = annotations_to_delete


class Confirm_Delete(forms.ModelForm):
    confirm_delete = forms.BooleanField(
        widget=CheckboxInput(attrs={'class': 'confirm_delete'}),
        required=False,
        initial=False
    )

    class Meta:
        model = starship_models.Annotation
        fields = (
            # 'annotation',
        )


class Upload_Annotation(forms.Form):
    upload = forms.FileField(
        validators=[validate_excel_file],
        required=True,
        help_text='Must be a excel file'
    )


class AccessionForm(CrispyModelForm):
    """Form for creating/editing Accessions"""
    class Meta:
        model = starship_models.Accessions
        fields = ['ship_name', 'accession_tag', 'version_tag']


class Confirm_Upload_Annotation(forms.Form):

    def __init__(self, *args, **kwargs):
        super(Confirm_Upload_Annotation, self).__init__(*args, **kwargs)
        # self.initial['user_annotation'] = self.kwargs['user_annotation']
        # self.initial['user_public_note'] = self.kwargs['user_public_note']
        # self.initial['user_private_note'] = self.kwargs['user_private_note']
        # self.initial['user_flag'] = self.kwargs['user_flag']

    # select = forms.BooleanField(
    #
    # )
    choices = [('New', 'New'),
               ('Original', 'Original'),
               ('Custom', 'Custom')]

    select_annotation = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )
    select_public_note = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )
    select_private_note = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )
    select_flag = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )

    custom_annotation = forms.CharField(
        widget=forms.TextInput,
        required=False
    )
    custom_public_note = forms.CharField(
        widget=forms.TextInput,
        required=False
    )
    custom_private_note = forms.CharField(
        widget=forms.TextInput,
        required=False
    )
    custom_flag = forms.ChoiceField(
        choices=starship_models.Annotation.flag_options,
        required=False
    )

    user_annotation = forms.CharField(
        # hidden=True
        required=False
        # widget=forms.TextInput
    )
    user_flag = forms.ChoiceField(
        choices=starship_models.Annotation.flag_options,
        # hidden=True
        required=False
    )
    user_public_note = forms.CharField(
        # hidden=True
        required=False
    )
    user_private_note = forms.CharField(
        # hidden=True
        required=False
    )

    db_annotation = forms.CharField(
        # hidden=True
        required=False
        # widget=forms.TextInput
    )
    db_flag = forms.ChoiceField(
        choices=starship_models.Annotation.flag_options,
        # hidden=True
        required=False
    )
    db_public_note = forms.CharField(
        # hidden=True
        required=False
    )
    db_private_note = forms.CharField(
        # hidden=True
        required=False
    )
    db_pk = forms.IntegerField(
        required=True
    )

    class Meta:
        readonly = (
            'db_pk',
            'db_private_note',
            'db_public_note',
            'db_flag',
            'db_annotation',
            'user_private_note',
            'user_public_note',
            'user_flag',
            'user_annotation',
        )
        # model = starship_models.Annotation
        # fields = (
        #     # 'annotation'
        #     # 'private_notes'
        #     # 'public_notes'
        #     # 'flag'
        # )
