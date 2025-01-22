from django.http import HttpResponse
from django.views import generic
from django.urls import reverse_lazy
from django.forms import modelformset_factory, formset_factory
from django.shortcuts import render, redirect
from django.contrib.auth.models import User
from django.conf import settings
from django.dispatch import receiver
from django.db.models.signals import post_save, post_delete
from django.core.cache import cache
from django.db.models import Q
from django.contrib.auth.mixins import LoginRequiredMixin, PermissionRequiredMixin
from django.db import transaction

import os
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from tempfile import TemporaryDirectory
import re
import pandas as pd
from pandas import ExcelWriter
import argparse
import tarfile
import sys
import subprocess
import multiprocessing

from result_viewer.views import MixinForBaseTemplate, add_context_for_starship_viz

from starship import gene_calling
from starship.genomic_loci_conversions import *
from starship import create_deliverables
from starship.forms import get_file_handle
from starship import forms as starship_forms
from starship import models as starship_models
from starship.tasks import create_CDS_annotations, create_trna_annotations, add_annotations_and_features_to_db, create_custom_CDS_annotations

# Annotation history information
class Annotation_History(LoginRequiredMixin, MixinForBaseTemplate, generic.View):
    template_name = 'starship/annotation_history.html'

    def get(self, request, annotation_pk):
        context = {}
        current_annotation = starship_models.Annotation.objects.get(id=annotation_pk)
        history = current_annotation.history.all().order_by('history_date')
        context = {'current_annotation': current_annotation, 'history': history}
        return render(request, self.template_name, context)


# return users annotations based on user that made the request
class My_Annotations(LoginRequiredMixin, MixinForBaseTemplate, generic.View):
    def get(self, request):
        user = request.user
        annotations = starship_models.Annotation.objects.filter(assigned_to=user).\
            prefetch_related('feature_set', 'feature_set__starship').select_related('assigned_to')
        context = {'user': user, 'annotations': annotations}
        return render(request, 'starship/my_annotations.html', context)


# In the name, used when uploading annotations
class Upload_Annotation(LoginRequiredMixin, PermissionRequiredMixin, MixinForBaseTemplate, generic.View):
    confirm_annotation_formset_factory = formset_factory(
        form=starship_forms.Confirm_Upload_Annotation,
        extra=0,
    )

    template_name = 'starship/upload_annotation.html'
    permission_required = 'starship.add_annotation'
    permission_denied_message = 'You do not have permission to access this page. Please contact your administrator.'
    login_url = reverse_lazy('login')

    def get(self, request):
        upload_form = starship_forms.Upload_Annotation

        context = self.get_context_data()
        context['annotation_upload_form'] = upload_form
        return render(request, self.template_name, context)

    def post(self, request):
        context = {}
        upload_form = starship_forms.Upload_Annotation(request.POST, request.FILES)

        if upload_form.is_valid():
            upload_annotation_df = pd.read_excel(upload_form.cleaned_data['upload'])

            # Clean up columns names to be lowercase, no trailing spaces, and no intenal spaces
            upload_annotation_df.columns = upload_annotation_df.columns.str.lower()
            upload_annotation_df.columns = upload_annotation_df.columns.str.strip()
            upload_annotation_df.columns = upload_annotation_df.columns.str.replace(' ', '_')

            upload_annotation_df = upload_annotation_df.fillna('')

            annotation_form_data_list = []

            unannotated_upload_count = 0
            skipped_upload_count = 0
            need_review_count = 0
            failed_count = 0
            failed_list = []

            # Derived from stackoverflow.com/questions/4727327/
            flag_options_reverse = dict((v, k) for k, v in starship_models.Annotation.flag_options)

            # Correct column names for backward compatibility to upload annotations
            # Match phate output
            if 'manual_annotation' in upload_annotation_df:
                upload_annotation_df = upload_annotation_df.rename(columns={'manual_annotation': 'final_annotation'})
            # To match database and match phate output
            if 'internal_notes' in upload_annotation_df:
                upload_annotation_df = upload_annotation_df.rename(columns={'internal_notes': 'private_notes'})
            # From phate output
            if 'gene_note' in upload_annotation_df:
                upload_annotation_df = upload_annotation_df.rename(columns={'gene_note': 'public_notes'})

            # Confirm all columns are present
            if 'final_annotation' not in upload_annotation_df \
                    or 'public_notes' not in upload_annotation_df \
                    or 'private_notes' not in upload_annotation_df \
                    or 'flag' not in upload_annotation_df \
                    or 'protein_sequence' not in upload_annotation_df:
                upload_form.errors.update(Names=': Please change the headers of your excel file to "Final_Annotation",'
                                                ' "Public_Notes", "Private_Notes", "Flag", and "Protein_Sequence"')
                return render(request, self.template_name, {'annotation_upload_form': upload_form})

            for index, row in upload_annotation_df.iterrows():
                # check if the sequence provided is in the database
                flag = row.flag
                flag = flag.upper()
                flag = flag.strip()
                if flag == 'TRNA':
                    flag = 'tRNA'
                flag = flag.replace('_', ' ')
                if flag not in flag_options_reverse:
                    # upload_form.errors.update(': The flag %s is not a valid flag.' % flag)
                    failed_count += 1
                    failed_list.append(
                        {
                            'annotation': row.final_annotation,
                            'public_notes': row.public_notes,
                            'private_notes': row.private_notes,
                            'flag': row.flag,
                            'sequence': row.protein_sequence,
                        }
                    )
                    continue
                    # return render(request, self.template_name, {'annotation_upload_form': upload_form})

                if starship_models.Annotation.objects.filter(sequence=row.protein_sequence):
                    db_annotation = starship_models.Annotation.objects.get(sequence=row.protein_sequence)
                    if db_annotation.public_notes == 'nan':
                        db_annotation.public_notes = ''
                        db_annotation.save()
                    if db_annotation.private_notes == 'nan':
                        db_annotation.private_notes = ''
                        db_annotation.save()
                else:
                    upload_form.errors.update(
                        Sequence=': The annotation %s with Protein Sequence %s is not in the database' % (
                            row.final_annotation,
                            row.protein_sequence
                        )
                    )
                    failed_count += 1
                    failed_list.append(
                        {
                            'annotation': row.final_annotation,
                            'public_notes': row.public_notes,
                            'private_notes': row.private_notes,
                            'flag': row.flag,
                            'sequence': row.protein_sequence,
                        }
                    )
                    continue
                    # return render(request, self.template_name, {'annotation_upload_form': upload_form})

                if row.final_annotation == db_annotation.annotation \
                        and row.public_notes == db_annotation.public_notes \
                        and row.private_notes == db_annotation.private_notes \
                        and flag == db_annotation.get_flag_display():
                    skipped_upload_count += 1
                elif db_annotation.get_flag_display() == 'UNANNOTATED':
                    # automatically save the user annotations to the unannotated database annotation
                    annotation = starship_models.Annotation.objects.get(sequence=row.protein_sequence)
                    annotation.annotation = row.final_annotation
                    annotation.public_notes = row.public_notes
                    annotation.private_notes = row.private_notes
                    annotation.flag = flag_options_reverse[flag]
                    annotation.save()
                    unannotated_upload_count += 1

                else:
                    annotation_form_data = {}

                    annotation_form_data['user_annotation'] = row.final_annotation
                    annotation_form_data['user_flag'] = flag_options_reverse[flag]
                    annotation_form_data['user_public_note'] = row.public_notes
                    annotation_form_data['user_private_note'] = row.private_notes

                    annotation_form_data['db_annotation'] = db_annotation.annotation
                    annotation_form_data['db_flag'] = db_annotation.flag
                    annotation_form_data['db_public_note'] = db_annotation.public_notes
                    annotation_form_data['db_private_note'] = db_annotation.private_notes
                    annotation_form_data['db_pk'] = db_annotation.pk

                    annotation_form_data['select_annotation'] = 'New'
                    annotation_form_data['select_flag'] = 'New'
                    annotation_form_data['select_public_note'] = 'New'
                    annotation_form_data['select_private_note'] = 'New'

                    annotation_form_data_list.append(annotation_form_data)
                    need_review_count += 1

            confirm_upload_annotation_formset = self.confirm_annotation_formset_factory(
                initial=annotation_form_data_list
            )

            context['confirm_upload_annotation_formset'] = confirm_upload_annotation_formset
            context['skipped_upload_count'] = skipped_upload_count
            context['unannotated_upload_count'] = unannotated_upload_count
            context['need_review_count'] = need_review_count
            context['failed_count'] = failed_count
            context['failed_list'] = failed_list

            return render(request, 'starship/confirm_upload_annotation.html', context)
        else:
            error_context = self.get_context_data()
            error_context['annotation_upload_form'] = upload_form
            return render(request, self.template_name, error_context)


# used when user needs to check with annotations they want to update
class Confirm_Upload_Annotation(LoginRequiredMixin, PermissionRequiredMixin, MixinForBaseTemplate, generic.View):
    confirm_upload_annotation_formset_factory = formset_factory(
        form=starship_forms.Confirm_Upload_Annotation,
        extra=0,
    )
    permission_required = ('starship.add_annotation', 'starship.change_annotation')
    permission_denied_message = 'You do not have permission to access this page. Please contact your administrator.'
    login_url = reverse_lazy('login')

    def get(self, request):
        context = self.get_context_data()
        upload_form = starship_forms.Upload_Annotation(
            data=request.GET
        )
        # context['starship_form'] = deepcopy(starship_form)
        if not upload_form.is_valid():
            context['upload_form'] = upload_form
            return render(request, 'starship/upload_annotation.html', context)

        return render(request, 'starship/confirm_upload_annotation.html', context)

    def post(self, request):
        confirm_upload_annotation_formset = self.confirm_upload_annotation_formset_factory(data=request.POST)
        if not confirm_upload_annotation_formset.is_valid():
            return redirect('starship:upload_annotations')

        with transaction.atomic():
            for upload_annotation_form in confirm_upload_annotation_formset:
                db_pk = upload_annotation_form.cleaned_data['db_pk']
                annotation = starship_models.Annotation.objects.get(pk=db_pk)

                if upload_annotation_form.cleaned_data['select_annotation'] == 'New':
                    annotation.annotation = upload_annotation_form.cleaned_data['user_annotation']
                elif upload_annotation_form.cleaned_data['select_annotation'] == 'Custom':
                    annotation.annotation = upload_annotation_form.cleaned_data['custom_annotation']

                if upload_annotation_form.cleaned_data['select_private_note'] == 'New':
                    annotation.private_notes = upload_annotation_form.cleaned_data['user_private_note']
                elif upload_annotation_form.cleaned_data['select_private_note'] == 'Custom':
                    annotation.private_notes = upload_annotation_form.cleaned_data['custom_private_note']

                if upload_annotation_form.cleaned_data['select_public_note'] == 'New':
                    annotation.public_notes = upload_annotation_form.cleaned_data['user_public_note']
                elif upload_annotation_form.cleaned_data['select_public_note'] == 'Custom':
                    annotation.public_note = upload_annotation_form.cleaned_data['custom_public_note']

                if upload_annotation_form.cleaned_data['select_flag'] == 'New':
                    annotation.flag = upload_annotation_form.cleaned_data['user_flag']
                elif upload_annotation_form.cleaned_data['select_flag'] == 'Custom':
                    annotation.flag = upload_annotation_form.cleaned_data['custom_flag']

                annotation.save()

        return redirect('starship:starship_list')


class Upload_Starship(LoginRequiredMixin, PermissionRequiredMixin, MixinForBaseTemplate, generic.View):
    model = starship_models.Starship
    permission_required = 'starship.add_starship'
    permission_denied_message = 'You do not have permission to access this page. Please contact your administrator.'
    login_url = reverse_lazy('login')


class Upload_Starship(Upload_Starship):
    template_name = 'starship/upload_starship.html'
    # context_object_name = 'starship'

    def get(self, request):
        upload_form = starship_forms.Starship_Upload_Form()

        context = self.get_context_data()
        context['upload_form'] = upload_form

        return render(request, self.template_name, context)

    def post(self, request):
        new_annotations = {}  # Annotation objects keyed by sequence
        new_features = []
        upload_form = starship_forms.Starship_Upload_Form(request.POST, request.FILES)

        if upload_form.is_valid():

            # Test if blast database exists
            for p in ["{db}.{ext}".format(db=settings.TERMINASE_DATABASE, ext=x) for x in
                      ['phr', 'pin', 'psq']]:
                if not os.path.isfile(p):
                    upload_form.errors.update(
                        circularly_permutated=': Terminase blast database does not exist. Looked for {}'.format(p)
                    )
                    return render(request, self.template_name, {'upload_form': upload_form})

            with TemporaryDirectory() as tempdir:
                output_dest = tempdir

                file = get_file_handle(request.FILES['upload'], mode='r')
                genome = SeqIO.read(file, 'fasta').seq.__str__().upper()
                name = request.POST['name']

                terminal_repeat = request.POST['terminal_repeat']

                starship_path = request.FILES['upload'].temporary_file_path()
                syspath = sys.path
                print(syspath)
                starship_cds = gene_calling.run_glimmer(starship_path, name, output_dest)

                starship = None
                feature_saved_list = []
                protein_list = []
                repeat = int(terminal_repeat)

                try:
                    with transaction.atomic():
                        starship = starship_models.Starship(starship_name=name, starship_sequence=genome, species='starship')
                        # needs to be saved before features, so you can access the pk
                        starship.save()

                        # if the starship is circularly permuted then get all the proteins from the starship and place them in temp file.
                        if upload_form.cleaned_data['circularly_permuted'] is True:
                            if starship_cds:
                                for cds in gene_calling.parse_glimmer_results(starship_cds):
                                    sequence = Seq(genome, IUPAC.ambiguous_dna)
                                    protein = get_protein_sequence(cds.start, cds.stop, cds.strand, sequence)
                                    record = SeqRecord(Seq(protein._data, generic_protein), id="%s | " % starship.starship_name,
                                                       description="%s | %s | %s | CDS" % (cds.start, cds.stop, cds.strand))
                                    protein_list.append(record)
                                file_path = os.path.join(tempdir, "proteins.fasta")
                                SeqIO.write(protein_list, file_path, "fasta")

                                # Test if blast database exists
                                for p in ["{db}.{ext}".format(db=settings.TERMINASE_DATABASE, ext=x) for x in ['phr', 'pin', 'psq']]:
                                    if not os.path.isfile(p):
                                        upload_form.errors.update(
                                            circularly_permutated=': Terminase blast database does not exist. Looked for {}'.format(p)
                                        )
                                        starship.delete()
                                        return render(request, self.template_name, {'upload_form': upload_form})

                                # Run blast search for terminase
                                blast_results_path = os.path.join(tempdir, '{}.blastp_results.xml'.format(starship.starship_name))
                                blast_result = subprocess.run([
                                    'blastp',
                                    '-query', file_path,
                                    '-db', settings.TERMINASE_DATABASE,
                                    '-evalue', '0.00001',
                                    '-outfmt', '5',
                                    '-out', blast_results_path,
                                    '-num_threads', str(multiprocessing.cpu_count())
                                ], stderr=subprocess.PIPE)

                                # If blast search fails return error
                                if blast_result.returncode != 0:
                                    upload_form.errors.update(
                                        circularly_permutated=': Blastp search failed. Return code = {}\nstderr = {}'.format(
                                            blast_result.returncode,
                                            blast_result.stderr.decode()
                                        )
                                    )
                                    starship.delete()
                                    return render(request, self.template_name, {'upload_form': upload_form})

                                passing_record = []
                                for rec_i, record in enumerate(NCBIXML.parse(open(blast_results_path, 'r'))):
                                    if len(record.alignments) > 3:
                                        passing_record.append((rec_i, record))
                                length = len(passing_record)

                                if length == 0:
                                    upload_form.errors.update(circularly_permutated=': No proteins were found to have a significan'
                                                                                    't terminase signal. Please resolve manually. '
                                                                                    'Starship was not uploaded.')
                                    starship.delete()
                                    return render(request, self.template_name, {'upload_form': upload_form})

                                if length > 2:
                                    upload_form.errors.update(circularly_permutated=': Number of called proteins with significant '
                                                                                    'terminase signal is greater than 2. Please'
                                                                                    ' resolve this manually. Starship was not uploaded')
                                    starship.delete()
                                    return render(request, self.template_name, {'upload_form': upload_form})

                                if length == 1:
                                    new_fasta_path = self.make_starship_start_from_protein(passing_record[0], starship_path, tempdir, starship)

                                else:
                                    strand1 = passing_record[0][1].query.split('|')[3].strip()
                                    strand2 = passing_record[1][1].query.split('|')[3].strip()

                                    if strand1 != strand2:
                                        upload_form.errors.update(circularly_permutated=': Small and large terminase on different'
                                                                                        'strands. Resolve manually. Starship not uploaded.')
                                        starship.delete()
                                        return render(request, self.template_name, {'upload_form': upload_form})

                                    if strand1 == '+':
                                        new_fasta_path = self.make_starship_start_from_protein(passing_record[0], starship_path, tempdir,
                                                                            starship)
                                    else:
                                        new_fasta_path = self.make_starship_start_from_protein(passing_record[1], starship_path, tempdir,
                                                                            starship)

                                # created so the user does not have to overwrite a file
                                new_output_dest = tempdir

                                new_starship_cds = gene_calling.run_glimmer(new_fasta_path, name, new_output_dest)
                                new_starship_t_rna = gene_calling.run_trnascan_se(new_fasta_path, name, new_output_dest)

                                if new_starship_cds:
                                    create_CDS_annotations(
                                        new_starship_cds,
                                        starship,
                                        upload_form.cleaned_data['assign_to'],
                                        new_annotations,
                                        new_features
                                    )

                                if new_starship_t_rna:
                                    create_trna_annotations(
                                        new_starship_t_rna,
                                        starship,
                                        upload_form.cleaned_data['assign_to'],
                                        new_annotations,
                                        new_features
                                    )

                        else:
                            # cds features
                            if starship_cds:
                                create_CDS_annotations(
                                    starship_cds, starship, upload_form.cleaned_data['assign_to'], new_annotations, new_features
                                )

                            # tRNA features and annotations
                            starship_t_rna = gene_calling.run_trnascan_se(starship_path, name, output_dest)
                            if starship_t_rna:
                                create_trna_annotations(
                                    starship_t_rna, starship, upload_form.cleaned_data['assign_to'], new_annotations, new_features
                                )

                            # terminal repeat features
                            if repeat > 0:
                                repeat_seq = genome[:repeat]

                                if starship_models.Annotation.objects.filter(sequence=repeat_seq).count() > 0:
                                    repeat_annotation = starship_models.Annotation.objects.get(sequence=repeat_seq)
                                else:
                                    repeat_annotation = starship_models.Annotation()
                                    repeat_annotation.sequence = repeat_seq
                                    repeat_annotation.annotation = 'None'
                                    repeat_annotation.public_notes = 'Direct terminal repeat. Detected with sequencing data ' \
                                                                     'via coverage based methods.'
                                    repeat_annotation.private_notes = 'This annotation was automatically generated.'
                                    repeat_annotation.flag = 9
                                    repeat_annotation.assigned_to = None
                                    new_annotations[repeat_annotation.sequence] = repeat_annotation

                                first_feature_repeat = starship_models.Feature(genome=starship, start=0, stop=repeat,
                                                                             type='Repeat Region', strand='+',
                                                                             annotation=repeat_annotation)
                                last_feature_repeat = starship_models.Feature(genome=starship, start=len(genome) - repeat,
                                                                            stop=len(genome), type='Repeat Region',
                                                                            strand='+',
                                                                            annotation=repeat_annotation)
                                new_features.append(first_feature_repeat)
                                new_features.append(last_feature_repeat)

                        add_annotations_and_features_to_db(new_annotations, new_features)
                        starship_models.starship_upload_complete.send(sender=None)

                except Exception as e:
                    upload_form.errors.update(Conflict='Exception occured upon upload: {}'.format(e))
                    context = self.get_context_data()
                    context['upload_form'] = upload_form

                    return render(request, self.template_name, context)

            return redirect('starship:starship_list')
        else:
            context = self.get_context_data()
            context['upload_form'] = upload_form

            return render(request, self.template_name, context)

    def make_starship_start_from_protein(self, record, fasta_path, output_dest, starship):
        """
        :param record:
        :param fasta_path:  file path to your temp file with the fasta info
        :param output_dest
        :param starship_name
        :return:
        """
        starship_rec = SeqIO.read(open(fasta_path), 'fasta')

        # Find location to start from
        loci_info = record[1].query.split('|')[1].strip()
        strand = record[1].query.split('|')[3].strip()
        start = int(loci_info)

        new_seq = starship_rec.seq[start:] + starship_rec.seq[:start]

        if strand == '-':
            new_seq = new_seq.reverse_complement()

        # Save fasta with new sequence
        new_fasta_path = os.path.join(output_dest, '{}.fasta'.format(starship.starship_name))
        starship_rec.seq = new_seq
        SeqIO.write(starship_rec, new_fasta_path, 'fasta')

        starship.starship_sequence = str(new_seq)
        starship.save()

        return new_fasta_path


class Upload_Custom_Starship(Upload_Starship):
    template_name = 'starship/upload_custom_starship.html'

    def get(self, request):
        upload_form = starship_forms.Custom_Starship_Upload_Form()

        context = self.get_context_data()
        context['upload_form'] = upload_form

        return render(request, self.template_name, context)

    def post(self, request):
        upload_form = starship_forms.Custom_Starship_Upload_Form(request.POST, request.FILES)

        if upload_form.is_valid():
            new_annotations = {}  # Annotation objects keyed by sequence
            new_features = []
            starship_rec = SeqIO.read(get_file_handle(upload_form.cleaned_data['starship_upload'], mode='r'), 'fasta')
            starship_sequence = starship_rec.seq.__str__().upper()
            cds_fh = get_file_handle(upload_form.cleaned_data['cds_upload'], mode='r')
            name = upload_form.cleaned_data['name']
            assign_to = upload_form.cleaned_data['assign_to']
            species = upload_form.cleaned_data['species']
            translation_table = upload_form.cleaned_data['translation_table']

            with transaction.atomic():
                with TemporaryDirectory() as tempdir:
                    # Create Starship object
                    starship_obj = starship_models.Starship(
                        starship_name=name,
                        starship_sequence=starship_sequence,
                        species=species
                    )
                    starship_obj.save()

                    # Create genome fasta in tempdir
                    record = SeqRecord(
                        Seq(starship_sequence),
                        id=name,
                        name=name,
                        description="",
                    )
                    fasta_path = os.path.join(tempdir, '{}.fasta'.format(name))
                    SeqIO.write(record, fasta_path, 'fasta')

                    # Run trnascan-se to find tRNAs if option selected
                    if upload_form.cleaned_data['run_trnascan']:
                        trnascan_output = gene_calling.run_trnascan_se(fasta_path, name, tempdir)
                        if trnascan_output:
                            create_trna_annotations(
                                trnascan_output, starship_sequence, assign_to, new_annotations, new_features
                            )

                    # Parse coordinate file and create annotation + feature objects for each CDS
                    create_custom_CDS_annotations(
                        cds_fh, translation_table, starship_obj, assign_to, new_annotations, new_features
                    )

                    add_annotations_and_features_to_db(new_annotations, new_features)

            return redirect('starship:starship_list')
        else:
            context = self.get_context_data()
            context['upload_form'] = upload_form

            return render(request, self.template_name, context)


# returns the list of genomes
class Starship_List(LoginRequiredMixin, MixinForBaseTemplate, generic.ListView):
    model = starship_models.Starship
    context_object_name = 'genomes'
    template_name = 'starship/starship_list.html'

    def get_context_data(self,  **kwargs):
        context = super(Starship_List, self).get_context_data(**kwargs)
        genomes = starship_models.Starship.objects.all().prefetch_related(
            'feature_set__annotation')

        # calculate number of unpolished CDS in starship
        context['starship_info'] = get_starship_data_dicts(genomes)

        return context


class Starship_List_SS(LoginRequiredMixin, MixinForBaseTemplate, generic.TemplateView):
    template_name = 'starship/starship_list_serverside.html'


# used for starship detail page
# @cache_page(60 * 15)
class Starship_Detail(LoginRequiredMixin, MixinForBaseTemplate, generic.DetailView):
    model = starship_models.Starship
    context_object_name = 'starship'
    template_name = 'starship/starship_detail.html'
    starship_dict = {}
    features_dict = {}

    def get_context_data(self, **kwargs):
        #get default context data
        context = super(Starship_Detail, self).get_context_data(**kwargs)

        context = add_context_for_starship_viz(context, context['starship'])

        context['features'] = starship_models.Feature.objects.filter(
            genome=context['starship']
        ).prefetch_related('annotation', 'starship')

        context['annotations'] = starship_models.Annotation.objects.filter(
            feature__in=context['features']
        ).prefetch_related('feature_set__starship', 'feature_set').select_related('assigned_to')

        self.starship_dict['starship_name'] = context['starship'].starship_name
        self.starship_dict['starship'] = context['starship'].starship_sequence
        context['starship_data'] = self.starship_dict

        upload_form = starship_forms.Starship_Upload_Form
        context['upload_form'] = upload_form

        # try:  TODO: Implement in BDRD version only
        #     sample_source_name_re = re.compile(r'AMD_(\w_.*_Phi_\d+).*$')
        #     parsed_starship_name = sample_source_name_re.match(self.starship_dict['starship_name']).group(1)
        #     context['sample_sources'] = sample_models.Sample_Source.objects.filter(name=parsed_starship_name)
        # except:
        #     context['sample_sources'] = None

        starship = starship_models.Starship.objects.filter(starship_name__exact=self.starship_dict['starship_name'])
        context['starship_info'] = get_starship_data_dicts(starship)

        return context


# used for feature list page
class Feature_List(LoginRequiredMixin, MixinForBaseTemplate, generic.ListView):
    model = starship_models.Feature
    context_object_name = 'starship'
    template_name = 'starship/feature_list.html'
    # test

    def get_context_data(self, **kwargs):
        context = super(Feature_List, self).get_context_data(**kwargs)
        context['features'] = starship_models.Feature.objects.all().prefetch_related('starship', 'annotation')\
            .order_by('start')
        return context


# used for feature detail page
class Feature_Detail(LoginRequiredMixin, generic.DetailView):
    model = starship_models.Feature
    context_object_name = 'feature'
    template_name = 'starship/feature_detail.html'


# used for annotation list page
class Annotation_List(LoginRequiredMixin, MixinForBaseTemplate, generic.ListView):
    confirm_annotation_formset_factory = modelformset_factory(
        model=starship_models.Annotation,
        form=starship_forms.Confirm_Delete,
        extra=0,
    )
    model = starship_models.Annotation
    context_object_name = 'annotations'
    template_name = 'starship/annotation_list.html'

    def get_context_data(self, **kwargs):
        context = super(Annotation_List, self).get_context_data(**kwargs)
        context['upload_form'] = starship_forms.Starship_Upload_Form
        if 'flag' in self.kwargs:
            flag = str(self.kwargs['flag'])
            if flag != 'tRNA':
                flag = flag.upper()
            flag = flag.replace('_', ' ')
            context['flag'] = flag
            # Derived from stackoverflow.com/questions/4727327/
            flag_options_reverse = dict((v, k) for k, v in starship_models.Annotation.flag_options)
            if flag in flag_options_reverse:
                annotations = starship_models.Annotation.objects.filter(
                    flag=flag_options_reverse[flag]
                ).prefetch_related(
                    'feature_set',
                    'feature_set__starship'
                ).select_related('assigned_to')
            else:
                annotations = starship_models.Annotation.objects.none()
        else:
            annotations = starship_models.Annotation.objects.all().prefetch_related(
                'feature_set',
                'feature_set__starship'
            ).select_related('assigned_to')
        maximum_display_limit = 5000
        if annotations.count() > maximum_display_limit:
            context['large_annotation_count'] = annotations.count()
            context['maximum_display_limit'] = maximum_display_limit
            annotations = annotations[:maximum_display_limit]
        context['annotations'] = annotations
        return context


class Annotation_List_Serverside(LoginRequiredMixin, MixinForBaseTemplate, generic.TemplateView):
    template_name = 'starship/annotation_list_serverside.html'


# used for annotation detail page
class Annotation_Detail(LoginRequiredMixin, MixinForBaseTemplate, generic.DetailView):
    model = starship_models.Annotation
    context_object_name = 'annotation'
    template_name = 'starship/annotation_detail.html'

    def get_context_data(self, **kwargs):
        context = super(Annotation_Detail, self).get_context_data(**kwargs)
        genomes = starship_models.Starship.objects.none()
        for feature in context['annotation'].feature_set.all():
            genomes = genomes | starship_models.Starship.objects.filter(id=feature.starship_id)
        context['starship_info'] = get_starship_data_dicts(genomes)

        context['exact_names'] = starship_models.Annotation.objects.filter(
            annotation__iexact=context['annotation'].annotation
        )
        context['exact_names'] = context['exact_names'].exclude(id=context['annotation'].id)

        return context


# used to show the nucleotide sequence through starship list view
class Get_Starship(LoginRequiredMixin, generic.View):
    def get(self, request):
        starship_id = request.GET.get('starship_id', None)
        context = {}
        if starship_id:
            context['starship'] = starship_models.Starship.objects.get(pk=starship_id)
            # context['features'] = starship.feature_set.all()
            print(starship_id)
        return render(request, 'starship/starship_sequence.html', context)


# used to show amino acid sequence through annotation list view
class Get_AA_Sequence(LoginRequiredMixin, generic.View):
    def get(self, request):
        annotation_id = request.GET.get('annotation_id', None)
        context = {}
        if annotation_id:
            context['annotation'] = starship_models.Annotation.objects.get(pk=annotation_id)
            context['header'] = ">%s | %s | %s | %s | %s" % (
                context['annotation'].accession,
                context['annotation'].annotation,
                context['annotation'].public_notes,
                context['annotation'].private_notes,
                context['annotation'].get_flag_display()
            )
        return render(request, 'starship/aa_sequence.html', context)


# used to show feature sequence through feature list view
class Get_Feature_Sequence(LoginRequiredMixin, generic.View):
    def get(self, request):
        feature_id = request.GET.get('feature_id', None)
        context = {}
        if feature_id:
            feature = starship_models.Feature.objects.get(pk=feature_id)

            context['feature'] = feature
            context['feature_sequence'] = get_dna_sequence(
                feature.start, feature.stop, feature.strand, Seq(feature.starship.starship_sequence)
            )

        return render(request, 'starship/feature_sequence.html', context)


# loads the delete starship page
class Starship_Delete(LoginRequiredMixin, PermissionRequiredMixin, MixinForBaseTemplate, generic.ListView):
    model = starship_models.Starship
    context_object_name = 'starship'
    template_name = 'starship/starship_delete.html'
    permission_required = 'starship.delete_starship'
    permission_denied_message = 'You do not have permission to access this page. Please contact your administrator.'
    login_url = reverse_lazy('login')

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super(Starship_Delete, self).get_context_data()
        starships = context['starship']
        starship_form = starship_forms.Starship_Delete(

        )
        context['starship_form'] = starship_form
        return context


# used to confirm objects to be deleted and delets objects
class Confirm_Starship_Delete(LoginRequiredMixin, PermissionRequiredMixin, MixinForBaseTemplate, generic.View):
    confirm_annotation_formset_factory = modelformset_factory(
        model=starship_models.Annotation,
        form=starship_forms.Confirm_Delete,
        extra=0
    )
    permission_required = 'starship.delete_starship'
    permission_denied_message = 'You do not have permission to access this page. Please contact your administrator.'
    login_url = reverse_lazy('login')

    def get(self, request):
        context = {}
        starship_form = starship_forms.Starship_Delete(
            data=request.GET
        )
        context['starship_form'] = deepcopy(starship_form)
        if not starship_form.is_valid():
            context['starship_form'] = starship_form
            return render(request, 'starship/starship_delete.html', context)
        # context = {}

        genomes = starship_form.cleaned_data['starship']
        features = starship_models.Feature.objects.filter(starship__in=genomes)
        annotations = starship_models.Annotation.objects.filter(feature__in=features)

        starships_not_being_deleted = starship_models.Starship.objects.all().exclude(pk__in=genomes)
        annotations_to_keep = annotations.filter(feature__starship__in=list(starships_not_being_deleted)).distinct()

        annotations = annotations.exclude(pk__in=annotations_to_keep) # difference(annotations_to_keep) <- can't do this because can't filter resulting QS
        annotations_to_delete = annotations.filter(flag__in=[7, 8, 9])  # 7 = unannotated
        annotations_to_check = annotations.exclude(pk__in=annotations_to_delete)

        annotations_to_check_form = self.confirm_annotation_formset_factory(
            queryset=annotations_to_check
        )

        annotations_to_delete_form = starship_forms.Annotations_Delete(
            data={'annotations': annotations_to_delete},
            annotations_to_delete=annotations_to_delete
        )

        total_features = 0
        for genome in genomes:
            amount = starship.feature_set.all()
            total_features = total_features + amount.count()

        # display info of items being deleted
        context['total_features'] = total_features
        context['total_starships'] = genomes.count()
        context['anno_auto_delete'] = annotations_to_delete.count()

        context['annotations_to_check_form'] = annotations_to_check_form
        context['annotations_to_delete_form'] = annotations_to_delete_form

        return render(request, 'starship/confirm_delete.html', context)

    def post(self, request):

        starship_form = starship_forms.Starship_Delete(data=request.POST)
        if not starship_form.is_valid():
            return redirect('starship:starship_delete')

        annotations_to_delete_form = starship_forms.Annotations_Delete(
            data=request.POST,
            empty_permitted=True,
            use_required_attribute=False
        )
        if not annotations_to_delete_form.is_valid():
            return redirect('starship:starship_delete')

        annotations_to_confirm_form = self.confirm_annotation_formset_factory(data=request.POST)
        if not annotations_to_confirm_form.is_valid():
            return redirect('starship:starship_delete')

        genomes = starship_form.cleaned_data['starship']

        if 'annotations' in annotations_to_delete_form.cleaned_data:
            annotations_to_delete = annotations_to_delete_form.cleaned_data['annotations']
        else:
            annotations_to_delete = []

        annotations_to_confirm = annotations_to_confirm_form.cleaned_data

        with transaction.atomic():
            for genome in genomes:
                starship.delete()

            # delete annotations that were not checked
            for annotation in annotations_to_delete:
                if annotation.feature_set.all().count() == 0:
                    annotation.delete()

            for starship_selected in annotations_to_confirm:
                if starship_selected.get('confirm_delete') is True:
                    annotation = starship_selected.get('id')
                    annotation.delete()

        return redirect('starship:starship_list')


# Deals with the annotation list table buttons. Depending on the button selected an action will happen
class Annotation_Bulk(LoginRequiredMixin, generic.View):
    template_name = 'starship/starship_list.html'

    def post(self, request):
        context = {}

        if 'download' in request.POST:
            annotation_id_list = []
            annotation_re = re.compile('annotation-\d+')

            for key, value in request.POST.items():
                if annotation_re.match(key):
                    annotation_id_list.append(int(value))
            context = {}
            annotation_list = []
            context['annotations'] = starship_models.Annotation.objects.filter(pk__in=annotation_id_list)

            for annotation in context['annotations']:
                anno = annotation.annotation
                aa_sequence = annotation.sequence
                public_note = annotation.public_notes
                private_note = annotation.private_notes
                flag = annotation.get_flag_display()

                sequence = SeqRecord(Seq(aa_sequence, generic_protein), id=annotation.accession + " |",
                                     description="%s | %s | %s | %s" % (anno, public_note, private_note, flag))
                annotation_list.insert(0, sequence)

            with TemporaryDirectory() as temp:
                file_path = os.path.join(temp, "Annotations.faa")
                SeqIO.write(annotation_list, file_path, "fasta")

                file = open(file_path, 'rb')
                response = HttpResponse(file, content_type='text/fasta')

                response['Content-Disposition'] = 'attachment; filename=Annotations.faa'

            return response
        if 'assign_to_user' in request.POST:
            # starship = None
            annotation_id_list = []
            annotation_re = re.compile('annotation-\d+')

            assigned_user = request.POST['assign_to']
            # assigned_user = assigned_user.split()

            for key, value in request.POST.items():
                if annotation_re.match(key):
                    annotation_id_list.append(int(value))
            # if length of annotation_id_list is 0 then return page
            if len(annotation_id_list) == 0:
                return redirect(request.META['HTTP_REFERER'])

            # dummy_annotation = starship_models.Annotation.objects.get(id=annotation_id_list[0])
            # for feature in dummy_annotation.feature_set.all():
            #     starship = feature.starship

            context['annotations'] = starship_models.Annotation.objects.filter(pk__in=annotation_id_list)
            users = User.objects.all()
            user = None
            for us in users:
                if assigned_user == '':
                    continue
                for group in us.groups.all():
                    if us.pk == int(assigned_user):
                        user = us
            for annotation in context['annotations']:
                annotation.assigned_to = user
                annotation.save()

            return redirect(request.META['HTTP_REFERER'])

        if 'excel' in request.POST:
            annotation_id_list = []
            annotation_re = re.compile('annotation-\d+')

            for key, value in request.POST.items():
                if annotation_re.match(key):
                    annotation_id_list.append(int(value))

            context['annotations'] = starship_models.Annotation.objects.filter(pk__in=annotation_id_list)
            with TemporaryDirectory() as temp:
                file_path = os.path.join(temp, 'Annotations_excel.xlsx')

                final_annotation = []
                public_notes = []
                private_notes = []
                flag = []
                protein_sequence = []
                for annotation in context['annotations']:
                    final_annotation.append(annotation.annotation)
                    public_notes.append(annotation.public_notes)
                    private_notes.append(annotation.private_notes)
                    flag.append(annotation.get_flag_display())
                    protein_sequence.append(annotation.sequence)

                df = pd.DataFrame({
                    'Final_Annotation': final_annotation,
                    'Public_Notes': public_notes,
                    'Private_Notes': private_notes,
                    'Flag': flag,
                    'Protein_Sequence': protein_sequence
                })
                writer = ExcelWriter(file_path)
                df.to_excel(writer, 'Annotations_excel', index=False)
                writer.save()

                file = open(file_path, 'rb')
                response = HttpResponse(file, content_type='text/fasta')

                response['Content-Disposition'] = 'attachment; filename=Annotations_excel.xlsx'

            return response

        # return render(request, self.template_name, context)


# used to show user the correct template for uploading annotations
def download_excel_template(request):
    if request.user.is_authenticated:
        with TemporaryDirectory() as temp:
            file_path = os.path.join(temp, 'excel_template.xlsx')
            df = pd.DataFrame({
                'Final_Annotation': [''],
                'Public_Notes':  [''],
                'Private_Notes':  [''],
                'Flag': [''],
                'Protein_Sequence': ['']
            })

            writer = ExcelWriter(file_path)
            df.to_excel(writer, 'excel_template', index=False)
            writer.save()

            file = open(file_path, 'rb')
            response = HttpResponse(file, content_type='text/fasta')

            response['Content-Disposition'] = 'attachment; filename=excel_template.xlsx'
        return response


# Will download all of the users annotations into an excel file with the correct headers to re-upload them
def download_excel_annotations(request):
    user = request.user
    if user.is_authenticated():
        annotations = starship_models.Annotation.objects.filter(assigned_to=user)
        with TemporaryDirectory() as temp:
            file_path = os.path.join(temp, '%s_excel.xlsx' % user.first_name)

            final_annotation = []
            public_notes = []
            private_notes = []
            flag = []
            protein_sequence = []
            for annotation in annotations:
                final_annotation.append(annotation.annotation)
                public_notes.append(annotation.public_notes)
                private_notes.append(annotation.private_notes)
                flag.append(annotation.get_flag_display())
                protein_sequence.append(annotation.sequence)

            df = pd.DataFrame({
                'Final_Annotation': final_annotation,
                'Public_Notes': public_notes,
                'Private_Notes': private_notes,
                'Flag': flag,
                'Protein_Sequence': protein_sequence
            })
            writer = ExcelWriter(file_path)
            df.to_excel(writer, '%s_excel' % user.first_name, index=False)
            writer.save()

            file = open(file_path, 'rb')
            response = HttpResponse(file, content_type='text/fasta')

            response['Content-Disposition'] = 'attachment; filename=%s_excel.xlsx' % user.first_name
        return response


# Will download all the users unannotated annotations into an excel file with the correct headers to re-upload them
def download_unannotated_annotations(request):
    user = request.user
    if user.is_authenticated():
        annotations = starship_models.Annotation.objects.filter(assigned_to=user,
                                                              flag=7,)
        with TemporaryDirectory() as temp:
            file_path = os.path.join(temp, '%s_excel.xlsx' % user.first_name)

            final_annotation = []
            public_notes = []
            private_notes = []
            flag = []
            protein_sequence = []
            for annotation in annotations:
                final_annotation.append(annotation.annotation)
                public_notes.append(annotation.public_notes)
                private_notes.append(annotation.private_notes)
                flag.append(annotation.get_flag_display())
                protein_sequence.append(annotation.sequence)

            df = pd.DataFrame({
                'Final_Annotation': final_annotation,
                'Public_Notes': public_notes,
                'Private_Notes': private_notes,
                'Flag': flag,
                'Protein_Sequence': protein_sequence
            })
            writer = ExcelWriter(file_path)
            df.to_excel(writer, '%s_excel' % user.first_name, index=False)
            writer.save()

            file = open(file_path, 'rb')
            response = HttpResponse(file, content_type='text/fasta')

            response['Content-Disposition'] = 'attachment; filename=%s_excel.xlsx' % user.first_name
        return response

    else:
        raise


# used in starship list class
def get_starship_data_dicts(genomes):
    objects = []
    for genome in genomes:
        starship_dict = {}
        starship_dict['starship_name'] = starship.starship_name
        starship_dict['species'] = starship.species
        starship_dict['contigID'] = starship.contigID
        starship_dict['elementBegin'] = starship.elementBegin
        starship_dict['elementEnd'] = starship.elementEnd
        starship_dict["starship_family"] = starship.starship_family
        starship_dict["starship_navis"] = starship.starship_navis
        starship_dict["starship_haplotype"] = starship.starship_haplotype
        gene = starship.feature_set.filter(type='gene').count()
        gene_features = starship.feature_set.filter(type='gene')
        # Derived from stackoverflow.com/questions/4727327/
        flag_options_reverse = dict((v, k) for k, v in starship_models.Annotation.flag_options)
        annotations = starship_models.Annotation.objects.filter(feature__in=gene_features)
        starship_dict['unpolished_gene_count'] = annotations.filter(flag=flag_options_reverse['UNANNOTATED']).count()
        starship_dict['green_gene_count'] = annotations.filter(flag=flag_options_reverse['GREEN']).count()
        starship_dict['yellow_gene_count'] = annotations.filter(flag=flag_options_reverse['YELLOW']).count()
        starship_dict['red_gene_count'] = annotations.filter(flag=flag_options_reverse['RED']).count()
        starship_dict['review_name_gene_count'] = annotations.filter(flag=flag_options_reverse['REVIEW NAME']).count()
        starship_dict['gene_count'] = gene
        starship_dict['starship_length'] = len(starship.starship_sequence)
        objects.append(starship_dict)
    return objects

# Clear cache on save of any database
@receiver(post_save)
def post_save_delete(sender, **kwargs):
    # print('post save')
    cache.clear()

# Clear Cache on delete of any database
@receiver(post_delete)
def post_delete(sender, **kwargs):
    cache.clear()


def get_annotation_editors():
    return User.objects.filter(
        Q(is_superuser=True) |
        Q(user_permissions__codename='change_annotation') |
        Q(groups__permissions__codename='change_annotation')
    ).distinct()


# Will download the fasta file of users desired starship
def starship_download_fasta(request, starship_id):
    if request.user.is_authenticated:
        context = {}
        genome = starship_models.Starship.objects.get(pk=starship_id)
        starship_name = starship.starship_name
        nucleotide = starship.starship_sequence
        sequence = SeqRecord(Seq(nucleotide, IUPAC.ambiguous_dna), id=starship_name, description=starship_name)

        with TemporaryDirectory() as temp:
            file_path = os.path.join(temp, "%s.fsa" % starship_name)
            SeqIO.write(sequence, file_path, "fasta")
            file = open(file_path, 'rb')
            response = HttpResponse(file, content_type='text/fasta')

            response['Content-Disposition'] = 'attachment; filename=%s.fsa' % starship_name
        return response


# Will download the fasta file of users desired starship
def download_deliverables(request, starship_id):
    if request.user.is_authenticated:
        starship_name = starship_models.Starship.objects.get(pk=starship_id).starship_name

        with TemporaryDirectory() as tempdir:
            # Create deliverables in temp directory
            create_deliverables_args = argparse.Namespace(
                output_folder=tempdir,
                starship_name=starship_name,
                deliverables=[],
            )
            create_deliverables.main(create_deliverables_args)

            # Add each deliverable to a tarball
            def rename_extension(file):
                if file.endswith('.gbf'):
                    new_name = '{filename}.{extension}'.format(filename=file[:-4], extension='gbk')
                    os.rename(os.path.join(tempdir, file), os.path.join(tempdir, new_name))
                    return new_name
                return file

            file_list = [rename_extension(f) for f in os.listdir(tempdir) if os.path.isfile(os.path.join(tempdir, f)) and f[-4:] not in ['.txt', '.sqn', '.tbl']]
            tar_path = os.path.join(tempdir, '%s_deliverables.tar.gz' % starship_name)

            with tarfile.open(tar_path, mode='w:gz') as tar:
                for f in file_list:
                    tar.add(os.path.join(tempdir, f), arcname=f)

            # Send to user
            file = open(tar_path, 'rb')
            response = HttpResponse(file, content_type='text/tar')
            response['Content-Disposition'] = 'attachment; filename=%s_deliverables.tar.gz' % starship_name

        return response


# Will download a single annotation for the user
def annotation_download(request, annotation_id):
    if request.user.is_authenticated:
        annotation_obj = starship_models.Annotation.objects.get(pk=annotation_id)
        annotation = annotation_obj.annotation
        aa_sequence = annotation_obj.sequence
        public_note = annotation_obj.public_notes
        private_note = annotation_obj.private_notes
        flag = annotation_obj.get_flag_display()

        sequence = SeqRecord(Seq(aa_sequence, generic_protein), id=annotation_obj.accession + " |",
                             description="%s | %s | %s | %s" % (annotation, public_note, private_note, flag))

        with TemporaryDirectory() as temp:
            file_path = os.path.join(temp, "%s.faa" % annotation_obj.accession)
            SeqIO.write(sequence, file_path, "fasta")

            file = open(file_path, 'rb')
            response = HttpResponse(file, content_type='text/fasta')

            response['Content-Disposition'] = 'attachment; filename=%s.faa' % annotation_obj.accession

        return response
