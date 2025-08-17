from django.shortcuts import render
from django.views import generic
from django.urls import reverse, reverse_lazy
from django.contrib.auth.mixins import LoginRequiredMixin, PermissionRequiredMixin
from django.core.exceptions import ObjectDoesNotExist
from django.http import Http404, HttpResponse
from kombu.exceptions import OperationalError

import math
from copy import copy

from Bio.Blast import NCBIXML

from MAS.celery import app
from result_viewer.forms import AnnotationForm, GenomeSearchForm
from result_viewer.hhsuite2_text import Hhsuite2TextParser
from result_viewer.models import *
from result_viewer.navigator import FlagNavigator, GenomeNavigator, AssignmentNavigator
from starship.models import *
from result_viewer.utils import interproscan_xml_to_dict


def blastp_alignment_to_str(alignment):
    '''
    The __str__ method for HSP shows only part of the alignment
    :param hsp: biopython hsp object with a single HSPFragment
    :return: String representation of it
    '''
    text = ''
    line_length = 50

    for i, hsp in enumerate(alignment.hsps):
        spacer = max([len(str(hsp.query_end)), len(str(hsp.sbjct_end)),
                      len(str(hsp.query_start)), len(str(hsp.sbjct_start))]) + 2
        sstrand = hsp.sbjct_start < hsp.sbjct_end
        qstrand = hsp.query_start < hsp.query_end
        text += '{:*^{width}}\n'.format(
            ' HSP {}, e-value={:.2e}, bit-score={} '.format(i, hsp.expect, hsp.bits),
            width=line_length + spacer * 2 + 8
        )

        for row_i in range(math.ceil(len(hsp.query) / line_length)):
            start = row_i * line_length
            stop = (row_i+1) * line_length

            query_seq = hsp.query[start:stop]
            query_align_length = len(query_seq) - query_seq.count('-')
            qline_start = hsp.query_start + start if qstrand else hsp.query_start - start
            text += 'query {:<{width}} {} ({})\n'.format(
                '({})'.format(qline_start),
                query_seq,
                qline_start + query_align_length if qstrand else qline_start - query_align_length,
                width=spacer
            )

            text += '      {:<{width}} {}\n'.format('', hsp.match[start:stop], width=spacer)

            sbjct_seq = hsp.sbjct[start:stop]
            sbjct_align_length = len(sbjct_seq) - sbjct_seq.count('-')
            sline_start = hsp.sbjct_start + start if sstrand else hsp.sbjct_start - start
            text += 'sbjct {:<{width}} {} ({})\n\n'.format(
                '({})'.format(sline_start),
                sbjct_seq,
                sline_start + sbjct_align_length if sstrand else sline_start - sbjct_align_length,
                width=spacer
            )

    return text


def hhsearch_hsp_to_dict(hsp):
    '''
    :param hsp: biopython hsp object from hhsearch parsing
    :return: A dict wich can be used in jinja templating and serialized to json
    '''
    d = {
        'hit_id': hsp.hit_id,
        'pdb_chain_name': None,
        'hit_seq_len': hsp.hit_seq_len,
        'hit_description': hsp.hit_description,
        'top_score': hsp.score,
        'sum_score': hsp.score,
        'lowest_evalue': hsp.evalue,
        'query_start': hsp.query_start,
        'query_end': hsp.query_end,
        'prob': hsp.prob,
        'text': hsp.text,
        'hits': [{
            'hit_start': hsp.hit_start,
            'hit_end': hsp.hit_end,
            'query_start': hsp.query_start,
            'query_end': hsp.query_end,
            'evalue': hsp.evalue,
            'score': hsp.score,
            'prob': hsp.prob
        }]
    }

    try:
        mapping = PDB_Accession_Mapping.objects.get(pdb_accession=hsp.hit_id)
        if mapping.pdb_chain_name:
            d['pdb_chain_name'] = mapping.pdb_chain_name
    except ObjectDoesNotExist:
        pass

    return d


def blastp_alignment_to_dict(alignment):
    '''
    :param hsp: biopython Alignment object from NCBIXML parsing (Don't use SearchIO)
    :return: A dict which can be used in jinja templating and serialized to json
    '''
    hits = []
    top_score = 0
    sum_score = 0
    lowest_evalue = 100
    query_start = None
    query_end = None

    for hsp in alignment.hsps:
        hits.append({
            'hit_start': hsp.sbjct_start,
            'hit_end': hsp.sbjct_end,
            'query_start': hsp.query_start,
            'query_end': hsp.query_end,
            'evalue': hsp.expect,
            'score': hsp.bits
        })
        lowest_evalue = hsp.expect if hsp.expect < lowest_evalue else lowest_evalue
        top_score = hsp.bits if hsp.bits > top_score else top_score
        query_start = hsp.query_start if query_start is None or hsp.query_start < query_start else query_start
        query_end = hsp.query_end if query_end is None or hsp.query_end > query_end else query_end
        sum_score += hsp.bits

    return {
        'hit_id': alignment.accession,
        'hit_description': alignment.title,
        'hit_seq_len': alignment.length,
        'top_score': top_score,
        'sum_score': sum_score,
        'lowest_evalue': lowest_evalue,
        'query_start': query_start,
        'query_end': query_end,
        'hits': hits,
        'text': blastp_alignment_to_str(alignment)
    }


def add_context_for_starship_viz(context, starship, current_annotation_id=None):
    """
    Add visualization context for a starship
    :param context: The context dictionary to update
    :param starship: The starship object to visualize
    :param current_annotation_id: Optional ID of the current annotation
    :return: Updated context dictionary
    """
    context['starship_id'] = starship.id
    context['starship_length'] = len(starship.starship_sequence)
    feature_objects = Feature.objects.filter(starship__id=starship.id)
    # current_annotation_id = int(self.kwargs['accession'], 36)
    features = []
    features_dict = {}
    for feature in feature_objects:
        f = {}
        f['id'] = feature.id
        f['strand'] = feature.strand
        try:
            f['flag'] = feature.annotation.get_flag_display()
        except AttributeError:
            f['flag'] = 7
        f['type'] = feature.type
        try:
            f['public_note'] = feature.annotation.public_notes
        except AttributeError:
            f['public_note'] = "nan"

        try:
            f['private_note'] = feature.annotation.private_notes
        except AttributeError:
            f['private_note'] = "nan"

        try:
            f['private_note'] = feature.annotation.private_notes
        except AttributeError:
            f['private_note'] = "nan"

        try:
            f['accession'] = feature.annotation.accession
        except AttributeError:
            f['accession'] = "nan"
        try:
            f['annotation_id'] = feature.annotation.id
        except AttributeError:
            f['annotation_id'] = "nan"
        try:
            f['annotation'] = feature.annotation.annotation
        except AttributeError:
            f['annotation'] = "nan"
        f['href'] = reverse('view-results', args=(f['accession'], 'GenomeNavigator', starship.starship_name))
        if f['annotation_id'] == current_annotation_id:
            features_dict['feature_id'] = f['id']

        # Deal with genes split over genome's ends
        if feature.start > feature.stop and feature.strand == '+':
            f2 = copy(f)
            f['start'] = feature.start
            f['stop'] = context['starship_length']
            f2['start'] = 0
            f2['stop'] = feature.stop
            features.append(f2)
        else:
            f['start'] = feature.start
            f['stop'] = feature.stop

        features.append(f)
    features_dict['features'] = features
    context['feature_data'] = features_dict
    return context


class MixinForBaseTemplate(generic.base.ContextMixin):

    def get_context_data(self, **kwargs):
        context = super(MixinForBaseTemplate, self).get_context_data(**kwargs)

        # Determine if Dark mode activated
        if self.request.user.is_authenticated and hasattr(self.request.user, 'userpreferences'):
                context['dark_mode'] = self.request.user.userpreferences.dark_mode_activated
        else:
            context['dark_mode'] = self.request.session.get('dark_mode', False)

        context['starship_search_form'] = GenomeSearchForm()
        # Safely check worker status
        try:
            context['worker_active'] = bool(app.control.inspect(settings.CELERY_WORKERS).ping())
        except (OperationalError, AttributeError, ConnectionError):
            # If we can't connect to the message broker or CELERY_WORKERS isn't defined
            context['worker_active'] = False
        return context


class ChangeTheme(MixinForBaseTemplate, generic.View):
    def post(self, request):
        if self.request.user.is_authenticated:
            if hasattr(self.request.user, 'userpreferences'):
                prefs = self.request.user.userpreferences
                prefs.dark_mode_activated = not prefs.dark_mode_activated
                prefs.save()
            else:
                prefs = UserPreferences(
                    user=self.request.user,
                    dark_mode_activated=not self.request.session.get('dark_mode', False)
                )
                prefs.save()
        else:
            self.request.session['dark_mode'] = not self.request.session.get('dark_mode', False)
        return HttpResponse('Theme changed')


class ViewNoResults(LoginRequiredMixin, MixinForBaseTemplate, generic.TemplateView):
    login_url = reverse_lazy('login')
    template_name = 'result_viewer/no_results.html'

    def get_context_data(self, **kwargs):
        context = super(ViewNoResults, self).get_context_data(**kwargs)
        context['nav_arg'] = self.kwargs['nav_arg']
        context['navigator'] = self.kwargs['navigator']

        return context


# Create your views here.
class ViewResults(LoginRequiredMixin, MixinForBaseTemplate, generic.UpdateView):
    login_url = reverse_lazy('login')
    template_name = 'result_viewer/view_results.html'
    model = Annotation
    # pk_url_kwarg = 'aa'
    form_class = AnnotationForm

    def get_object(self, queryset=None):
        return self.model.objects.get(pk=int(self.kwargs['accession'], 36))

    def get_context_data(self, **kwargs):
        context = super(ViewResults, self).get_context_data(**kwargs)

        # URLs
        context['run_search_url'] = reverse('run_search')
        context['run_search_for_starship_url'] = reverse('run_search_for_starship')

        # Default tool and database
        context['tool'] = 'blastp'
        context['database'] = 'swissprot'

        # Query info
        context['accession'] = self.kwargs['accession']
        context['query_len'] = len(context['object'].sequence)

        # Navigation context
        context['navigator'] = self.nav.as_context()
        context['navigator_index'] = self.nav.idx
        context['navigator_index'] = 3

        # History info
        context['history'] = context['annotation'].history.all().order_by('history_date')

        # Starship visualization
        if context['navigator']['type'] == 'GenomeNavigator':
            starship = JoinedShips.objects.get(starship_name=self.kwargs['nav_arg'])
            context['starship_id'] = starship.id
            if starship.feature_set.count() < 1000:
                current_annotation_id = int(self.kwargs['accession'], 36)
                context = add_context_for_starship_viz(context, starship, current_annotation_id)

        # Disable form fields if user has no permission to edit
        if not self.request.user.has_perm('starship.change_annotation'):
            for field_name in context['form'].fields:
                context['form'].fields[field_name].disabled = True

        # ######  Try to load hhsearch results  #######
        context['hhsearch_alignments'] = {}

        for db in HHSearch_Result.database_options:
            db = db[0]
            context['hhsearch_alignments'][db] = {'alignments': None, 'date_ran': None, 'status': None}

            try:
                hhsearch_result = HHSearch_Result.objects.get(annotation=context['annotation'], database=db)

                if hhsearch_result.result:
                    for hhr in Hhsuite2TextParser(hhsearch_result.result.open(mode='r')):
                        break
                    context['hhsearch_alignments'][db]['alignments'] = {x.hit_id: hhsearch_hsp_to_dict(x) for x in hhr.hsps}

                else:
                    context['hhsearch_alignments'][db]['alignments'] = []

                context['hhsearch_alignments'][db]['date_ran'] = hhsearch_result.run_date
                context['hhsearch_alignments'][db]['status'] = hhsearch_result.status

            except ObjectDoesNotExist:
                pass

        # ######  Try to load blastp results  #######
        context['blastp_alignments'] = {}

        for db in Blastp_Result.database_options:
            db = db[0]
            context['blastp_alignments'][db] = {'alignments': None, 'date_ran': None, 'status': None}

            try:
                blastp_result = Blastp_Result.objects.get(annotation=context['annotation'], database=db)

                if blastp_result.result:
                    record = NCBIXML.read(blastp_result.result.open(mode='r'))
                    context['blastp_alignments'][db]['alignments'] = {x.accession: blastp_alignment_to_dict(x) for x in
                                                                      record.alignments}

                else:
                    context['blastp_alignments'][db]['alignments'] = []

                context['blastp_alignments'][db]['date_ran'] = blastp_result.run_date
                context['blastp_alignments'][db]['status'] = blastp_result.status

            except ObjectDoesNotExist:
                pass

        # ######  Try to load rpsblast results  #######
        context['rpsblast_alignments'] = {}

        for db in RPSBlast_Result.database_options:
            db = db[0]
            context['rpsblast_alignments'][db] = {'alignments': None, 'date_ran': None, 'status': None}

            try:
                rpsblast_result = RPSBlast_Result.objects.get(annotation=context['annotation'], database=db)

                if rpsblast_result.result:
                    record = NCBIXML.read(rpsblast_result.result.open(mode='r'))
                    context['rpsblast_alignments'][db]['alignments'] = {x.accession: blastp_alignment_to_dict(x) for x in record.alignments}

                else:
                    context['rpsblast_alignments'][db]['alignments'] = []

                context['rpsblast_alignments'][db]['date_ran'] = rpsblast_result.run_date
                context['rpsblast_alignments'][db]['status'] = rpsblast_result.status

            except ObjectDoesNotExist:
                pass

        # Add InterProScan results
        context['interpro_alignments'] = {}
        
        for db in Interpro_Result.database_options:
            db = db[0]
            context['interpro_alignments'][db] = {'matches': None, 'date_ran': None, 'status': None}
            
            try:
                interpro_result = Interpro_Result.objects.get(annotation=context['annotation'], database=db)
                
                if interpro_result.result:
                    context['interpro_alignments'][db]['matches'] = interproscan_xml_to_dict(
                        interpro_result.result.open(mode='r')
                    )
                else:
                    context['interpro_alignments'][db]['matches'] = []
                
                context['interpro_alignments'][db]['date_ran'] = interpro_result.run_date
                context['interpro_alignments'][db]['status'] = interpro_result.status
                
            except ObjectDoesNotExist:
                pass

        return context

    def dispatch(self, request, *args, **kwargs):
        # Set up Navigator
        if self.kwargs['navigator'] == 'FlagNavigator':
            self.nav = FlagNavigator(int(self.kwargs['nav_arg']), self.kwargs['accession'])

        elif self.kwargs['navigator'] == 'GenomeNavigator':
            self.nav = GenomeNavigator(self.kwargs['nav_arg'], self.kwargs['accession'])

        elif self.kwargs['navigator'] == 'AssignmentNavigator':
            self.nav = AssignmentNavigator(self.kwargs['nav_arg'], self.kwargs['accession'])

        else:
            raise Http404('Invalid Navigator')

        return super().dispatch(request, *args, **kwargs)

    def get_success_url(self):
        if 'go_to_next' not in self.request.POST or self.nav.next() is None:
            if self.kwargs['navigator'] == 'FlagNavigator':
                navarg = self.object.flag
            elif self.kwargs['navigator'] == 'GenomeNavigator':
                navarg = self.kwargs['nav_arg']
            elif self.kwargs['navigator'] == 'AssignmentNavigator':
                if self.request.POST['assigned_to'] == self.kwargs['nav_arg']:
                    navarg = self.kwargs['nav_arg']
                else:
                    return reverse('accession-redirect', args=(self.kwargs['accession'],))

            return reverse('view-results', args=(self.kwargs['accession'], self.kwargs['navigator'], navarg))
        else:
            return self.nav.next()

    def post(self, request, *args, **kwargs):
        if request.user.has_perm('starship.change_annotation'):
            return super().post(request, *args, **kwargs)
        else:
            raise PermissionError('Current user does not have permissions to change annotations')


# class ProteinDoesNotExistView(MixinForBaseTemplate, generic.TemplateView):
#     template_name = 'home/index.html'
#
#     def get_context_data(self, **kwargs):
#         context = super().get_context_data(**kwargs)
#
#         return context


class FlagNavRedirect(generic.RedirectView):
    def get_redirect_url(self, *args, **kwargs):
        qs = Annotation.objects.filter(flag=kwargs['flag']).order_by('pk')

        if qs.count() > 0:
            first_obj = qs[0]
            return reverse('view-results', args=(first_obj.accession, 'FlagNavigator', kwargs['flag']))

        else:
            return reverse('no-results', args=('FlagNavigator', kwargs['flag']))
            # raise ObjectDoesNotExist('Starship has no annotations.')


class GenomeNavRedirect(generic.RedirectView):
    def get_redirect_url(self, *args, **kwargs):
        try:
            first_obj = GenomeNavigator(kwargs['starship_name']).get_current_object()
            return reverse('view-results', args=(first_obj.accession, 'GenomeNavigator', kwargs['starship_name']))
        except:
            return reverse('no-results', args=('GenomeNavigator', kwargs['starship_name']))


class AssignmentNavRedirect(generic.RedirectView):
    def get_redirect_url(self, *args, **kwargs):
        nav = AssignmentNavigator(kwargs['user'])

        if nav.size > 0:
            first_obj = nav.get_current_object()
            return reverse('view-results', args=(first_obj.accession, 'AssignmentNavigator', kwargs['user']))
        else:
            return reverse('no-results', args=('AssignmentNavigator', kwargs['user']))

            # raise ObjectDoesNotExist('Error: User %s does not have any annotations assigned' % kwargs['user'])


class AccessionRedirect(generic.RedirectView):
    def get_redirect_url(self, *args, **kwargs):
        '''
        If only the accession is provided we will use the flag navigator becuase one annotation
        could belong to multiple starships and an annotation may not be assigned to anyone.
        '''
        obj = Annotation.objects.get(pk=int(kwargs['accession'], 36))
        return reverse('view-results', args=(kwargs['accession'], 'FlagNavigator', obj.flag))