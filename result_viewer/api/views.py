from datetime import datetime
import subprocess
import json
from tempfile import NamedTemporaryFile

from django.http import Http404, HttpResponseBadRequest
from django.conf import settings
from django.db.models import Q, F, Func, Count, Value, When, Case, CharField
from django.db.models.functions import Length
from django.urls import reverse

from rest_framework import viewsets, permissions, decorators, parsers, response, status, exceptions
from rest_framework.views import APIView
from rest_framework.parsers import FormParser

from result_viewer.api.serializers import *
from result_viewer.api.tasks import run_single_search, run_multiple_search

from starship.models import Annotation, Feature
from starship.starbase_models import JoinedShips

class RunSearchAjaxView(APIView):
    '''
    This view is called when a user selects the 'Run' button
    '''

    def post(self, request, format=None):
        s = RunSearchAjaxSerializer(data=request.data)

        if s.is_valid():
            accession = s.data['accession']
            tool = s.data['tool']
            database = s.data['database']

            uri_scheme = 'https://' if request.is_secure() else 'http://'
            mas_server_uri = uri_scheme + request.site.domain
            
            # Start the search task
            run_single_search.delay(accession, tool, database, mas_server_uri)

            # Return success response with status 1 (running)
            return response.Response({'status': 1, 'message': 'Search started successfully'})

        return response.Response(s.errors, status=status.HTTP_400_BAD_REQUEST)

    def get(self, request, format=None):
        '''
        Check the status of a search task
        '''
        accession = request.GET.get('accession')
        tool = request.GET.get('tool')
        database = request.GET.get('database')

        if not all([accession, tool, database]):
            return response.Response({'error': 'Missing required parameters'}, status=status.HTTP_400_BAD_REQUEST)

        try:
            annotation = Annotation.objects.get(pk=int(accession, 36))
            
            # Get the appropriate result model
            result_model = None
            if tool == 'blastp':
                result_model = Blastp_Result
            elif tool == 'hhsearch':
                result_model = HHSearch_Result
            elif tool == 'rpsblast':
                result_model = RPSBlast_Result
            else:
                return response.Response({'error': 'Invalid tool'}, status=status.HTTP_400_BAD_REQUEST)

            # Check if result exists
            try:
                result_obj = result_model.objects.get(annotation=annotation, database=database)
                return response.Response({
                    'status': result_obj.status,
                    'run_date': result_obj.run_date,
                    'message': 'Status retrieved successfully'
                })
            except result_model.DoesNotExist:
                # No result entry exists, so search was never run
                return response.Response({
                    'status': None,
                    'run_date': None,
                    'message': 'Search never run'
                })

        except Annotation.DoesNotExist:
            return response.Response({'error': 'Invalid accession'}, status=status.HTTP_400_BAD_REQUEST)
        except ValueError:
            return response.Response({'error': 'Invalid accession format'}, status=status.HTTP_400_BAD_REQUEST)

    def check_permissions(self, request):
        super().check_permissions(request)
        if not (request.user.is_superuser or request.user.groups.filter(name='Data Editors').exists()):
            self.permission_denied(
                request, message='Permission Denied: Only data editors/superusers can use this API view'
            )


class RunAllStarshipProteinsAjaxView(APIView):

    def post(self, request, format=None):
        s = RunAllStarshipProteinsAjaxSerializer(data=json.loads(request.data['data']))

        if s.is_valid():
            # Get data from serializer
            starship_name = s.data['starship']
            rerun = s.data['rerun']
            tools_and_databases = s.data['tools_and_databases']

            uri_scheme = 'https://' if request.is_secure() else 'http://'
            mas_server_uri = uri_scheme + request.site.domain

            run_multiple_search.delay(starship_name, rerun, tools_and_databases, mas_server_uri)

            return response.Response(s.data)

        return response.Response(s.errors, status=status.HTTP_400_BAD_REQUEST)

    def check_permissions(self, request):
        super().check_permissions(request)
        if not (request.user.is_superuser or request.user.groups.filter(name='Data Editors').exists()):
            self.permission_denied(
                request, message='Permission Denied: Only data editors/superusers can use this API view'
            )


class PipelineAPIMixin(object):
    '''
    Mixin for adding correct permissions to views for use by luigi pipeline
    '''
    def check_permissions(self, request):
        super().check_permissions(request)
        if not (request.user.is_superuser or request.user.username == 'luigi'):
            self.permission_denied(
                request, message='Permission Denied: Only Luigi/superusers can use this API view'
            )


class TestConnectionView(PipelineAPIMixin, APIView):
    def get(self, request):
        return response.Response()


class GetProtSeqView(PipelineAPIMixin, APIView):
    permission_classes = [permissions.DjangoModelPermissions]
    queryset = Annotation.objects.all()

    def get(self, request, accession, format=None):
        try:
            annotation = self.queryset.get(pk=int(accession, 36))

        except Annotation.DoesNotExist:
            raise Http404

        s = ProteinSeqSerializer(annotation)
        return response.Response(s.data)


class GetStarshipView(PipelineAPIMixin, APIView):
    permission_classes = [permissions.DjangoModelPermissions]
    queryset = JoinedShips.objects.all()

    def get(self, request, starship_name, format=None):
        try:
            starship = self.queryset.get(starshipID=starship_name)
            starship_features = Feature.objects.filter(starship=starship)
            counts = starship_features.aggregate(
                tRNAs=Count('id', filter=Q(type='tRNA')),
                repeats=Count('id', filter=Q(type='Repeat Region')),
                gene=Count('id', filter=Q(type='gene')),
                cds=Count('id', filter=Q(type='CDS'))
            )
            if counts['repeats'] > 0:
                DTR_length = len(starship_features.get(type='Repeat Region').annotation.sequence)
            else:
                DTR_length = 0

        except JoinedShips.DoesNotExist:
            raise Http404

        s = StarshipSeqSerializer({
            'starship_name': starship.starshipID,
            'starship_sequence': starship.starship_sequence,
            'num_cds': counts['cds'],
            'num_gene': counts['gene'],
            'num_trna': counts['tRNAs'],
            'len_dtr': DTR_length
        })
        return response.Response(s.data)


class UploadResultsView(PipelineAPIMixin, APIView):
    parser_classes = (parsers.MultiPartParser, parsers.FormParser,)

    def post(self, request):
        s = UploadResultsSerializer(data=request.data)

        if s.is_valid():
            tool = s.data['tool']
            annotation = Annotation.objects.get(pk=int(s.data['accession'], 36))

            def save_result(result_class):
                obj = result_class.objects.get(annotation=annotation, database=s.data['database'])

                if 'result' in request.data:
                    obj.result.save(name=s.data['accession'], content=request.data['result'])

                obj.status = s.data['status']
                obj.save()

            if tool == 'blastp':
                save_result(Blastp_Result)

            elif tool == 'hhsearch':
                save_result(HHSearch_Result)

            elif tool == 'rpsblast':
                save_result(RPSBlast_Result)

            return response.Response(s.data)

        return response.Response(s.errors, status=status.HTTP_400_BAD_REQUEST)


class CombinedStarshipData:
    """A simple object that provides attribute access to combined starship data"""
    def __init__(self, starship_dict, feature_counts):
        # Store the starship data and feature counts
        self._data = starship_dict
        self._counts = feature_counts
    
    def __getattr__(self, name):
        # First check if it's in the feature counts
        if name in self._counts:
            return self._counts[name]
        # Then check if it's in the starship data
        elif name in self._data:
            return self._data[name]
        # Return None for missing attributes
        else:
            return None


class StarshipData:
    def __init__(self, starship):
        self.starship_name = '<a href="{url}">{name}</a>'.format(
            url=reverse('starship:starship_detail', kwargs={'pk': starship.id}),
            name=starship.starshipID
        )
        
        # Now we can use the aggregated counts from the query
        self.species = getattr(starship, 'taxonomy_species', 'N/A')
        self.starship_length = getattr(starship, 'starship_length', 0)
        self.num_cds = getattr(starship, 'num_cds', 0)
        self.num_gene = getattr(starship, 'num_gene', 0)
        self.num_unannotated = getattr(starship, 'num_unannotated', 0)
        self.num_review = getattr(starship, 'num_review', 0)
        self.num_green = getattr(starship, 'num_green', 0)
        self.num_yellow = getattr(starship, 'num_yellow', 0)
        self.num_red = getattr(starship, 'num_red', 0)
        self.contigID = getattr(starship, 'contigID', 'N/A')
        self.elementBegin = getattr(starship, 'elementBegin', 0)
        self.elementEnd = getattr(starship, 'elementEnd', 0)
        self.starship_family = getattr(starship, 'starship_family', 'N/A')
        self.starship_navis = getattr(starship, 'starship_navis', 'N/A')
        self.starship_haplotype = getattr(starship, 'starship_haplotype', 'N/A')
        
        self.download = '<a href="{}">download fasta</a>'.format(
            reverse('starship:starship_download_fasta', kwargs={'starship_id': starship.id})
        )
        self.navigator = '<a href="{}"><div class="glyphicon glyphicon-hand-right"></div></a>'.format(
            reverse('starship-nav-redirect', kwargs={'starship_name': starship.starshipID})
        )


class GetStarshipDataView(APIView):
    def get(self, request):
        params = self.read_parameters(request.GET)

        # Query the starbase database for JoinedShips
        starships = JoinedShips.objects.select_related(
            'ship_family_id', 'tax_id', 'ship_navis_id', 'ship_haplotype_id', 'ship_id'
        ).annotate(
            # Get basic starship data
            starship_id=F('starshipID'),
            taxonomy_species=F('tax_id__name'),
            starship_family=F('ship_family_id__familyName'),
            starship_navis=F('ship_navis_id__navis_name'),
            starship_haplotype=F('ship_haplotype_id__haplotype_name'),
            starship_length=Length('ship_id__sequence'),
            contig_id=F('ship_id__starshipfeatures__contigID'),
            element_begin=F('ship_id__starshipfeatures__elementBegin'),
            element_end=F('ship_id__starshipfeatures__elementEnd'),
        )

        total_num_starships = starships.count()

        # Filter starships
        starships = starships.filter(
            Q(starshipID__icontains=params['search_val']) |
            Q(tax_id__species__icontains=params['search_val'])
        ).order_by(self.get_order_by_arg(params['order_col'], params['order_dir']))

        filtered_num_starships = starships.count()
        starships = starships[params['start'] : params['start']+params['length']]

        # Now get the feature counts from the mas database
        # Get the actual JoinedShips IDs from the starbase database
        starship_joined_ids = [s.id for s in starships]  # These are the JoinedShips IDs
        
        # Query the mas database for features using the JoinedShips IDs directly
        from starship.models import Feature, Annotation
        
        # Get feature counts for each starship
        feature_counts = {}
        for joined_id in starship_joined_ids:
            # Count features by type - use the JoinedShips ID directly
            cds_count = Feature.objects.filter(
                starship_id=joined_id,
                type='CDS'
            ).count()
            
            gene_count = Feature.objects.filter(
                starship_id=joined_id,
                type='gene'
            ).count()
            
            # Count annotations by flag
            unannotated_count = Annotation.objects.filter(
                feature__starship_id=joined_id,
                flag=7
            ).count()
            
            review_count = Annotation.objects.filter(
                feature__starship_id=joined_id,
                flag=3
            ).count()
            
            green_count = Annotation.objects.filter(
                feature__starship_id=joined_id,
                flag=0
            ).count()
            
            yellow_count = Annotation.objects.filter(
                feature__starship_id=joined_id,
                flag=1
            ).count()
            
            red_count = Annotation.objects.filter(
                feature__starship_id=joined_id,
                flag=2
            ).count()
            
            feature_counts[joined_id] = {
                'num_cds': cds_count,
                'num_gene': gene_count,
                'num_unannotated': unannotated_count,
                'num_review': review_count,
                'num_green': green_count,
                'num_yellow': yellow_count,
                'num_red': red_count,
            }

        # Create starship data with feature counts
        starship_data = []
        for starship in starships:
            counts = feature_counts.get(starship.id, {  # Use starship.id instead of starship.starshipID
                'num_cds': 0, 'num_gene': 0, 'num_unannotated': 0,
                'num_review': 0, 'num_green': 0, 'num_yellow': 0, 'num_red': 0
            })
            
            # Create a combined object with both starbase and mas data
            combined_starship = CombinedStarshipData(starship.__dict__, counts)
            starship_data.append(StarshipData(combined_starship))

        return response.Response(StarshipDataListSerializer({
            'data': starship_data,
            'draw': params['draw'],
            'recordsTotal': total_num_starships,
            'recordsFiltered': filtered_num_starships
        }).data)

    def read_parameters(self, query_dict):
        """ Converts and cleans up the GET parameters. """
        return {
            'start': int(query_dict.get('start')),
            'length': int(query_dict.get('length')),
            'draw': int(query_dict.get('draw')),
            'order_col': int(query_dict.get('order[0][column]')),
            'order_dir': query_dict.get('order[0][dir]'),
            'search_val': query_dict.get('search[value]')
        }

    def get_order_by_arg(self, order_col, order_dir):
        if order_col == 0:
            return 'starshipID' if order_dir == 'asc' else '-starshipID'
        elif order_col == 1:
            return 'taxonomy_species' if order_dir == 'asc' else '-taxonomy_species'
        elif order_col == 2:
            return Length('starship_sequence').asc() if order_dir == 'asc' else Length('starship_sequence').desc()
        elif order_col == 3:
            return 'num_gene' if order_dir == 'asc' else '-num_gene'
        elif order_col == 4:
            return 'num_unannotated' if order_dir == 'asc' else '-num_unannotated'
        elif order_col == 5:
            return 'num_review' if order_dir == 'asc' else '-num_review'
        elif order_col == 6:
            return 'num_green' if order_dir == 'asc' else '-num_green'
        elif order_col == 7:
            return 'num_yellow' if order_dir == 'asc' else '-num_yellow'
        elif order_col == 8:
            return 'num_red' if order_dir == 'asc' else '-num_red'
        elif order_col == 9:
            return 'contigID' if order_dir == 'asc' else '-contigID'
        elif order_col == 10:
            return 'elementBegin' if order_dir == 'asc' else '-elementBegin'
        elif order_col == 11:
            return 'elementEnd' if order_dir == 'asc' else '-elementEnd'
        elif order_col == 12:
            return 'starship_family' if order_dir == 'asc' else '-starship_family'
        elif order_col == 13:
            return 'starship_navis' if order_dir == 'asc' else '-starship_navis'
        elif order_col == 14:
            return 'starship_haplotype' if order_dir == 'asc' else '-starship_haplotype'
        

class AnnotationData:
    def __init__(self, annotation, starship_name=None):
        self.sequence = '<div class="glyphicon glyphicon-menu-right aa-control details-control"></div>' \
                        '<input type="hidden" class="annotation" value="{pk}">'.format(pk=annotation.pk)
        self.feature = '<ul class="list-group">'
        for feature in annotation.feature_set.all():
            self.feature += '<li class="list-group-item list-group-item-text">'
            self.feature += '<a href="{url}">{description}</a>'.format(
                url=reverse('starship:feature_detail', kwargs={'pk': feature.pk}),
                description=feature
            )
        self.feature += '</ul>'
        self.accession = '<a href="{url}">{accession}</a>'.format(
            url=reverse('starship:annotation_detail', kwargs={'pk': annotation.pk}),
            accession=annotation.accession_sql
        )
        self.annotation = annotation.annotation
        self.length = len(annotation.sequence)
        self.public_notes = annotation.public_notes
        self.private_notes = annotation.private_notes
        self.flag = annotation.get_flag_display()
        self.assigned_to = annotation.assigned_to
        self.download = '<a href="{url}"><div class="glyphicon glyphicon-download"></div></a>'.format(
            url=reverse('starship:annotation_download', kwargs={'annotation_id': annotation.pk})
        )
        self.history = '<a href="{url}"><i class="fa fa-history"></i></a>'.format(
            url=reverse('starship:annotation_history', kwargs={'annotation_pk': annotation.pk})
        )

        if starship_name:
            view_results_url = reverse('view-results', kwargs={
                'navigator': 'StarshipNavigator',
                'accession': annotation.accession,
                'nav_arg': starship_name
            })
        else:
            view_results_url = reverse('accession-redirect', kwargs={'accession': annotation.accession})

        self.view_results = '<a href="{url}"><div class="glyphicon glyphicon-hand-right"></div></a>'.format(
            url=view_results_url
        )


class GetAnnotationListView(APIView):
    def get(self, request):
        params = self.read_parameters(request.GET)

        if params['starship_id']:
            annotations = Annotation.objects.filter(feature__starship_id=params['starship_id']).distinct()
            starship_name = JoinedShips.objects.get(id=params['starship_id']).starshipID
        else:
            annotations = Annotation.objects
            starship_name = None

        # Annotate each annotation with accession (Only works in mysql/mariadb)
        annotations = annotations.annotate(
            accession_sql=Func(Func(F('id'), 10, 36, function='CONV'), 5, 0, function='LPAD')
        )
        # Annotate with flag strings
        when_clauses = [When(flag=i, then=Value(flag)) for i, flag in Annotation.flag_options]
        annotations = annotations.annotate(
            flag_str=Case(*when_clauses, output_field=CharField())
        )

        total_num_annotations = annotations.count()

        annotations = annotations.filter(
            Q(flag_str__icontains=params['search_val']) |
            Q(sequence__icontains=params['search_val']) |
            Q(accession_sql__icontains=params['search_val']) |
            Q(annotation__icontains=params['search_val']) |
            Q(public_notes__icontains=params['search_val']) |
            Q(private_notes__icontains=params['search_val']) |
            Q(assigned_to__username__icontains=params['search_val'])
        ).order_by(self.get_order_by_arg(params['order_col'], params['order_dir']))

        filtered_num_annotations = annotations.count()

        annotations = annotations[params['start'] : params['start'] + params['length']]
        annotation_data = [AnnotationData(a, starship_name) for a in annotations]

        return response.Response(AnnotationDataListSerializer({
            'data': annotation_data,
            'draw': params['draw'],
            'recordsTotal': total_num_annotations,
            'recordsFiltered': filtered_num_annotations
        }).data)

    def read_parameters(self, query_dict):
        return {
            'start': int(query_dict.get('start')),
            'length': int(query_dict.get('length')),
            'draw': int(query_dict.get('draw')),
            'order_col': int(query_dict.get('order[0][column]')),
            'order_dir': query_dict.get('order[0][dir]'),
            'search_val': query_dict.get('search[value]'),
            'starship_id': None if query_dict.get('starship_id') == 'NaN' else int(query_dict.get('starship_id'))
        }

    def get_order_by_arg(self, order_col, order_dir):
        if order_col == 2:
            return 'accession_sql' if order_dir == 'asc' else '-accession_sql'
        elif order_col == 3:
            return 'annotation' if order_dir == 'asc' else '-annotation'
        elif order_col == 4:
            return Length('sequence').asc() if order_dir == 'asc' else Length('sequence').desc()
        elif order_col == 5:
            return 'public_notes' if order_dir == 'asc' else '-public_notes'
        elif order_col == 6:
            return 'private_notes' if order_dir == 'asc' else '-private_notes'
        elif order_col == 7:
            return 'flag_str' if order_dir == 'asc' else '-flag_str'
        elif order_col == 8:
            return 'assigned_to__username' if order_dir == 'asc' else '-assigned_to__username'
        else:
            return 'accession_sql' if order_dir == 'asc' else '-accession_sql'
