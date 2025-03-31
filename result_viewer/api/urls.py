from django.urls import include, path

from result_viewer.api.views import *


urlpatterns = [
    path('run-search/', RunSearchAjaxView.as_view(), name='run_search'),
    path('run-search-for-starship/', RunAllStarshipProteinsAjaxView.as_view(), name='run_search_for_starship'),
    path('upload-results/', UploadResultsView.as_view(), name='upload_results'),
    path('test', TestConnectionView.as_view(), name='test_connection'),
    path('get-protein/<str:accession>', GetProtSeqView.as_view(), name='get_protein'),
    path('get-starship/<str:starship_name>', GetStarshipView.as_view(), name='get_starship'),
    path('get-starship-data', GetStarshipDataView.as_view(), name='get_starship_data'),
    path('get-annotation-data', GetAnnotationListView.as_view(), name='get_annotation_data'),
    path('api-auth/', include('rest_framework.urls', namespace='rest_framework'))
]