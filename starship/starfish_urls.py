# Starfish-specific URLs without the starfish/ prefix
from django.urls import path
from starship import views

app_name = 'starfish'

urlpatterns = [
    # Starfish run management
    path('', views.StarfishRunListView.as_view(), name='starfish_run_list'),
    path('create/', views.StarfishRunCreateView.as_view(), name='starfish_run_create'),
    path('<int:pk>/', views.StarfishRunDetailView.as_view(), name='starfish_run_detail'),
    path('<int:pk>/edit/', views.StarfishRunUpdateView.as_view(), name='starfish_run_update'),
    path('<int:pk>/delete/', views.StarfishRunDeleteView.as_view(), name='starfish_run_delete'),
    path('<int:pk>/start/', views.StarfishRunStartView.as_view(), name='starfish_run_start'),
    path('<int:pk>/cancel/', views.StarfishRunCancelView.as_view(), name='starfish_run_cancel'),
    path('<int:pk>/rerun/', views.StarfishRunRerunView.as_view(), name='starfish_run_rerun'),
    path('<int:pk>/resume/', views.StarfishRunResumeView.as_view(), name='starfish_run_resume'),
    path('<int:pk>/import/', views.StarfishImportToMasView.as_view(), name='starfish_import_to_mas'),
    path('<int:pk>/status/', views.StarfishRunStatusView.as_view(), name='starfish_run_status'),
    path('<int:pk>/log/', views.StarfishRunLogView.as_view(), name='starfish_run_log'),
    
    # Starfish genome management
    path('<int:run_id>/genomes/', views.StarfishGenomeListView.as_view(), name='starfish_genome_list'),
    path('<int:run_id>/genomes/add/', views.StarfishGenomeCreateView.as_view(), name='starfish_genome_create'),
    path('<int:run_id>/samplesheet/', views.StarfishSamplesheetUploadView.as_view(), name='starfish_samplesheet_upload'),
    
    # Starfish elements
    path('<int:run_id>/elements/', views.StarfishElementListView.as_view(), name='starfish_element_list'),
    path('elements/<int:pk>/', views.StarfishElementDetailView.as_view(), name='starfish_element_detail'),
]
