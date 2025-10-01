# from django.conf.urls import url, include
# from django.views.decorators.cache import cache_page
from django.urls import path
# from django.views.decorators.csrf import csrf_exempt
from starship import views

app_name = 'starship'

urlpatterns = [
    path('starship/detail/<int:pk>', views.Starship_Detail.as_view(), name='starship_detail'),
    path('starship/list/', views.Starship_List.as_view(), name='starship_list'),
    path('starship/download/fasta/<int:starship_id>', views.starship_download_fasta, name='starship_download_fasta'),
    path('starship/download/deliverables/<int:starship_id>', views.download_deliverables, name='download_deliverables'),
    path('starship/upload/', views.StarshipUpload.as_view(), name='upload_starship'),
    path('starship/<int:pk>/delete/', views.Starship_Delete.as_view(), name='starship_delete'),
    path('starship/delete/confirm/', views.Confirm_Starship_Delete.as_view(), name='confirm_starship_delete'),
    path('feature/detail/<int:pk>', views.Feature_Detail.as_view(), name='feature_detail'),
    path('annotation/list/', views.Annotation_List_Serverside.as_view(), name='annotation_list'),
    path('annotation/detail/<int:pk>', views.Annotation_Detail.as_view(), name='annotation_detail'),
    path('annotation/download/<int:annotation_id>', views.annotation_download, name='annotation_download'),
    path('annotation/upload/', views.Upload_Annotation.as_view(), name='upload_annotations'),
    path('annotation/upload/confirm/', views.Confirm_Upload_Annotation.as_view(), name='confirm_upload_annotations'),
    path('annotation/history/<int:annotation_pk>', views.Annotation_History.as_view(), name='annotation_history'),
    path('annotation/download/excel/', views.download_excel_annotations, name='download_excel_annotations'),
    path('annotation/download/excel_template/', views.download_excel_template, name='download_excel_template'),
    path('annotation/download/excel/unannotated_annotations/', views.download_unannotated_annotations, name='download_unannotated_annotations'),

    # ajax url used for nucleotide sequence drop down
    path('ajax/starship/get', views.Get_Starship.as_view(), name='get_starship'),
    # ajax url used for amino acid sequence drop down
    path('ajax/annotation/get/aa_sequence/', views.Get_AA_Sequence.as_view(), name='get_aa_sequence'),
    path('ajax/feature/get/nucleotide_sequence/', views.Get_Feature_Sequence.as_view(), name='get_feature_sequence'),
    
    # URLs for comprehensive data models from SQLite schema
    path('accessions/', views.AccessionsList.as_view(), name='accessions_list'),
    path('accessions/<int:pk>/', views.AccessionDetail.as_view(), name='accession_detail'),
    path('accessions/create/', views.AccessionCreate.as_view(), name='accession_create'),
    path('joined-ships/', views.JoinedShipsList.as_view(), name='joined_ships_list'),
    path('joined-ships/<int:pk>/', views.JoinedShipDetail.as_view(), name='joined_ship_detail'),
    path('taxonomy/', views.TaxonomyList.as_view(), name='taxonomy_list'),
    path('papers/', views.PapersList.as_view(), name='papers_list'),
    path('bulk-upload/', views.BulkDataUpload.as_view(), name='bulk_data_upload'),
    path('quality-overview/', views.StarshipQualityOverview.as_view(), name='quality_overview'),
]
