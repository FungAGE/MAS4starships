# from django.conf.urls import url, include
# from django.views.decorators.cache import cache_page
from django.urls import path
# from django.views.decorators.csrf import csrf_exempt
from starship import views

app_name = 'starship'

urlpatterns = [
    path('starship/detail/<int:pk>', views.Genome_Detail.as_view(), name='phage_detail'),
    path('starship/list/', views.Genome_List_SS.as_view(), name='phage_list'),
    path('starship/download/fasta/<int:starship_id>', views.starship_download_fasta, name='phage_download_fasta'),
    path('starship/download/deliverables/<int:starship_id>', views.download_deliverables, name='download_deliverables'),
    path('custom-starship/upload/', views.Upload_Custom_Genome.as_view(), name='upload_custom_starship'),
    path('starship/delete/', views.Genome_Delete.as_view(), name='phage_delete'),
    path('starship/delete/confirm/', views.Confirm_Genome_Delete.as_view(), name='confirm_phage_delete'),
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
    path('ajax/phage/get', views.Get_Genome.as_view(), name='get_starship'),
    # ajax url used for amino acid sequence drop down
    path('ajax/annotation/get/aa_sequence/', views.Get_AA_Sequence.as_view(), name='get_aa_sequence'),
    path('ajax/feature/get/nucleotide_sequence/', views.Get_Feature_Sequence.as_view(), name='get_feature_sequence'),
]
