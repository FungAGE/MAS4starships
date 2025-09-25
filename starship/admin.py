from django.contrib import admin

# Register your models here.
from django.contrib import admin

from . import models as starship_models

admin.site.register(starship_models.Accessions)
admin.site.register(starship_models.Ships)
admin.site.register(starship_models.Captains)
admin.site.register(starship_models.Taxonomy)
admin.site.register(starship_models.Genome)
admin.site.register(starship_models.Papers)
admin.site.register(starship_models.FamilyNames)
admin.site.register(starship_models.StarshipFeatures)
admin.site.register(starship_models.Navis)
admin.site.register(starship_models.Haplotype)
admin.site.register(starship_models.Gff)
admin.site.register(starship_models.JoinedShips)
admin.site.register(starship_models.Feature)
admin.site.register(starship_models.Annotation)

# admin.site.register(starship_models.HHSearch_Result)
# admin.site.register(starship_models.Blastp_Result)
# admin.site.register(starship_models.RPSBlast_Result)
