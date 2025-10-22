from django.contrib import admin

from starship.starbase_models import (
    Accessions, Ships, Captains, Taxonomy, Genome, Papers,
    FamilyNames, StarshipFeatures, Navis, Haplotype, Gff, JoinedShips
)

admin.site.register(Accessions)
admin.site.register(Ships)
admin.site.register(Captains)
admin.site.register(Taxonomy)
admin.site.register(Genome)
admin.site.register(Papers)
admin.site.register(FamilyNames)
admin.site.register(StarshipFeatures)
admin.site.register(Navis)
admin.site.register(Haplotype)
admin.site.register(Gff)
admin.site.register(JoinedShips)

# admin.site.register(starship_models.HHSearch_Result)
# admin.site.register(starship_models.Blastp_Result)
# admin.site.register(starship_models.RPSBlast_Result)
