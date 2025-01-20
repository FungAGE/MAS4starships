from django.contrib import admin

# Register your models here.
from django.contrib import admin

from . import models as starship_models


admin.site.register(starship_models.Starship)
admin.site.register(starship_models.Feature)
admin.site.register(starship_models.Annotation)
# admin.site.register(starship_models.HHSearch_Result)
# admin.site.register(starship_models.Blastp_Result)
# admin.site.register(starship_models.RPSBlast_Result)
