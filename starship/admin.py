from django.contrib import admin

from starship.models import ValidationIssue, ValidationRun
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


class ValidationIssueInline(admin.TabularInline):
    model = ValidationIssue
    extra = 0
    readonly_fields = ("issue_type", "table_name", "record_id", "severity", "details", "resolved")


@admin.register(ValidationRun)
class ValidationRunAdmin(admin.ModelAdmin):
    list_display = ("id", "started_at", "finished_at", "status", "dry_run", "db_version")
    list_filter = ("status", "dry_run")
    readonly_fields = ("started_at", "finished_at", "error_message")
    inlines = [ValidationIssueInline]


# admin.site.register(starship_models.HHSearch_Result)
# admin.site.register(starship_models.Blastp_Result)
# admin.site.register(starship_models.RPSBlast_Result)
