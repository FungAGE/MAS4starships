from django.contrib import admin

from starship.models import (
    GeneAnnotationRun,
    StagingGffFeature,
    ValidationIssue,
    ValidationRun,
)
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


class StagingGffFeatureInline(admin.TabularInline):
    model = StagingGffFeature
    extra = 0
    fields = (
        "accession_tag",
        "feature_type",
        "start",
        "end",
        "source",
        "review_status",
        "reviewed_by",
        "promoted_at",
    )
    readonly_fields = ("accession_tag", "feature_type", "start", "end", "source", "promoted_at")
    raw_id_fields = ("reviewed_by",)


@admin.register(GeneAnnotationRun)
class GeneAnnotationRunAdmin(admin.ModelAdmin):
    list_display = ("id", "accession_tag", "tool", "status", "created_at", "completed_at")
    list_filter = ("status", "tool")
    search_fields = ("accession_tag",)
    readonly_fields = ("created_at", "completed_at", "celery_task_id")
    inlines = [StagingGffFeatureInline]


@admin.register(StagingGffFeature)
class StagingGffFeatureAdmin(admin.ModelAdmin):
    list_display = ("id", "run", "accession_tag", "feature_type", "start", "end", "review_status")
    list_filter = ("review_status", "source")
    search_fields = ("accession_tag",)
    raw_id_fields = ("run", "reviewed_by")


# admin.site.register(starship_models.HHSearch_Result)
# admin.site.register(starship_models.Blastp_Result)
# admin.site.register(starship_models.RPSBlast_Result)
