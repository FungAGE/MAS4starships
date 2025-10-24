from datetime import datetime, timedelta

from django.views import generic

from result_viewer.views import MixinForBaseTemplate
from starship.models import Annotation, StagingStarship
from starship.starbase_models import JoinedShips, Accessions, Taxonomy, FamilyNames, ShipQualityTags


def get_starship_breakdown():
    """
    Get breakdown statistics for starships similar to get_database_stats().
    Returns a dictionary with various counts and statistics.
    """
    try:
        total_starships = JoinedShips.objects.count()
        curated_starships = JoinedShips.objects.filter(curated_status='curated').count()
        uncurated_starships = total_starships - curated_starships
        unique_accessions = JoinedShips.objects.values('accession_id').distinct().count()
        species_count = JoinedShips.objects.filter(tax_id__isnull=False).values('tax_id').distinct().count()
        family_count = JoinedShips.objects.filter(
            ship_family_id__isnull=False
        ).exclude(
            ship_family_id__familyName__in=['NA', 'None', 'NULL', '']
        ).values('ship_family_id').distinct().count()
        
        evidence_breakdown = {}
        evidence_types = JoinedShips.objects.exclude(evidence__isnull=True).exclude(evidence='').values_list('evidence', flat=True).distinct()
        for evidence_type in evidence_types:
            if evidence_type:
                evidence_breakdown[evidence_type] = JoinedShips.objects.filter(evidence=evidence_type).count()
        
        source_breakdown = {}
        source_types = JoinedShips.objects.exclude(source__isnull=True).exclude(source='').values_list('source', flat=True).distinct()
        for source_type in source_types:
            if source_type:
                source_breakdown[source_type] = JoinedShips.objects.filter(source=source_type).count()
        
        return {
            'total_starships': total_starships,
            'curated_starships': curated_starships,
            'uncurated_starships': uncurated_starships,
            'unique_accessions': unique_accessions,
            'species_count': species_count,
            'family_count': family_count,
            'evidence_breakdown': evidence_breakdown,
            'source_breakdown': source_breakdown,
        }
    except Exception as e:
        return {
            'total_starships': 0,
            'curated_starships': 0,
            'uncurated_starships': 0,
            'unique_accessions': 0,
            'species_count': 0,
            'family_count': 0,
            'evidence_breakdown': {},
            'source_breakdown': {},
        }


class HomePageView(MixinForBaseTemplate, generic.TemplateView):
    template_name = "home/index.html"

    def get_context_data(self, **kwargs):
        context = super(HomePageView, self).get_context_data(**kwargs)

        if self.request.user.is_authenticated:
            weekstart = datetime.now() - timedelta(
                days=datetime.now().weekday(), hours=datetime.now().hour
            )
            this_weeks_annotations = Annotation.history.filter(
                history_date__gte=weekstart,
                history_user=self.request.user,
                history_type="~",
            )
            context["annotated_this_week"] = (
                this_weeks_annotations.values("id").distinct().count()
            )

            last_weeks_annotations = Annotation.history.filter(
                history_date__lte=weekstart,
                history_date__gte=weekstart - timedelta(days=7),
                history_user=self.request.user,
                history_type="~",
            )
            context["annotated_last_week"] = (
                last_weeks_annotations.values("id").distinct().count()
            )

            starship_breakdown = get_starship_breakdown()
            context["ship_count"] = starship_breakdown['total_starships']
            context["starship_breakdown"] = starship_breakdown
            
            # Add staging submission counts
            pending_submissions = StagingStarship.objects.filter(status='pending').count()
            context["pending_submissions"] = pending_submissions
            
            # Check if user has any submissions pending review
            user_pending_submissions = StagingStarship.objects.filter(
                status='pending',
                submitted_by=self.request.user
            ).count()
            context["user_pending_submissions"] = user_pending_submissions

        return context
