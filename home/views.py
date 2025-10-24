from datetime import datetime, timedelta

from django.views import generic

from result_viewer.views import MixinForBaseTemplate
from starship.models import Annotation, StagingStarship
from starship.starbase_models import JoinedShips, ShipQualityTags


def get_starship_breakdown():
    """
    Get breakdown statistics for starships similar to get_database_stats().
    Returns a dictionary with various counts and statistics.
    """
    accepted_quality_tags = ["missing_direct_repeats","missing_tir","missing_boundaries","missing_genome_context","unannotated","missing_empty_site"]
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
        
        # Get quality tag breakdown
        quality_tag_breakdown = {}
        quality_tag_types = ShipQualityTags.objects.values_list('tag_type', flat=True).distinct()
        for tag_type in quality_tag_types:
            if tag_type and tag_type in accepted_quality_tags:
                quality_tag_breakdown[tag_type] = ShipQualityTags.objects.filter(tag_type=tag_type).count()
                
        # Get total quality tags count
        total_quality_tags = ShipQualityTags.objects.count()
        
        # Get ships with quality tags count
        ships_with_quality_tags = JoinedShips.objects.filter(quality_tags__isnull=False).distinct().count()
        
        return {
            'total_starships': total_starships,
            'curated_starships': curated_starships,
            'uncurated_starships': uncurated_starships,
            'unique_accessions': unique_accessions,
            'species_count': species_count,
            'family_count': family_count,
            'evidence_breakdown': evidence_breakdown,
            'source_breakdown': source_breakdown,
            'quality_tag_breakdown': quality_tag_breakdown,
            'total_quality_tags': total_quality_tags,
            'ships_with_quality_tags': ships_with_quality_tags,
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
            'quality_tag_breakdown': {},
            'total_quality_tags': 0,
            'ships_with_quality_tags': 0,
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
