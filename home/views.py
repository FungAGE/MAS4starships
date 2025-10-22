from datetime import datetime, timedelta

from django.views import generic

from result_viewer.views import MixinForBaseTemplate
from starship.models import Annotation, StagingStarship
from starship.starbase_models import JoinedShips


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

            ship_count = JoinedShips.objects.count()
            context["ship_count"] = ship_count
            
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
