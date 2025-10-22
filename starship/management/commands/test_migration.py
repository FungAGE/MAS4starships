from django.core.management.base import BaseCommand
from starship.models import StagingStarship, JoinedShips, Accessions, Ships
from django.contrib.auth.models import User

class Command(BaseCommand):
    help = 'Test the staging submission migration process'

    def handle(self, *args, **options):
        # Get the first staging submission
        try:
            submission = StagingStarship.objects.filter(status='pending').first()
            if not submission:
                self.stdout.write(self.style.WARNING('No pending submissions found'))
                return
            
            self.stdout.write(f'Testing migration for submission: {submission.starshipID}')
            
            # Test the migration
            main_starship = submission.migrate_to_main_database()
            self.stdout.write(self.style.SUCCESS(f'Migration successful! Created starship {main_starship.id}'))
            
            # Check if it was created in the main database
            main_count = JoinedShips.objects.filter(starshipID=submission.starshipID).count()
            self.stdout.write(f'Found {main_count} matching starships in main database')
            
        except Exception as e:
            self.stdout.write(self.style.ERROR(f'Migration failed: {str(e)}'))
            import traceback
            self.stdout.write(traceback.format_exc())
