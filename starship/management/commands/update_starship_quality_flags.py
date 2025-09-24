from django.core.management.base import BaseCommand
from starship.models import JoinedShips


class Command(BaseCommand):
    help = 'Update quality flags for all Starship sequences based on available data'

    def add_arguments(self, parser):
        parser.add_argument(
            '--dry-run',
            action='store_true',
            help='Show what would be updated without making changes',
        )

    def handle(self, *args, **options):
        dry_run = options['dry_run']
        
        if dry_run:
            self.stdout.write(self.style.WARNING('DRY RUN MODE - No changes will be made'))
        
        starships = JoinedShips.objects.all()
        total_count = starships.count()
        
        self.stdout.write(f'Processing {total_count} starships...')
        
        flag_counts = {
            0: 0,  # COMPLETE
            1: 0,  # MISSING_BOUNDARIES
            2: 0,  # MISSING_CAPTAIN
            3: 0,  # MISSING_CLASSIFICATION
            4: 0,  # INCOMPLETE
            5: 0,  # UNANNOTATED
        }
        
        for i, starship in enumerate(starships):
            old_flag = starship.quality_flag
            new_flag = starship.calculate_quality_flag()
            
            flag_counts[new_flag] += 1
            
            if old_flag != new_flag:
                flag_display = dict(JoinedShips.QUALITY_FLAG_CHOICES)[new_flag]
                self.stdout.write(
                    f'  {starship.starshipID}: {old_flag} -> {new_flag} ({flag_display})'
                )
                
                if not dry_run:
                    starship.quality_flag = new_flag
                    starship.save(update_fields=['quality_flag'])
            
            # Progress indicator
            if (i + 1) % 100 == 0:
                self.stdout.write(f'  Processed {i + 1}/{total_count} starships')
        
        # Summary
        self.stdout.write('\nSummary of quality flags:')
        for flag_value, count in flag_counts.items():
            flag_display = dict(JoinedShips.QUALITY_FLAG_CHOICES)[flag_value]
            percentage = (count / total_count * 100) if total_count > 0 else 0
            self.stdout.write(f'  {flag_display}: {count} ({percentage:.1f}%)')
        
        if dry_run:
            self.stdout.write(self.style.WARNING('\nDRY RUN completed. Run without --dry-run to apply changes.'))
        else:
            self.stdout.write(self.style.SUCCESS(f'\nSuccessfully updated quality flags for {total_count} starships'))




