#!/bin/bash

# Script to completely reset the development environment
# WARNING: This will remove ALL media files and reset task states

# Get directory of script and project root
SCRIPT_DIR="$(dirname "$0")"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT" || exit

echo "========================================="
echo "DEVELOPMENT ENVIRONMENT RESET"
echo "========================================="
echo ""
echo "This will:"
echo "1. Stop all services"
echo "2. Remove ALL media files (not just old ones)"
echo "3. Clean up stale running tasks in database"
echo "4. Clean up PID files"
echo ""

read -p "Are you sure you want to continue? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Reset cancelled."
    exit 0
fi

echo ""
echo "Starting reset process..."

# Stop services
echo "1. Stopping services..."
"$SCRIPT_DIR/manage-services.sh" stop

# Remove ALL media files
echo "2. Removing ALL media files..."
if [ -d "$PROJECT_ROOT/media" ]; then
    rm -rf "$PROJECT_ROOT/media"
    echo "   Media directory removed."
else
    echo "   No media directory found."
fi

# Clean up database (reset running tasks)
echo "3. Cleaning up database..."
if [ -f "$PROJECT_ROOT/.env" ]; then
    source "$PROJECT_ROOT/.env"
    
    # Activate conda environment
    eval "$(conda shell.bash hook)"
    conda activate mas-env
    
    # Set environment variables
    export DEVELOPER_MODE=TRUE
    export DB_HOST=127.0.0.1
    export DB_PORT=3307
    
    # Run database cleanup
    python manage.py shell << 'EOF'
from search_manager.models import SearchTask
from django.utils import timezone

# Reset all running tasks to failed status
running_tasks = SearchTask.objects.filter(status=1)
count = running_tasks.count()
if count > 0:
    running_tasks.update(status=2, date_finished=timezone.now())
    print(f"Reset {count} running tasks to failed status.")
else:
    print("No running tasks found to reset.")
EOF
    
    echo "   Database cleanup completed."
else
    echo "   Warning: .env file not found, skipping database cleanup."
fi

# Clean up PID files
echo "4. Cleaning up PID files..."
"$SCRIPT_DIR/manage-services.sh" cleanup

echo ""
echo "========================================="
echo "RESET COMPLETED!"
echo "========================================="
echo ""
echo "Your development environment has been reset."
echo "You can now start fresh with:"
echo "  ./scripts/development/manage-services.sh start"
echo "" 