#!/bin/bash

# Script to clean up media files created during development
# These files are search results that get stored temporarily

# Get directory of script and project root
SCRIPT_DIR="$(dirname "$0")"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT" || exit

MEDIA_DIR="$PROJECT_ROOT/media"

echo "Cleaning up development media files..."

if [ -d "$MEDIA_DIR" ]; then
    echo "Found media directory: $MEDIA_DIR"
    
    # Show current size
    if command -v du &> /dev/null; then
        CURRENT_SIZE=$(du -sh "$MEDIA_DIR" 2>/dev/null | cut -f1)
        echo "Current media directory size: $CURRENT_SIZE"
    fi
    
    # Count files before cleanup
    FILE_COUNT=$(find "$MEDIA_DIR" -type f 2>/dev/null | wc -l)
    echo "Files before cleanup: $FILE_COUNT"
    
    # Remove all files older than 7 days
    echo "Removing files older than 7 days..."
    find "$MEDIA_DIR" -type f -mtime +7 -delete 2>/dev/null
    
    # Remove empty directories
    echo "Removing empty directories..."
    find "$MEDIA_DIR" -type d -empty -delete 2>/dev/null
    
    # Count files after cleanup
    FILE_COUNT_AFTER=$(find "$MEDIA_DIR" -type f 2>/dev/null | wc -l)
    echo "Files after cleanup: $FILE_COUNT_AFTER"
    
    # Show new size
    if command -v du &> /dev/null; then
        NEW_SIZE=$(du -sh "$MEDIA_DIR" 2>/dev/null | cut -f1)
        echo "New media directory size: $NEW_SIZE"
    fi
    
    echo "Media cleanup completed!"
else
    echo "No media directory found - nothing to clean up."
fi 