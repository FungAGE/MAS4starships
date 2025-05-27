#!/bin/bash

# * note that this script is not meant for testing celery worker/luigi tasks

# Get directory of script and project root
SCRIPT_DIR="$(dirname "$0")"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT" || exit

# Then use PROJECT_ROOT for all paths
source "$PROJECT_ROOT/.env"

# Start required Docker services
echo "Starting database container..."
echo "Please enter your sudo password:"
sudo -S docker compose up -d messagebroker db < /dev/tty

# Check if conda environment exists and create if it doesn't
if ! conda env list | grep -q "^mas-env "; then
    echo "Creating conda environment 'mas-env'..."
    conda env create -f container_setup_files/mas-dev-environment.yml
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate mas-env

# Load environment variables from .env file
if [ -f .env ]; then
    echo "Loading environment variables..."
    set -a
    source .env
    set +a
    # Override database host for local development
    export DB_HOST=127.0.0.1
    export DB_PORT=3307
else
    echo "Error: .env file not found!"
    exit 1
fi

# Run database migrations
# echo "Running database migrations..."
# python manage.py migrate

# Create logs directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/logs"

# Start Luigi scheduler in the background
echo "Starting Luigi scheduler..."
$SCRIPT_DIR/dev-luigi.sh &

# Start Celery worker in the background
echo "Starting Celery worker..."
$SCRIPT_DIR/dev-worker.sh &

# Start development server
echo "Starting development server..."
python manage.py runserver localhost:8000 