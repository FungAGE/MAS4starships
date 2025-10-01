#!/bin/bash

# This script is used to set up the development environment

# Exit on any error
set -e

# Get directory of script and project root
SCRIPT_DIR="$(dirname "$0")"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT" || exit

# Function to check if Docker is running
check_docker() {
    if ! sudo docker info >/dev/null 2>&1; then
        echo "Error: Docker is not running or not accessible!"
        echo "Please start Docker and try again."
        exit 1
    fi
}

# Function to check if conda is available
check_conda() {
    if ! command -v conda &> /dev/null; then
        echo "Error: conda is not installed or not in PATH!"
        echo "Please install conda and try again."
        exit 1
    fi
}

# Function to wait for database to be ready
wait_for_db() {
    echo "Waiting for database to be ready..."
    for i in {1..30}; do
        if sudo docker exec mas-sql-server mysqladmin ping -h"127.0.0.1" --silent; then
            echo "Database is ready!"
            return 0
        fi
        echo "Waiting for database... ($i/30)"
        sleep 2
    done
    echo "Error: Database did not become ready in time"
    exit 1
}

# Check prerequisites
check_docker
check_conda

# Load environment variables
if [ -f "$PROJECT_ROOT/.env" ]; then
    echo "Loading environment variables..."
    set -a
    source "$PROJECT_ROOT/.env"
    set +a
    # Override database host for local development
    export DB_HOST=127.0.0.1
    export DB_PORT=3307
else
    echo "Error: .env file not found!"
    exit 1
fi

# Start required Docker services
echo "Starting required Docker services..."
echo "You may be prompted for your sudo password..."
sudo docker compose up -d messagebroker db

# Wait for database to be ready
wait_for_db

# Check if conda environment exists and create if it doesn't
if ! conda env list | grep -q "^mas-env "; then
    echo "Creating conda environment 'mas-env'..."
    conda env create -f container_setup_files/mas-dev-environment.yml
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate mas-env

# Create logs directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/logs"

# Stop any existing services and clean up
echo "Stopping any existing services..."
"$SCRIPT_DIR/manage-services.sh" stop

# Start services fresh
echo "Starting services..."
"$SCRIPT_DIR/manage-services.sh" start

# Verify services are running
echo "Verifying services..."
"$SCRIPT_DIR/manage-services.sh" status

# Start development server
echo "Starting development server..."
echo "You can access the application at http://localhost:8000"
echo "Use Ctrl+C to stop the server"
echo ""
echo "Available management commands:"
echo "  ./scripts/development/manage-services.sh {start|stop|restart|status|cleanup|media-cleanup}"
echo "  ./scripts/development/cleanup-media.sh  - Clean up old media files"
echo "  ./scripts/development/reset-dev-env.sh  - Complete environment reset"
echo ""
echo "When done, run './scripts/development/manage-services.sh stop' to stop all services"
python manage.py runserver 0.0.0.0:8000

# When the server is stopped with Ctrl+C, stop the services
echo "Stopping services..."
"$SCRIPT_DIR/manage-services.sh" stop 