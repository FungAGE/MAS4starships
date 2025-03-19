#!/bin/bash

# * note that this script is not meant for testing celery worker/luigi tasks

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

# Start development server
echo "Starting development server..."
python manage.py runserver localhost:8000 