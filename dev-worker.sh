#!/bin/bash

# Get directory of script
wd="$(dirname "$0")"
cd "$wd" || exit

# Load environment variables
source .env

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate mas-env

# Set development mode
export DEVELOPER_MODE=TRUE
export CELERY_WORKER=TRUE
export DB_HOST=127.0.0.1
export DB_PORT=3307
# Use the same broker URL as in development settings
export CELERY_BROKER_URL="amqp://mas:${RABBITMQ_DEFAULT_PASS}@127.0.0.1:${RABBITMQ_PORT}/"
# Add Luigi password for development
export LUIGI_USER_PASSWORD="changeme"

# Print debug information
echo "Starting Celery worker with:"
echo "CELERY_BROKER_URL: $CELERY_BROKER_URL"
echo "DB_HOST: $DB_HOST"
echo "DB_PORT: $DB_PORT"

# Start celery worker with specific name to match settings and more verbose logging
celery -A MAS worker --hostname=mas-worker@host -l DEBUG --concurrency=1 