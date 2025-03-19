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
# Override message broker settings for development using the same port as docker-compose
export CELERY_BROKER_URL="amqp://mas:${RABBITMQ_DEFAULT_PASS}@127.0.0.1:${RABBITMQ_PORT}/"

# Start celery worker with debug logging
celery -A MAS worker -l debug 