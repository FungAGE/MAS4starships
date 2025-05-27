#!/bin/bash

# Get project root directory
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate mas-env

# Start luigid daemon
echo "Starting Luigi scheduler..."
luigid --port 8082 --background --pidfile "$PROJECT_ROOT/luigi.pid" --logdir "$PROJECT_ROOT/logs"