#!/bin/bash

# Get directory of script and project root
SCRIPT_DIR="$(dirname "$0")"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT" || exit

# Source environment variables
if [ -f "$PROJECT_ROOT/.env" ]; then
    source "$PROJECT_ROOT/.env"
else
    echo "Error: .env file not found!"
    exit 1
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate mas-env

# Set development environment variables
export DEVELOPER_MODE=TRUE
export CELERY_WORKER=TRUE
export DB_HOST=127.0.0.1
export DB_PORT=3307
export CELERY_BROKER_URL="amqp://mas:${RABBITMQ_DEFAULT_PASS}@127.0.0.1:${RABBITMQ_PORT}/"
export LUIGI_USER_PASSWORD="changeme"

# Function to check if a process is running
check_process() {
    local pattern=$1
    pgrep -f "$pattern" >/dev/null
    return $?
}

# Function to start Luigi
start_luigi() {
    if check_process "luigid.*--port 8082"; then
        echo "Luigi is already running."
    else
        echo "Starting Luigi scheduler..."
        luigid --port 8082 --background --pidfile "$PROJECT_ROOT/luigi.pid" --logdir "$PROJECT_ROOT/logs"
        sleep 2
        if check_process "luigid.*--port 8082"; then
            echo "Luigi started successfully."
        else
            echo "Failed to start Luigi!"
            return 1
        fi
    fi
}

# Function to start Celery
start_celery() {
    if check_process "celery.*mas-worker@host"; then
        echo "Celery worker is already running."
    else
        echo "Starting Celery worker..."
        celery -A MAS worker --hostname=mas-worker@host -l DEBUG --concurrency=1 --pidfile="$PROJECT_ROOT/celery.pid" &
        sleep 2
        if check_process "celery.*mas-worker@host"; then
            echo "Celery worker started successfully."
        else
            echo "Failed to start Celery worker!"
            return 1
        fi
    fi
}

# Function to stop Luigi
stop_luigi() {
    if [ -f "$PROJECT_ROOT/luigi.pid" ]; then
        echo "Stopping Luigi scheduler..."
        kill "$(cat "$PROJECT_ROOT/luigi.pid")" 2>/dev/null
        rm -f "$PROJECT_ROOT/luigi.pid"
    fi
    pkill -f "luigid.*--port 8082"
    echo "Luigi scheduler stopped."
}

# Function to stop Celery
stop_celery() {
    if [ -f "$PROJECT_ROOT/celery.pid" ]; then
        echo "Stopping Celery worker..."
        kill "$(cat "$PROJECT_ROOT/celery.pid")" 2>/dev/null
        rm -f "$PROJECT_ROOT/celery.pid"
    fi
    pkill -f "celery.*mas-worker@host"
    echo "Celery worker stopped."
}

# Function to show status
show_status() {
    echo "Checking service status..."
    
    if check_process "luigid.*--port 8082"; then
        echo "Luigi: RUNNING"
    else
        echo "Luigi: STOPPED"
    fi
    
    if check_process "celery.*mas-worker@host"; then
        echo "Celery: RUNNING"
    else
        echo "Celery: STOPPED"
    fi
}

# Function to clean up stale PID files
cleanup_pids() {
    echo "Cleaning up stale PID files..."
    for pid_file in "$PROJECT_ROOT"/*.pid; do
        if [ -f "$pid_file" ]; then
            if ! kill -0 "$(cat "$pid_file")" 2>/dev/null; then
                rm -f "$pid_file"
                echo "Removed stale PID file: $pid_file"
            fi
        fi
    done
}

# Function to clean up media files
cleanup_media() {
    echo "Running media cleanup..."
    "$SCRIPT_DIR/cleanup-media.sh"
}

# Function to show usage
show_usage() {
    echo "Usage: $0 {start|stop|restart|status|cleanup|media-cleanup}"
    echo ""
    echo "Commands:"
    echo "  start        - Start Luigi and Celery services"
    echo "  stop         - Stop Luigi and Celery services"
    echo "  restart      - Restart Luigi and Celery services"
    echo "  status       - Show status of Luigi and Celery services"
    echo "  cleanup      - Clean up stale PID files"
    echo "  media-cleanup - Clean up old media files from development"
    echo ""
}

# Main command handling
case "$1" in
    start)
        cleanup_pids
        start_luigi
        start_celery
        ;;
    stop)
        stop_celery
        stop_luigi
        cleanup_pids
        ;;
    restart)
        stop_celery
        stop_luigi
        cleanup_pids
        sleep 2
        start_luigi
        start_celery
        ;;
    status)
        show_status
        ;;
    cleanup)
        cleanup_pids
        ;;
    media-cleanup)
        cleanup_media
        ;;
    *)
        show_usage
        exit 1
        ;;
esac

exit 0 