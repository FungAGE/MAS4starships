#!/bin/bash
# Watch for changes in static files and collect them
while true; do
    sudo docker exec mas /home/daemon/miniconda/envs/mas/bin/python manage.py collectstatic --noinput
    sleep 2
done