from celery.app.control import Control
from django.conf import settings
from MAS.celery import app

def worker_status(request):
    """Check if celery worker is active"""
    control = Control(app)
    active_workers = False
    try:
        if control.inspect().active():
            active_workers = True
    except Exception:
        active_workers = False
    
    return {
        'worker_active': active_workers
    } 