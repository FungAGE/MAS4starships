from celery.app.control import Control
from django.conf import settings
from MAS.celery import app
import logging

logger = logging.getLogger(__name__)

def worker_status(request):
    """Check if celery worker is active"""
    control = Control(app)
    active_workers = False
    try:
        # Try multiple inspection methods
        logger.debug("Checking worker status...")
        
        # Check ping first
        ping_response = control.inspect().ping()
        logger.debug(f"Ping response: {ping_response}")
        
        # Check active workers
        active = control.inspect().active()
        logger.debug(f"Active workers: {active}")
        
        # Check registered workers
        registered = control.inspect().registered()
        logger.debug(f"Registered workers: {registered}")
        
        # Consider worker active if any of these checks pass
        if ping_response or active or registered:
            active_workers = True
            logger.info("Celery worker is active")
        else:
            logger.warning("No active Celery workers detected")
            
    except Exception as e:
        logger.error(f"Error checking worker status: {str(e)}", exc_info=True)
        active_workers = False
    
    return {
        'worker_active': active_workers
    } 