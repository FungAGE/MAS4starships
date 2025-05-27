import luigi
import os

from MAS.settings import BASE_DIR

class Globals(luigi.Config):
    '''
    Global variables. Set using luigi configuration file.
    '''
    # Directory to write output to
    OUTPUT_DIR = luigi.Parameter(
        default=os.path.join(os.path.split(os.path.dirname(os.path.realpath(__file__)))[0], 'output')
    )
    MAS_USERNAME = luigi.Parameter(default='luigi')
    MAS_PASSWORD = luigi.Parameter(default=os.environ.get('LUIGI_USER_PASSWORD', 'changeme'))
    MAS_CRT = luigi.Parameter(default=None)
    ERROR_LOG = os.path.join(BASE_DIR, 'logs', 'error.log')  # Instead of just 'logs'
    CLUSTER = luigi.Parameter(default=False)
    NUM_WORKERS = luigi.IntParameter(default=20)
    # CONDA_ENVIRONMENT = luigi.Parameter(default='mas-worker')

global_config = Globals()