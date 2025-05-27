import os

if os.getenv("DEVELOPER_MODE") == "TRUE":
    from .settings_files.development import *

else:
    from .settings_files.production import *

BASE_DIR = BASE_DIR

# KEEP SECRET
SECRET_KEY = (
    open(os.path.join(BASE_DIR, "MAS", "settings_files", "secret_key.txt"), "r")
    .read()
    .strip()
)

if os.getenv("CELERY_WORKER") == "TRUE":
    DATABASES = {
        "default": {
            "ENGINE": "django.db.backends.mysql",
            "OPTIONS": {
                "host": "0.0.0.0",
                "port": 3307,
                "database": "mas",
                "user": "root",
                "password": os.getenv("MYSQL_ROOT_PASSWORD"),
            },
        }
    }
else:
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.mysql',
            'OPTIONS': {
                'host': os.getenv('DB_HOST', 'mas-sql-server'),
                'database': 'mas',
                'user': 'root',
                'password': os.getenv('MYSQL_ROOT_PASSWORD'),
                'port': int(os.getenv('DB_PORT', '3306')),
            }
        }
    }

ROOT_URLCONF = "MAS.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [os.path.join(BASE_DIR, "templates")],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

WSGI_APPLICATION = "MAS.wsgi.application"

DATABASE_ROUTERS = ["MAS.routing.DatabaseRouter"]

AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]

LANGUAGE_CODE = "en-us"

TIME_ZONE = "Europe/Stockholm"

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.0/howto/static-files/
STATIC_ROOT = os.path.join(BASE_DIR, "static-files")
STATIC_URL = "/static/"
STATICFILES_DIRS = [
    os.path.join(BASE_DIR, "static"),
]

LOGIN_REDIRECT_URL = "home"

MEDIA_URL = "/media/"
MEDIA_ROOT = os.path.join(BASE_DIR, "media")

FILE_UPLOAD_HANDLERS = ["django.core.files.uploadhandler.TemporaryFileUploadHandler"]

CRISPY_TEMPLATE_PACK = "bootstrap3"
CRISPY_FAIL_SILENTLY = not DEBUG

REST_FRAMEWORK = {
    "DEFAULT_PERMISSION_CLASSES": [
        "rest_framework.permissions.IsAuthenticated",
    ]
}

DATA_UPLOAD_MAX_NUMBER_FIELDS = 10000

INTERNAL_DIR = "/mnt/sda/johannesson_lab/adrian/bin/MAS/databases/internal_db"
NUCLEOTIDE_DB_PATH = os.path.join(INTERNAL_DIR, "nucleotide")
NUCLEOTIDE_FASTA_PATH = os.path.join(INTERNAL_DIR, "nucleotide.fa")
NUCLEOTIDE_DATABASE = os.path.join(INTERNAL_DIR, "nucleotide.db")
PROTEIN_DB_PATH = os.path.join(INTERNAL_DIR, "protein")
PROTEIN_FASTA_PATH = os.path.join(INTERNAL_DIR, "protein.fa")
PROTEIN_DATABASE = os.path.join(INTERNAL_DIR, "protein.db")

SWISSPROT_DIR = "/mnt/sda/johannesson_lab/adrian/bin/MAS/databases/swissprot"
SWISSPROT_DB_PATH = os.path.join(SWISSPROT_DIR, 'uniprot_sprot')
SWISSPROT_FASTA_PATH = os.path.join(SWISSPROT_DIR, 'uniprot_sprot.fasta')

CDD_DIR = "/mnt/sda/johannesson_lab/cdd/cdd_delta"
UNICLUST_DIR = "/mnt/sda/johannesson_lab/uniclust"
PDB_DIR = "/mnt/sda/johannesson_lab/pdb"
INTERPRO_DIR = "/mnt/sda/johannesson_lab/adrian/bin/conda-envs/mas/share/InterProScan/data"

GIT_DIR = os.path.join(BASE_DIR, ".git")

if os.getenv('CELERY_WORKER') == 'TRUE':
    CELERY_BROKER_URL = 'amqp://mas:{password}@mas-message-broker:{port}'.format(
        password=os.getenv('RABBITMQ_DEFAULT_PASS'),
        port=os.getenv('RABBITMQ_PORT', '5672')
    )
else:
    CELERY_BROKER_URL = 'amqp://mas:{password}@127.0.0.1:{port}'.format(
        password=os.getenv('RABBITMQ_DEFAULT_PASS'),
        port=os.getenv('RABBITMQ_PORT', '5672')
    )
CELERY_ACCEPT_CONTENT = ['json']
CELERY_TASK_SERIALIZER = 'json'
CELERY_WORKERS = ['mas-worker@host']

LUIGI_CFG = os.path.join(BASE_DIR, "AnnotationToolPipeline", "luigi.cfg")

STARSHIP_NAME_FORMAT = ""
# GENOME_NAME_FORMAT = r"(?:AMD|INT)_[A-Z]_[a-z]+_[0-9A-Z]+_Phi_[0-9]{3}$"
