#!/bin/bash

# Set conda env path
export CONDA_ENVS_PATH="/mnt/sda/johannesson_lab/adrian/bin/conda-envs"

# Activate conda
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)";

# Clean conda cache
conda clean --all -y

# Configure channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Create environment and install dependencies
mamba create -y --name mas-worker python=3.8 mysqlclient;
conda activate mas-worker;
mamba install -c bioconda -c conda-forge -y glimmer blast trnascan-se hhsuite;
pip install luigi celery django django-simple-history django-debug-toolbar django-crispy-forms djangorestframework biopython==1.77 pandas requests==2.24.0 celery;

# Generate a secret key
python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())' > "${0%/*}/MAS/settings_files/secret_key.txt"