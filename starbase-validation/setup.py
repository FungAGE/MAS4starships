"""Setup for starbase-validation package."""

from setuptools import setup, find_packages

setup(
    name="starbase-validation",
    version="0.1.0",
    description="Validation/ETL pipeline for Starbase database quality",
    packages=find_packages(),
    install_requires=[
        "PyYAML>=5.0",
    ],
    entry_points={
        "console_scripts": [
            "starbase-validate=starbase_validation.pipeline:main",
        ],
    },
    python_requires=">=3.8",
)
