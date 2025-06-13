#!/bin/bash

# load CONDA_ROOT
. .env

# load module for nextflow
. $CONDA_ROOT/etc/profile.d/conda.sh

conda activate nf-core