#!/bin/bash

condition=$1
feature=$2

Rscript --vanilla src/model_training_testing_unimodal.R condition feature

