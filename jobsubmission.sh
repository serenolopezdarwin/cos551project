#!/bin/bash
scriptname = $1
sbatch -D `pwd` job.sh --export=scriptname
