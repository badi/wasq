#!/usr/bin/env bash

rsync -r --progress --stats haldarfe:'/pscratch/cabdulwa/adaptive-sampling-simulations/t_*' .
