#!/bin/bash

IMAGE=peridigm-cuda

docker run --rm -i -t --gpus all \
       -u=${UID}\
       -v ${PWD}:${PWD}\
       -w ${PWD}\
       -it ${IMAGE}\
       /bin/bash
