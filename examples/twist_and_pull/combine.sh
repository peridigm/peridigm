#!/usr/bin/env bash

BASENAME=cylinder
NUMPROCS=4

epu -p ${NUMPROCS} ${BASENAME}QuasiStatic
epu -p ${NUMPROCS} ${BASENAME}Explicit
conjoin  -output ${BASENAME}.e ${BASENAME}QuasiStatic.e ${BASENAME}Explicit.e

