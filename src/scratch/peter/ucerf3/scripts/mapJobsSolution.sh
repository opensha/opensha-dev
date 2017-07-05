#!/bin/bash

JOBGROUP=UC33-brAvg-test

# Local config for script
SHA_LOCAL=/Users/pmpowers/projects/OpenSHA
DIST_LOCAL=$SHA_LOCAL/dist
LIB_LOCAL=$SHA_LOCAL/lib
TMP_LOCAL=$SHA_LOCAL/tmp/UC33/mapjobs
SCRIPT=$TMP_LOCAL/$JOBGROUP.pbs

# Remote config for jobs
BASEDIR=/home/scec-00/pmpowers
JAVADIR=$BASEDIR/lib
SRCDIR=$BASEDIR/UC33/src/bravg/FM-DM-MS
OUTDIR=$BASEDIR/UC33/maps/$JOBGROUP

# Calc config
SOL_FILE=$SRCDIR/UC33brAvg_FM31_ABM_ELLB.zip
IMR=WUS_FAULT_13
GRID='CA_RELM'
SPACING='0.1'
PERIOD='GM0P00'
HRS=1
NODES=3
QUEUE=nbns
BG=INCLUDE

java -cp $DIST_LOCAL/OpenSHA_complete.jar:$LIB_LOCAL/commons-cli-1.2.jar \
	scratch.peter.ucerf3.scripts.MapsFromSolution \
	$QUEUE $NODES $HRS $JAVADIR $SCRIPT \
	$SOL_FILE $IMR $GRID $SPACING $PERIOD $BG $OUTDIR

if [[ $? == 0 ]] ; then
	echo 'PBS script is here:'
	echo $SCRIPT
fi
	
