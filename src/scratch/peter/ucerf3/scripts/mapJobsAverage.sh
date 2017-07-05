#!/bin/bash

JOBGROUP=UC33avgSolTest

# Local config for script
SHA_LOCAL=/Users/pmpowers/projects/OpenSHA
DIST_LOCAL=$SHA_LOCAL/dist
LIB_LOCAL=$SHA_LOCAL/lib
TMP_LOCAL=$SHA_LOCAL/tmp/UC33/mapjobs
SCRIPT=$TMP_LOCAL/$JOBGROUP.pbs

# Remote config for jobs
BASEDIR=/home/scec-00/pmpowers
JAVADIR=$BASEDIR/lib
SRCDIR=$BASEDIR/UC33/src/conv
OUTDIR=$BASEDIR/UC33/maps/$JOBGROUP

# Calc config
SOL_FILE=$SRCDIR/FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip
SOL_COUNT=2
GRID='CA_RELM'
SPACING='0.1'
PERIOD='GM0P00'
HRS=3
NODES=72
QUEUE=nbns

java -cp $DIST_LOCAL/OpenSHA_complete.jar:$LIB_LOCAL/commons-cli-1.2.jar \
	scratch.peter.ucerf3.scripts.MapsFromAverage \
	$QUEUE $NODES $HRS $JAVADIR $SCRIPT $SOL_COUNT \
	$SOL_FILE $GRID $SPACING $PERIOD $OUTDIR

if [[ $? == 0 ]] ; then
	echo 'PBS script is here:'
	echo $SCRIPT
fi
