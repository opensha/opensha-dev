#!/bin/bash

JOBGROUP=UC33conv
QUEUE=nbns
NODES=64
HOURS=3
PERIODS=GM0P00,GM0P20,GM1P00,GM4P00
BG=INCLUDE

#Local config for script
SHA_LOCAL=/Users/pmpowers/projects/OpenSHA
DIST_LOCAL=$SHA_LOCAL/dist
LIB_LOCAL=$SHA_LOCAL/lib
TMP_LOCAL=$SHA_LOCAL/tmp/UC33/curvejobs
SCRIPT=$TMP_LOCAL/$JOBGROUP.pbs

# Remote config for jobs
BASEDIR=/home/scec-00/pmpowers
SRCDIR=$BASEDIR/UC33/src
JAVA_LIB=$BASEDIR/lib
JOBDIR=$BASEDIR/UC33/curvejobs
SITEFILE=$JOBDIR/sites/all.txt

# Convergence runs need to be called by index
OUTDIR=$BASEDIR/UC33/curves/conv/$JOBGROUP
FILENAME=/conv/FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate7.9_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip
SOLFILE=$SRCDIR$FILENAME
ERF_COUNT=200

# build pbs
java -cp $DIST_LOCAL/OpenSHA_complete.jar:$LIB_LOCAL/commons-cli-1.2.jar \
	scratch.peter.ucerf3.scripts.CurvesFromAverage \
	$QUEUE $NODES $HOURS $JAVA_LIB $SCRIPT \
	-exact=1 $SOLFILE $SITEFILE $ERF_COUNT $PERIODS $BG $OUTDIR

if [[ $? == 0 ]] ; then
	echo 'PBS script is here:'
	echo $SCRIPT
fi
