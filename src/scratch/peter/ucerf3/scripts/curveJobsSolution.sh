#!/bin/bash

JOBGROUP=UC33curveSolTest
PERIODS=GM0P00,GM0P20,GM1P00

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
SITEFILE=$BASEDIR/UC33/curvejobs/sites/test.txt

OUTDIR=$BASEDIR/UC33/curves/conv/$JOBGROUP
FILENAME=/conv/FM3_1_ZENGBB_Shaw09Mod_DsrTap_CharConst_M5Rate8.7_MMaxOff7.6_NoFix_SpatSeisU3_mean_sol.zip
SOLFILE=$SRCDIR$FILENAME

# build pbs
java -cp $DIST_LOCAL/OpenSHA_complete.jar:$LIB_LOCAL/commons-cli-1.2.jar \
	scratch.peter.ucerf3.scripts.CurvesFromSolution \
	$JAVA_LIB $SCRIPT $SOLFILE $SITEFILE $PERIODS $OUTDIR

if [[ $? == 0 ]] ; then
	echo 'PBS script is here:'
	echo $SCRIPT
fi
	