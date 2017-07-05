#!/bin/bash

JOBGROUP=run0
QUEUE=nbns
NODES=32
HOURS=6
PERIODS=GM0P00,GM3P00
BG=INCLUDE

#Local config for script
SHA_LOCAL=/Users/pmpowers/projects/OpenSHA
DIST_LOCAL=$SHA_LOCAL/dist
LIB_LOCAL=$SHA_LOCAL/lib
TMP_LOCAL=$SHA_LOCAL/tmp/UC33/curvejobs
SCRIPT=$TMP_LOCAL/invRuns/$JOBGROUP.pbs

# Remote config for jobs
BASEDIR=/home/scec-00/pmpowers
SRCDIR=$BASEDIR/UC33/src
JAVA_LIB=$BASEDIR/lib
JOBDIR=$BASEDIR/UC33/curvejobs
SITEFILE=$JOBDIR/sites/noPBR.txt
BRANCHFILE=$JOBDIR/branches/tree1440.txt

# Equation set weight runs need to be called by branchID
#OUTDIR=$BASEDIR/UC33/curves/vars/$JOBGROUP
#SOLFILE=$SRCDIR/vars/2013_05_09-ucerf3p3-branch-wt-test_COMPOUND_SOL.zip

#OUTDIR=$BASEDIR/UC33/curves/tree/$JOBGROUP
#FILENAME=/tree/2013_05_10-ucerf3p3-production-10runs_COMPOUND_SOL.zip
#SOLFILE=$SRCDIR$FILENAME

OUTDIR=$BASEDIR/UC33/curves/invRuns/$JOBGROUP
SOLFILE=/home/scec-02/kmilner/ucerf3/inversion_compound_plots/2013_05_10-ucerf3p3-production-10runs/run_compounds/2013_05_10-ucerf3p3-production-10runs_run0_COMPOUND_SOL.zip

# build pbs
java -cp $DIST_LOCAL/OpenSHA_complete.jar:$LIB_LOCAL/commons-cli-1.2.jar \
	scratch.peter.ucerf3.scripts.CurvesFromCompound \
	$QUEUE $NODES $HOURS $JAVA_LIB \
	$SCRIPT -exact=1 $SOLFILE $SITEFILE $BRANCHFILE $PERIODS $BG $OUTDIR

if [[ $? == 0 ]] ; then
	echo 'PBS script is here:'
	echo $SCRIPT
fi
	