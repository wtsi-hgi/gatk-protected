#!/bin/sh

SAMPLE_SET=$1
JOB_QUEUE=week
RUN_DIR=`dirname $0`
TMP_DIR=$RUN_DIR/tmp
RETRY=2
SCATTER_COUNT=400
SCATTER_MEMORY=2

shift

BSUB=''
RUN=''
if [ "$1" == "-run" ]; then
    BSUB="bsub -N -P Q -q week -R rusage[mem=2] -cwd $RUN_DIR "
    RUN=-run
    shift
fi

cd $RUN_DIR

$BSUB \
java -Xmx1g -Djava.io.tmpdir=$TMP_DIR -jar Queue.jar \
  -jobQueue $JOB_QUEUE -Y $RUN_DIR/$SAMPLE_SET.tsv \
  -S $RUN_DIR/HybridSelectionPipeline.scala -log $RUN_DIR/queue_log.txt \
  -retry $RETRY -varScatter $SCATTER_COUNT -varMemory $SCATTER_MEMORY \
  -bsub $RUN $@
