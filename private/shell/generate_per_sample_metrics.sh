#!/bin/sh

STING_HOME=${STING_HOME:-`pwd`/../..}
JYTHON=/humgen/gsa-hpprojects/software/bin/jython2.5.2/jython
SAM_JAR=`ls $STING_HOME/settings/repository/net.sf/sam-*.jar`
PICARD_PUBLIC_JAR=`ls $STING_HOME/settings/repository/net.sf/picard-*.jar`
PICARD_PRIVATE_JAR=`ls $STING_HOME/settings/repository/edu.mit.broad/picard-private-parts-*.jar`
ORACLE_JAR=${ORACLE_HOME:-/broad/software/nonfree/Linux/redhat_5_x86_64/pkgs/oracle_full_client/111_client}/jdbc/lib/ojdbc6.jar

TSV=`basename $1`
METRICS="${TSV%.*}.metrics"
PDF="${TSV%.*}.pdf"

$JYTHON \
  -J-classpath $SAM_JAR:$PICARD_PUBLIC_JAR:$PICARD_PRIVATE_JAR:$ORACLE_JAR \
  $STING_HOME/private/python/generate_per_sample_metrics.py $TSV true > $METRICS && \

Rscript $STING_HOME/private/R/exomePreQC.R $METRICS $PDF
