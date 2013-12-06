#!/bin/bash
# ####################################################################
# My wrapper for the bamboo launch script. This was needed when we 
# moved to a CentOS host system, where you need use statements in your
# environment.
# ####################################################################
# eval /broad/software/dotkit/init <- don't use this; old approach for the next line
. /broad/tools/scripts/useuse
reuse .combined_LSF_SGE
reuse Ruby-1.8
reuse R-2.11
reuse Oracle-full-client
reuse .cx-oracle-5.0.2-python-2.6.5-oracle-full-client-11.1
reuse Java-1.7
reuse Ant-1.8
reuse Git-1.7
reuse Python-2.7

if [[ !(${LD_LIBRARY_PATH} =~ bwa) ]]
then
	export LD_LIBRARY_PATH=/humgen/gsa-hpprojects/GATK/data/bwa:${LD_LIBRARY_PATH}
fi

echo "LD_LIBRARY_PATH:"
echo $LD_LIBRARY_PATH
echo ""
echo "JAVA_HOME:"
echo $JAVA_HOME
echo ""

# Silence BSUB submit messages so that they don't print
# to stderr and show up as "errors" in Bamboo builds.
export BSUB_QUIET=Y

BAMBOO_INSTALL_DIR="atlassian-bamboo-4.4.1"

echo "trying to run"
if [ -f "/local/gsa-engineering/bamboo/${BAMBOO_INSTALL_DIR}/bamboo.pid" ]
then
	echo "Shutting down previous Bamboo instance..."
	cd "${BAMBOO_INSTALL_DIR}"
	sh ./bamboo.sh stop
	cd ..
fi

git version
javac -version
ant -version

cd "${BAMBOO_INSTALL_DIR}"
sh ./bamboo.sh start
cd ..

