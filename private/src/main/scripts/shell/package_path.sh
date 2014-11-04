#!/bin/bash
#
# Other scripts that need the paths to package jars should use this script.
#
# Should be run from the root directory of your git clone.
#
# Example:
#   java -jar `private/src/main/scripts/shell/package_path.sh [ gatk | queue ] [ protected | private ]`  ...
#

case "$1$2" in

    gatkprivate)
        ls private/gatk-package-internal/target/gatk-package-internal-*.jar
        ;;

    queueprivate)
        ls private/gatk-queue-package-internal/target/gatk-queue-package-internal-*.jar
        ;;

    gatkprotected)
        ls protected/gatk-package-distribution/target/gatk-package-distribution-*.jar
        ;;

    queueprotected)
        ls protected/gatk-queue-package-distribution/target/gatk-queue-package-distribution-*.jar
        ;;

    *)
        echo "Usage: $0 [ gatk | queue ] [ protected | private ]" >&2
        exit 1

esac
