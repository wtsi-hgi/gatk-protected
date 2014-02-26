#!/bin/bash
#
# Downloads and parses GATK reports from our S3 buckets. Processes reports asynchronously
# in fixed-size chunks.
#

timestamp=`date +"%Y_%m_%d_%H:%M"`
s3_root_dir="/humgen/gsa-hpprojects/GATK/reports/s3"
s3_config_dir="${s3_root_dir}/config"
archive_dir="${s3_root_dir}/archive"
download_root="/local/gsa-engineering/GATKLogs"
# script_dir="/local/gsa-engineering/cron_clones/unstable/private/GATKLogs"
script_dir="${download_root}/scripts"
download_chunk_size=100000
global_timeout_in_seconds=80000
s3funnel_dir="${download_root}/s3funnel_mod"
s3funnel_threads=10
gsa_bucket="broad.gsa.gatk.run.reports"
mark_bucket="GATK_Run_Reports"

export PYTHONPATH="${s3funnel_dir}:${PYTHONPATH}"

log_event() {
    echo `date` " -- $1"
}

total_chunks=0

for account in gsa mark
do
    if [ $account == "gsa" ]
    then
        s3_bucket="${gsa_bucket}"
    else
        s3_bucket="${mark_bucket}"
    fi

    s3_config_file="${s3_config_dir}/s3cfg_${account}"
    download_master_id="${timestamp}_${account}"
    bucket_list_file="${download_root}/${download_master_id}_ls"

    export AWS_ACCESS_KEY_ID=`grep access_key "${s3_config_file}" | awk -F' = ' '{ print $2; }'`
    export AWS_SECRET_ACCESS_KEY=`grep secret_key "${s3_config_file}" | awk -F' = ' '{ print $2; }'`

    log_event "Retrieving file listing for account ${account}"
    "${s3funnel_dir}/s3funnel" --insecure -t "${s3funnel_threads}" "${s3_bucket}" LIST > "${bucket_list_file}"
    log_event "Done retrieving file listing for account ${account}"

    cd "${download_root}"
    rm -f *.done
    chunk_prefix="${download_master_id}_part_"
    split -d -a 3 -l "${download_chunk_size}" "${bucket_list_file}" "${chunk_prefix}"

    num_chunks=`ls -1 ${chunk_prefix}* | wc -l`
    total_chunks=`expr ${total_chunks} + ${num_chunks}`
    log_event "Processing ${num_chunks} chunks for account ${account}"

    for chunk in `ls -1 ${chunk_prefix}*`
    do
        log_event "Dispatching asynchronous downloader process for chunk ${chunk}"
        nohup "${script_dir}/downloadGATKReportsChunk.sh" "${chunk}" "${s3_bucket}" "${s3_config_file}" > "${chunk}.log" 2>&1 &
        sleep 1
    done

    rm "${bucket_list_file}"
done

num_done_chunks=0
elapsed_seconds=0
while [ ${num_done_chunks} -lt ${total_chunks} ]
do
    sleep 60
    elapsed_seconds=`expr ${elapsed_seconds} + 60`

    if [ ${elapsed_seconds} -gt ${global_timeout_in_seconds} ]
    then
        log_event "Global timeout of ${global_timeout_in_seconds} reached, giving up waiting for chunk downloaders to complete."
        exit 1
    fi

    num_done_chunks=`ls -1 ${timestamp}*.done | wc -l`
    log_event "${num_done_chunks} chunks done"
done

rm -f ${timestamp}*.done
log_event "${num_done_chunks} chunks finished processing. Echoing log output:"

for chunk_log in ${timestamp}*.log
do
    echo ""
    echo "${chunk_log}"
    echo "---------------------------------------------------------------"
    cat "${chunk_log}"
done

log_event "Archiving logs"
mv ${timestamp}*.log "${archive_dir}"

exit 0
