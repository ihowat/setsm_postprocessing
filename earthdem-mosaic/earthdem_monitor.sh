#!/bin/bash

usage() {
    echo "Aggregate slurm jobs by JobName and State."
    echo
    echo "Usage: $0 [--starttime YYYY-MM-DD]"
    echo
    echo "    --starttime YYYY-MM-DD    Show jobs started on the provided date."
    echo "    --zone TEXT               Show results only for the specified zone."
    exit 0
}

# Handle options
while (( "$#" )); do
    case "$1" in
        -h|--help)
            usage
            ;;
        --starttime)
            starttime_provided=true
            starttime=$2
            shift
            shift
            ;;
        --zone)
            zone_provided=true
            zone=$2
            shift
            shift
            ;;
        --) # end argument parsing
            shift
            break
            ;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1
            ;;
    esac
done

temp_csv="~/.earthdem_monitor.csv"

if [[ $zone_provided ]]; then
        where_clause=" where utm_zone = '${zone}' "
else
        where_clause=""
fi


SQL="PIVOT (select regexp_extract(JobName, '(.+)_(utm\d\d[n|s])', 2) as utm_zone, regexp_extract(JobName, '(.+)_(utm\d\d[n|s])', 1) as job_prefix, State from read_csv('${temp_csv}')${where_clause}) on State using count(job_prefix) group by utm_zone, job_prefix order by utm_zone, job_prefix"

sacct_cmd="sacct --user=$USER --format=JobName,State --allocations --parsable2 --delimiter=','"

duckdb_cmd="duckdb :memory: \"${SQL}\""

if [[ $starttime_provided ]]; then

    sacct_cmd="${sacct_cmd} --starttime=${starttime}"

fi

echo "Job states at: $(date)"
eval "${sacct_cmd} > ${temp_csv}"
eval "${duckdb_cmd}"