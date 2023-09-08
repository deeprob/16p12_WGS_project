


# jobs that finished successfully have a file size of 941
# get the job number for all other jobs
ll logs/3_regenotype/ | grep -v "941 F" | grep -v tot | cut -f9 -d' ' | cut -f1 -d. > 3_failed_jobs.txt



