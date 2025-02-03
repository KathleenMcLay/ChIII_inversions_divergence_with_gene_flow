#run revigo - submit job to web based tool
# Shell script for programatic access to Revigo. Run it with:
# revigo.sh revigo_GOterms.csv <- assuming this is list of GO terms in a csv file (csv file is $1)

# Submit job to Revigo
jobid=$(curl "http://revigo.irb.hr/StartJob" -X POST --silent --data-urlencode "cutoff=0.7" --data-urlencode "valueType=pvalue" --data-urlencode "speciesTaxon=0" --data-urlencode "measure=SIMREL" --data-urlencode "goList@$1" --header "Content-Type: application/x-www-form-urlencoded" | jq '.jobid')

# Check job status
running=1
while [ $running -ne 0 ]
do
    running=$(curl "http://revigo.irb.hr/QueryJob" -X GET --silent --data-urlencode "jobid=$jobid" --data-urlencode "type=jstatus" | jq '.running')
    sleep 1
done

# Fetch results
curl "http://revigo.irb.hr/QueryJob" -X GET --silent --data-urlencode "jobid=$jobid" --data-urlencode "namespace=1" --data-urlencode "type=table" > revigo_results.tsv
