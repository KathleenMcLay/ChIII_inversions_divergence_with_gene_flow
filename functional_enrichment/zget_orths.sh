# get the Arabidopsis protein file 
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-59/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.pep.all.fa.gz

# blast using docker (must have docker app installed and running) 
docker run -it --rm -v $PWD:/data ncbi/blast \
  bash -c "
    makeblastdb -in /data/Arabidopsis_thaliana.TAIR10.pep.all.fa -dbtype prot -out /data/arabidopsis_db &&
    blastx -query /data/annotated_without_contam.fnn \
           -db /data/arabidopsis_db \
           -out /data/blast_results.tsv \
           -evalue 1e-5 \
           -outfmt 6 \
           -num_threads 4"
