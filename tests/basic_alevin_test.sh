ALVBIN=$1
OUT=$PWD

tfile=$(mktemp /tmp/foo.XXXXXXXXX)

/usr/bin/time -o $tfile $ALVBIN alevin -1 /mnt/scratch4/mohsen/strt_reads/umi/SRR1547890_umi.fastq.gz -2 /mnt/scratch4/mohsen/strt_reads/seqs/SRR1547890.fastq.gz -la --end 5 --umilength 6 --barcodelength 10 --nobarcode -o $OUT/prediction -i /mnt/scratch5/avi/alevin/strt_testing/salData/salIndex/ -p 20 --nosoftmap --dedup --tgMap /mnt/scratch5/avi/alevin/data/mouse/gtf/txp2gene.tsv --dumpbarcodeeq --dumpcsvcounts

cat $tfile

tar -xvzf alevin_test_data.tar.gz

python2 alevin_test_data/src-py/get_correlation.py --sf prediction

echo 'EM #reads'
awk '{sum += $3} END {print sum}' prediction/alevin/cell/AAA/quant.sf
