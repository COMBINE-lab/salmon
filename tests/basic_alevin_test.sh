ALVBIN=$1
OUT=$PWD

tfile=$(mktemp /tmp/foo.XXXXXXXXX)

if [ $2 = "em" ];
then
/usr/bin/time -o $tfile $ALVBIN alevin -1 /mnt/scratch5/avi/alevin/data/10x/mohu/100/all_bcs.fq -2 /mnt/scratch5/avi/alevin/data/10x/mohu/100/all_reads.fq --chromium --nosoftmap --dedup --dumpbarcodeeq --out $OUT/prediction --index /mnt/scratch5/avi/alevin/testing/salmonData/index -la --whitelist /mnt/scratch5/avi/alevin/data/10x/mohu/100/whitelist.tsv
else
/usr/bin/time -o $tfile $ALVBIN alevin -1 /mnt/scratch5/avi/alevin/data/10x/mohu/100/all_bcs.fq -2 /mnt/scratch5/avi/alevin/data/10x/mohu/100/all_reads.fq --chromium --nosoftmap --dedup --dumpbarcodeeq --out $OUT/prediction --index /mnt/scratch5/avi/alevin/testing/salmonData/index -la --whitelist /mnt/scratch5/avi/alevin/data/10x/mohu/100/whitelist.tsv --noem
fi

cat $tfile

tar -xvzf alevin_test_data.tar.gz

python2 alevin_test_data/src-py/get_correlation.py --sf prediction

echo 'EQCLASS #reads'
awk 'NF>1 {sum+=$NF} END {print sum}' prediction/alevin/cell/AAAGATGAGAAACGAG/cell_eq_classes.txt

echo 'EM #reads'
awk '{sum += $3} END {print sum}' prediction/alevin/cell/AAAGATGAGAAACGAG/quant.sf
