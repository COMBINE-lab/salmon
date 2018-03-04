ALVBIN=$1
OUT=$PWD

tfile=$(mktemp /tmp/foo.XXXXXXXXX)

/usr/bin/time -o $tfile $ALVBIN alevin -lISR  -1 /mnt/scratch5/avi/alevin/data/10x/mohu/100/all_bcs.fq -2 /mnt/scratch5/avi/alevin/data/10x/mohu/100/all_reads.fq --chromium -o $OUT/prediction -i /mnt/scratch5/avi/alevin/testing/salmonData/index/ -p 20 --dedup --tgMap /mnt/scratch5/avi/alevin/data/mohu/gtf/txp2gene.tsv --mrna /mnt/scratch5/avi/alevin/data/mohu/gtf/mrna.txt --rrna /mnt/scratch5/avi/alevin/data/mohu/gtf/rrna.txt  --dumpbarcodeeq --dumpfeatures --dumpbarcodemap

cat $tfile

tar -xvzf alevin_test_data.tar.gz

echo "Barcodes.txt"
sort prediction/alevin/barcodes.txt > 1.txt
sort alevin_test_data/alevin/barcodes.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "Above line should be Zero"

echo "BarcodeSoftMap.txt"
sort prediction/alevin/barcodeSoftMaps.txt > 1.txt
sort alevin_test_data/alevin/barcodeSoftMaps.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "Above line should be Zero"

echo "frequency.txt"
sort prediction/alevin/frequency.txt > 1.txt
sort alevin_test_data/alevin/frequency.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "Above line should be Zero"

echo "mappedUMI.txt"
sort prediction/alevin/MappedUmi.txt > 1.txt
sort alevin_test_data/alevin/MappedUmi.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "Above line should be Zero"

rm 1.txt 2.txt diff.txt

#python2 alevin_test_data/src-py/alevin.py --txps 322667 --one prediction/alevin/quants_mat.csv --b1 prediction/alevin/barcodes.txt --two alevin_test_data/alevin/quants_mat.csv --b2 alevin_test_data/alevin/barcodes.txt  --type csv

python2 alevin_test_data/src-py/alevin.py --txps 322667 --one prediction/alevin/quants_mat.gz --b1 prediction/alevin/barcodes.txt --two alevin_test_data/alevin/quants_mat.gz --b2 alevin_test_data/alevin/barcodes.txt --type sf

#python2 alevin_test_data/src-py/alevin.py --txps 322667 --one prediction/alevin/cell_eq_mat.gz --b1 prediction/alevin/barcodes.txt --two alevin_test_data/alevin/cell_eq_mat.gz --b2 alevin_test_data/alevin/barcodes.txt --type eq
