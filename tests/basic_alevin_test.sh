ALVBIN=$1
OUT=$PWD

tfile=$(mktemp /tmp/foo.XXXXXXXXX)

/usr/bin/time -o $tfile $ALVBIN alevin -lISR --chromium  -1 /mnt/scratch5/avi/alevin/data/10x/v2/mohu/100/all_bcs.fq -2 /mnt/scratch5/avi/alevin/data/10x/v2/mohu/100/all_reads.fq -o $OUT/prediction -i /mnt/scratch5/avi/alevin/data/mohu/salmon_index/ -p 20 --tgMap /mnt/scratch5/avi/alevin/data/mohu/gtf/txp2gene.tsv --numCellBootstraps 100 --dumpFeatures #--dumpBfh --dumpBarcodeEq --dumpBarcodeMap --dumpCsvCounts --chromium --dumpUmiGraph --expectCells 1001  --end 6 --umiLength 10 --barcodeLength 16 

cat $tfile

tar -xvzf alevin_test_data.tar.gz

echo "Barcodes.txt"
sort prediction/alevin/quants_mat_rows.txt > 1.txt
sort alevin_test_data/alevin/quants_mat_rows.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "FAILED if above line > Zero"

echo "BarcodeSoftMap.txt"
sort prediction/alevin/barcodeSoftMaps.txt > 1.txt
sort alevin_test_data/alevin/barcodeSoftMaps.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "FAILED if above line > Zero"

echo "frequency.txt"
sort prediction/alevin/raw_cb_frequency.txt > 1.txt
sort alevin_test_data/alevin/raw_cb_frequency.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "FAILED if above line > Zero"

echo "mappedUMI.txt"
sort prediction/alevin/MappedUmi.txt > 1.txt
sort alevin_test_data/alevin/MappedUmi.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "FAILED if above line > Zero"

echo "whitelist.txt"
sort prediction/alevin/whitelist.txt > 1.txt
sort alevin_test_data/alevin/whitelist.txt > 2.txt
diff 1.txt 2.txt  > diff.txt
wc -l diff.txt
echo "FAILED if above line > Zero"

rm 1.txt 2.txt diff.txt

#python2 alevin_test_data/src-py/alevin.py --txps 322667 --one prediction/alevin/quants_mat.csv --b1 prediction/alevin/barcodes.txt --two alevin_test_data/alevin/quants_mat.csv --b2 alevin_test_data/alevin/barcodes.txt  --type csv

python alevin_test_data/src-py/alevin.py --txps 107450  --one prediction/alevin/quants_mat.gz --b1 prediction/alevin/quants_mat_rows.txt --two alevin_test_data/alevin/quants_mat.gz --b2 alevin_test_data/alevin/quants_mat_rows.txt --type sf

python alevin_test_data/src-py/alevin.py --txps 107450 --one prediction/alevin/cell_eq_mat.gz --b1 prediction/alevin/cell_eq_order.txt --two alevin_test_data/alevin/cell_eq_mat.gz --b2 alevin_test_data/alevin/cell_eq_order.txt --type eq

#python2 alevin_test_data/src-py/alevin.py --txps 322667 --one prediction/alevin/cell_eq_mat.gz --b1 prediction/alevin/barcodes.txt --two alevin_test_data/alevin/cell_eq_mat.gz --b2 alevin_test_data/alevin/barcodes.txt --type eq
