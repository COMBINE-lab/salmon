ALVBIN=$1
OUT=$PWD

tfile=$(mktemp /tmp/foo.XXXXXXXXX) 

/usr/bin/time -o $tfile $ALVBIN alevin -lISR --chromium -1 /mnt/scratch5/avi/alevin/data/10x/v2/mohu/100/all_bcs.fq -2 /mnt/scratch5/avi/alevin/data/10x/v2/mohu/100/all_reads.fq -o $OUT/prediction -i /mnt/scratch5/avi/alevin/data/mohu/salmon_index -p 20 --tgMap /mnt/scratch5/avi/alevin/data/mohu/gtf/txp2gene.tsv  --whitelist ./alevin_test_data/alevin/whitelist.txt  && 
#--dumpUmiGraph --numCellBootstraps 100  --dumpBfh --dumpBarcodeEq  --dumpMtx --expectCells 1001  --end 6 --umiLength 10 --barcodeLength 16

cat $tfile  

tar -xvzf alevin_test_data.tar.gz  

echo "Barcodes"  
sort prediction/alevin/quants_mat_rows.txt > 1.txt  
sort alevin_test_data/alevin/quants_mat_rows.txt > 2.txt  
diff 1.txt 2.txt  > diff.txt  
wc -l diff.txt  
echo "FAILED if above line > Zero"  

echo "features"  
sort prediction/alevin/featureDump.txt > 1.txt  
sort alevin_test_data/alevin/featureDump.txt > 2.txt  
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
python alevin_test_data/src-py/alevin.py --one prediction/ --two alevin_test_data/
