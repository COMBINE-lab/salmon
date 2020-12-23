ALVBIN=$1 #"/mnt/scratch6/salmon_ci/COMBINE-lab/salmon/master/build/salmon-latest_linux_x86_64/bin/salmon"
BASEDIR="/mnt/scratch6/avi/alevin/alevin"
OUT=$PWD

tfile=$(mktemp /tmp/foo.XXXXXXXXX) 

/usr/bin/time -o $tfile $ALVBIN alevin -lISR --chromium -1 ${BASEDIR}/data/10x/v2/mohu/100/all_bcs.fq.gz -2 ${BASEDIR}/data/10x/v2/mohu/100/all_reads.fq.gz -o $OUT/prediction -i ${BASEDIR}/data/mohu/salmon_index -p 20 --tgMap ${BASEDIR}/data/mohu/gtf/txp2gene.tsv --dumpMtx --no-version-check --dumpFeatures --dumpArborescence --writeMappings=$OUT/prediction/with_bug.sam #--whitelist ./alevin_test_data/alevin/quants_mat_rows.txt
#--dumpBfh --whitelist /mnt/scratch5/avi/alevin/bin/salmon/tests/whitelist.txt
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
