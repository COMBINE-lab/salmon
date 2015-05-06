#! /bin/sh

cd tests
. ./compat.sh

sort -k2,2 > ${pref}.md5sum <<EOF 
864c0b0826854bdc72a85d170549b64b ${pref}_m15_s2M.histo
864c0b0826854bdc72a85d170549b64b ${pref}_m15_s16M.histo
41fd8408dde0ea14bec7425b1a877140 ${pref}_m15.stats
376761a6e273b57b3428c14e3b536edf ${pref}_binary.dump
376761a6e273b57b3428c14e3b536edf ${pref}_text.dump
9251799dd5dbd3f617124aa2ff72112a ${pref}_binary.histo
c30cba4fe2886cea4abb27f5c30ea35e ${pref}_binary.stats
9251799dd5dbd3f617124aa2ff72112a ${pref}_text.histo
c30cba4fe2886cea4abb27f5c30ea35e ${pref}_text.stats
94625cd2d59e278f08421a673eb0926a ${pref}_m15_s2M_L2_U3.histo
94625cd2d59e278f08421a673eb0926a ${pref}_m15_s2M_L2_U3_automerge.histo
45fb383344e0fb0b7540718339be4c03 ${pref}_query_one_count
EOF

# Count with in memory hash doubling
$JF count -t $nCPUs -o ${pref}_m15_s2M.jf -s 2M -C -m 15 seq10m.fa
$JF histo ${pref}_m15_s2M.jf > ${pref}_m15_s2M.histo
$JF stats ${pref}_m15_s2M.jf > ${pref}_m15.stats

# Count without size doubling
$JF count -t $nCPUs -o ${pref}_m15_s16M.jf -s 16M -C -m 15 seq10m.fa
$JF histo ${pref}_m15_s16M.jf > ${pref}_m15_s16M.histo

# Count large merges in binary and text. Should agree
$JF count -m 40 -t $nCPUs -o ${pref}_text.jf -s 2M --text seq1m_0.fa
$JF count -m 40 -t $nCPUs -o ${pref}_binary.jf -s 2M seq1m_0.fa
$JF histo ${pref}_text.jf > ${pref}_text.histo
$JF stats ${pref}_text.jf > ${pref}_text.stats
$JF histo ${pref}_binary.jf > ${pref}_binary.histo
$JF stats ${pref}_binary.jf > ${pref}_binary.stats
$JF dump -c ${pref}_text.jf | sort > ${pref}_text.dump
$JF dump -c ${pref}_binary.jf | sort > ${pref}_binary.dump

# Check the lower and upper count without merging
$JF count -t $nCPUs -o ${pref}_m15_s2M_L2_U3.jf -s 2M -C -m 15 -L2 -U3 seq10m.fa
$JF histo ${pref}_m15_s2M_L2_U3.jf > ${pref}_m15_s2M_L2_U3.histo

# Check the lower and upper count limits with merging
$JF count -t $nCPUs -o ${pref}_m15_s2M_L2_U3_automerge.jf -s 2M -C -m 15 -L2 -U3 --disk seq10m.fa
$JF histo ${pref}_m15_s2M_L2_U3_automerge.jf > ${pref}_m15_s2M_L2_U3_automerge.histo

# Check query
$JF query ${pref}_binary.jf -s seq1m_0.fa    | grep ' 1$' | wc -l | sed -e 's/ //g' > ${pref}_query_one_count
# $JF query ${pref}_binary.jf -s seq1m_0.fa -C | grep ' 1$' | wc -l | sed -e 's/ //g' > ${pref}_query_canonical_one_count

# $JF count -m 40 -t $nCPUs -o ${pref}_text -s 2M --text seq1m_0.fa
# $JF info -s ${pref}_text0 | sort >  ${pref}_text0_dump

# $JF count -m 40 -t $nCPUs -o ${pref}_bin -s 2M seq1m_0.fa
# $JF dump -c ${pref}_bin0 | sort > ${pref}_bin0_dump

# $JF histo ${pref}_bin0 > ${pref}.histo

# Count all k-mers
# $JF count --matrix seq10m_matrix_22 -m 22 -t $nCPUs -o $pref \
#     -s 10000000 --timing ${pref}.timing --stats ${pref}.stats seq10m.fa
# # Output only counts >= 2
# $JF count --matrix seq10m_matrix_22 -m 22 -t $nCPUs -o ${pref}_L \
#     -s 10000000 --timing ${pref}.timing -L 2 seq10m.fa
# # Stream output
# $JF count --matrix seq10m_matrix_22 -m 22 -t $nCPUs -o /dev/fd/1 -O \
#     -s 10000000 --timing ${pref}.timing --stats ${pref}.stats seq10m.fa | \
#     cat > ${pref}_S

# $JF histo -f ${pref}_0 > ${pref}.histo
# $JF dump -c ${pref}_L_0 > ${pref}_L.dump
# $JF dump -c -L 2 ${pref}_0 > ${pref}.dump
# cat ${pref}_0 | $JF dump -c -L 2 /dev/fd/0 > ${pref}_stream.dump
# echo "GCCATTTCGATTAAAGAATGAT TAGGCATGCAACGCTTCCCTTT" | $JF query ${pref}_0 > ${pref}.query
# $JF dump -c ${pref}_0 > ${pref}_all.dump
# $JF dump -c ${pref}_S > ${pref}_S_all.dump

check ${pref}.md5sum

# cat ${pref}.timing

