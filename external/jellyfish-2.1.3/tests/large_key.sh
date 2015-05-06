#! /bin/sh

cd tests
. ./compat.sh

sort -k2,2 > ${pref}.md5sum <<EOF
ded3925fe6bbaca10accc10d1bde11b5 ${pref}_m100_2M_ordered
ded3925fe6bbaca10accc10d1bde11b5 ${pref}_m100_2k_ordered
ded3925fe6bbaca10accc10d1bde11b5 ${pref}_m100_2k_disk_ordered
EOF

head -n 10001 seq1m_0.fa | time $JF count -t $nCPUs -o ${pref}_m100_2M.jf -s 2M -m 100 /dev/fd/0
$JF dump -c ${pref}_m100_2M.jf | cut -d\  -f 1 | sort > ${pref}_m100_2M_ordered
head -n 10001 seq1m_0.fa | time $JF count -t $nCPUs -o ${pref}_m100_2k.jf -s 2k -m 100 /dev/fd/0
$JF dump -c ${pref}_m100_2k.jf | cut -d\  -f 1 | sort > ${pref}_m100_2k_ordered
head -n 10001 seq1m_0.fa | time $JF count -t $nCPUs -o ${pref}_m100_2k_disk.jf -s 2k --disk -m 100 /dev/fd/0
$JF dump -c ${pref}_m100_2k_disk.jf | cut -d\  -f 1 | sort > ${pref}_m100_2k_disk_ordered
check ${pref}.md5sum
