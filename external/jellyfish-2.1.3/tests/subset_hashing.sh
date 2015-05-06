#! /bin/sh

cd tests/
. ./compat.sh

sort -k2,2 > ${pref}.md5sum <<EOF
bd7a5f6ba000b282cd79cb9f342e7ede ${pref}_m35_s2M_if.histo
8eb6d4a50aeba178e4847c2da71dbb70 ${pref}_m10_s2M_if.histo
EOF

# Partial count (in _0 and _2) with 35-mers
$JF count -t $nCPUs -o ${pref}_m35_s2M_if.jf -s 2M -C -m 35 --if seq1m_0.fa --if seq1m_2.fa seq1m_1.fa seq1m_0.fa seq1m_3.fa seq1m_2.fa
$JF histo ${pref}_m35_s2M_if.jf > ${pref}_m35_s2M_if.histo

# Idem with 10-mers
$JF count -t $nCPUs -o ${pref}_m10_s2M_if.jf -s 6M -C -m 10 --if seq1m_0.fa --if seq1m_2.fa seq1m_1.fa seq1m_0.fa seq1m_3.fa seq1m_2.fa
$JF histo ${pref}_m10_s2M_if.jf > ${pref}_m10_s2M_if.histo

check ${pref}.md5sum
