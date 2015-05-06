#! /bin/sh

cd tests
. ./compat.sh

sort -k2,2 > ${pref}.md5sum <<EOF
72f1913b3503114c7df7a4dcc68ce867 ${pref}_m40_s16m.histo
72f1913b3503114c7df7a4dcc68ce867 ${pref}_automerge_m40_s1m.histo
72f1913b3503114c7df7a4dcc68ce867 ${pref}_m40_s1m_merged.histo
72f1913b3503114c7df7a4dcc68ce867 ${pref}_m40_s1m_text.histo
EOF

FILES="seq1m_0.fa seq1m_1.fa seq1m_0.fa seq1m_2.fa seq1m_2.fa"
echo $FILES | xargs $JF count -t $nCPUs -o ${pref}_m40_s16m.jf -s 4M -C -m 40
$JF histo ${pref}_m40_s16m.jf > ${pref}_m40_s16m.histo

ls | grep "^${pref}_m40_s1m[0-9].*" | xargs rm -f
echo $FILES | xargs $JF count -t $nCPUs -o ${pref}_m40_s1m -s 1M --disk --no-merge -C -m 40
$JF merge -o ${pref}_m40_s1m_merged.jf ${pref}_m40_s1m[0-9]*
ls | grep "^${pref}_m40_s1m[0-9].*" | xargs rm -f

echo $FILES | xargs $JF count -t $nCPUs -o ${pref}_automerge_m40_s1m.jf -s 1M --disk -C -m 40

echo $FILES | xargs $JF count -t $nCPUs -o ${pref}_m40_s1m_text.jf -s 1M --text --disk -C -m 40

$JF histo ${pref}_automerge_m40_s1m.jf > ${pref}_automerge_m40_s1m.histo
$JF histo ${pref}_m40_s1m_merged.jf > ${pref}_m40_s1m_merged.histo
$JF histo ${pref}_m40_s1m_text.jf > ${pref}_m40_s1m_text.histo

check ${pref}.md5sum
