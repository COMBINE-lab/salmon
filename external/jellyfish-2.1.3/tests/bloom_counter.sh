#! /bin/sh

cd tests
. ./compat.sh

sort -k2,2 > ${pref}.md5sum <<EOF
9251799dd5dbd3f617124aa2ff72112a ${pref}.histo
9251799dd5dbd3f617124aa2ff72112a ${pref}_filtered.histo
EOF

cat > ${pref}_commands <<EOF
gunzip -c seq1m_0.fa.gz
gunzip -c seq1m_0.fa.gz
EOF
$JF bc -t $nCPUs -o ${pref}.bc -s 1M -C -m 40 --timing ${pref}.timing -g ${pref}_commands -G 2
# Counting without filtering
$JF count -t $nCPUs -o ${pref}.jf -s 2M -C -m 40 seq1m_0.fa
# Filtereing should not do anything here: all mers are loaded twice in the bc
$JF count -t $nCPUs -o ${pref}_filtered.jf --bc ${pref}.bc -s 2M -C -m 40 seq1m_0.fa

$JF bc -t $nCPUs -o ${pref}_none.bc -s 2M -C -m 40 seq1m_0.fa seq1m_1.fa seq1m_1.fa
$JF count -t $nCPUs -o ${pref}_none.jf --bc ${pref}_none.bc -s 1M -C -m 40 seq1m_0.fa
$JF histo ${pref}_none.jf > ${pref}_none.histo

$JF histo ${pref}.jf > ${pref}.histo
$JF histo ${pref}_filtered.jf > ${pref}_filtered.histo

TOTAL=$(cut -d\  -f2 ${pref}.histo)
COLLISION=$(cut -d\  -f2 ${pref}_none.histo)
# FPR is 1 in 1000. Should not get more than 1/500 collisions.
[ $((TOTAL / 500 > COLLISION)) = 1 ] || {
    echo >&2 "Too many collisions"
    false
}

QUERY_TOT=$($JF query -s seq1m_0.fa ${pref}.bc | grep -c ' 2$')
[ $QUERY_TOT = $TOTAL ] || {
    echo >&2 "Queried count 2 mers should be all mers"
    false
}
QUERY_COL=$($JF query -s seq1m_0.fa ${pref}_none.bc | grep -c ' 2$')
[ $QUERY_COL = $COLLISION ] || {
    echo >&2 "Queried count 2 mers should equal collisions"
    false
}

check ${pref}.md5sum
