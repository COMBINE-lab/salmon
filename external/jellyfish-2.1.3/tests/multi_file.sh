#! /bin/sh

cd tests
. ./compat.sh

sort -k2,2 > ${pref}.md5sum <<EOF
d93b7678037814c256d1d9120a0e6422 ${pref}_m15_s2M.histo
d93b7678037814c256d1d9120a0e6422 ${pref}_m15_s2M_zip.histo
EOF

# Count multiple files with many readers
$JF count -t $nCPUs -F 4 -o ${pref}_m15_s2M.jf -s 2M -C -m 15 seq1m_0.fa seq1m_1.fa seq1m_2.fa seq10m.fa seq1m_3.fa seq1m_4.fa
$JF histo ${pref}_m15_s2M.jf > ${pref}_m15_s2M.histo

cat > ${pref}_gunzip_cmds <<EOF

  
  # Empty lines and comments just for fun

EOF
find . -name 'seq1m_*.fa.gz' | xargs -n 1 echo gunzip -c >> ${pref}_gunzip_cmds
$JF count -t $nCPUs -g ${pref}_gunzip_cmds -G 2 -C -m 15 -s 2M -o ${pref}_m15_s2M_zip.jf seq10m.fa
$JF histo ${pref}_m15_s2M_zip.jf > ${pref}_m15_s2M_zip.histo

# Test failure of generator
echo false > ${pref}_fail_cmds
STATUS=
$JF count -t $nCPUs -g ${pref}_fail_cmds -G 2 -C -m 15 -s 2M -o ${pref}_fail.jf || STATUS=$?
if [ -z "$STATUS" ]; then
    echo >&2 "Jellyfish did not report failure in generator command"
    false
fi

# Test failure to open file
STATUS=
$JF count -t $nCPUs -C -m 15 -s 2M -o ${pref}_fail.jf non_existent_sequence.fa || STATUS=$?
if [ -z "$STATUS" ]; then
    echo >&2 "Jellyfish did not report failure in opening a file"
    false
fi

check ${pref}.md5sum
