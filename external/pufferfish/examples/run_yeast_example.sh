### Run it from the examples directory. Addresses are relative

mkdir -p yeastchr01
wget -O yeastchr01/yeastchr01.fasta https://downloads.yeastgenome.org/sequence/S288C_reference/chromosomes/fasta/chr01.fsa
../build/src/pufferfish index -r yeastchr01/yeastchr01.fasta -o yeastchr01/yeast_pi -s -e 8
../build/src/pufferfish validate -i yeastchr01/yeast_pi
