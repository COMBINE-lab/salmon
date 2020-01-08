"""Mutate Reference.

Usage:
  MutateReference.py --in=<input> --out=<output> --rate=<rate>

Options:
  -h --help               Show this message.
  --in=<input>            The input reference sequence.
  --out=<output>          The output [mutated] sequences.
  --rate=<rate>           The rate of mutation.
"""
from docopt import docopt
from Bio import SeqIO
import random
import copy
import math
import sys

mutList = {'A' : ['T', 'C', 'G'],
           'T' : ['A', 'C', 'G'],
           'C' : ['A', 'G', 'T'],
           'G' : ['A', 'C', 'T'],
           'a' : ['c', 'g', 't'],
           't' : ['a', 'c', 'g'],
           'c' : ['a', 'g', 't'],
           'g' : ['a', 'c', 't'],
           'N' : ['A', 'T', 'C', 'G'],
           'n' : ['a', 't', 'c', 'g']}

def mutate(base):
    ml = mutList[base]
    return ml[random.randint(0, len(ml)-1)]

def randCeilFloor(x):
    if x < 1.0 or random.random() < 0.5:
        return int(math.ceil(x))
    else:
        return int(math.floor(x))

def main(args):
    mutRate = float(args['--rate'])
    refIn = args['--in']
    refOut = args['--out']

    ofile = open(refOut,'a')
    seqs = SeqIO.parse(refIn, 'fasta')
    currBase = 0
    nMut = 0
    totLen = 0
    nextMut = currBase + randCeilFloor(random.expovariate(mutRate))
    for i, s in enumerate(seqs):
        lenS = len(s.seq)
        firstBase = currBase
        lastBase = currBase + lenS
        totLen += lenS
        o = copy.copy(s.seq)
        x = s.seq.tomutable()
        while nextMut < lastBase:
            offset = nextMut - firstBase
            orig = o[offset]
            x[offset] = mutate(orig)
            assert(x[offset] != orig)
            nMut += 1
            currBase = nextMut
            nextMut = currBase + randCeilFloor(random.expovariate(mutRate))
        currBase = lastBase
        #assert(len(o) == len(x))
        #for a,b in itertools.izip(o,x):
        #    if a != b: nMut += 1
        s.seq = x
        SeqIO.write(s, ofile, 'fasta')
        if (i % 1000 == 0):
            sys.stderr.write("processed {} records; performed {} mutations; rate = {:.2f}%\r\r".format(i, nMut, (100.0 * nMut) / totLen))
        # print("\n")/tmp/hg18_transcripts_mut0.001.fa/tmp/hg18_transcripts_mut0.001.fa/tmp/hg18_transcripts_mut0.001.fa
        # print('Num Rec = {}, Num Mut = {}, Tot. length = {}, rate = {}'.format(i, nMut, totLen, (100.0 * nMut) / totLen))
    print("\ndone.")

if __name__ == "__main__":
    arguments = docopt(__doc__, version="Mutate Reference 1.0")
    main(arguments)
