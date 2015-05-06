"""Compute Mutation Rate

Usage:
  ComputeMutationRate.py <input1> <input2>


Options:
  -h --help               Show this message.
"""
from docopt import docopt
from Bio import SeqIO
import itertools


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def main(args):
    in1 = args['<input1>']
    in2 = args['<input2>']

    seq1 = SeqIO.parse(in1, 'fasta')
    seq2 = SeqIO.parse(in2, 'fasta')
    nMut = 0
    totLen = 0
    i = 0

    for s1, s2 in itertools.izip(seq1, seq2):
        o1 = ""
        o2 = ""
        newMut = True
        assert(s1.name == s2.name)
        totLen += len(s1.seq)
        for a,b in itertools.izip(s1.seq, s2.seq):
            if a != b:
                nMut+= 1
                newMut = True
                o1 += bcolors.OKGREEN + a + bcolors.ENDC
                o2 += bcolors.FAIL + b + bcolors.ENDC
            else:
                if newMut:
                    o1 += "***"
                    o2 += "***"
                    newMut = False

        i += 1
        print(o1)
        print(o2)
        if i % 1000 == 0:
            print("There were {} mutations in {} bases; total rate = {:0.2f}\r\r".format(\
                   nMut, totLen, (100.0 * nMut) / totLen))
    print("There were {} mutations in {} bases; total rate = {:0.2f}\r\r".format(\
           nMut, totLen, (100.0 * nMut) / totLen))

if __name__ == "__main__":
    args = docopt(__doc__, version="Compute Mutation Rate v1.0")
    main(args)
