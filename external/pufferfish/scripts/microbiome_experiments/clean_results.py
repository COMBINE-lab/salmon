import sys
import pandas as pd
import numpy as np

def extract_perf(infile, outfile, memPow=20):
    if (memPow == 20):
        memColName = 'memory(G)'
    elif (memPow == 10):
        memColName = 'memory(M)'
    else:
        print('Memmory conversion is neither to Gig nor to Meg')
        exit()
    ssp = pd.read_csv('{}'.format(infile), sep=' ', header=None,
                               names=['cold-warm', 'option', 'sample', 'tool', memColName, 'total_time', 'alignment_time'])
    ssp[memColName] = np.round(ssp[memColName]/2**memPow, 2)
    sspg = ssp.groupby(['cold-warm', 'sample', 'option', 'tool']).first()
    sspg.to_excel('{}'.format(outfile))

def extract_quant(infile, outfile):
	with open('{}'.format(infile)) as f:
	    alignment_info = []
	    alignments = []
	    for l in f:
	        sp=l.split()
	        if ((sp[0].strip() == 'covid19') or (sp[0].strip() == 'sars') or (sp[0].strip() == 'bat2008')):
	            data = sp[0].strip()
	            tool = sp[1].strip()
	            option = sp[2].strip()
	        elif (len(sp) == 1):
	            aligned = int(sp[0].strip())
	        else:
	            alignments += [[data, tool, option, aligned, int(sp[0].strip())/2, sp[1].strip()]]
	ss = pd.DataFrame.from_records(alignments)
	ss.columns = ['data', 'tool', 'option', 'aligned_reads', 'read_cnt', 'ref']
	sameSpecies = pd.pivot_table(ss, index=['data', 'option', 'tool', 'aligned_reads'], 
                             columns='ref', values='read_cnt', aggfunc=np.sum).fillna(0)
	sameSpecies.to_excel('{}'.format(outfile))


if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 4:
        print("""\
        
        This script is part of the metagenomic bulk experiment for Puffaligner and follows the required directory structure for it.
        Converts performance results into excel format with groups and merged-cells

            Required packages:
                numpy
                pandas
                openpyxl

        Usage:  python <script> <perf/quant> <input_file> <output_file(must have xlsx extension)>
        """)
        sys.exit(0)
    isPerf = False
    if (sys.argv[1] == "perf"):
        isPerf=True
    elif (sys.argv[1] != "quant"):
        print("ERROR! First argument can only have values \"perf\" or \"quant\".")
        sys.exit(0)
    
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    if (isPerf):
        extract_perf(input_file, output_file, 20)
    else:
        extract_quant(input_file, output_file)

