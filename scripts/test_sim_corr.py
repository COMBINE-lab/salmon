import argparse
import errno    
import os


def mkdir_p(path):
    """
    http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def computeSimMetrics(args):
    import pandas as pd
    import json
    simFile = args.sim
    estFile = args.est
    sd = pd.read_table(simFile).set_index('Name')
    ed = pd.read_table(estFile).set_index('Name')
    d = sd.join(ed, rsuffix='_est')

    res = {}
    res['pearson'] = d['TPM'].corr(d['TPM_est'])
    res['spearman'] = d['TPM'].corr(d['TPM_est'], method='spearman')

    import os.path
    pdir = os.path.dirname(args.out)
    if not os.path.exists(pdir):
        mkdir_p(pdir)

    with open(args.out, 'w') as outfile:
        json.dump(res, outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test simulated accuracy')
    parser.add_argument('--sim', type=str, help='path to simulated tpm')
    parser.add_argument('--est', type=str, help='path to estimated tpm')
    parser.add_argument('--out', type=str, help='where to write the output')
    args = parser.parse_args()
    computeSimMetrics(args)


