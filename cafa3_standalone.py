import os
import config
from itertools import groupby
import requests
import json
import logging
import numpy as np

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

targets_dir = config.paths['targets_dir']
endpoint = config.paths['predictor_endpoint']
file_header_str = config.file_header_str

headers = {"Content-Type": "application/json"}


def iter_directory_files(maindir=targets_dir):
    for dirname, dirnames, filenames in os.walk(maindir):
        for filename in filenames:
            yield filename, os.path.join(dirname, filename)


def fasta_iter(fasta_name):
    """
    Snippet from https://github.com/brentp
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def sort_list(scores_list=None):
    scores_list = sorted(scores_list, key=lambda tup: tup[1], reverse=True)
    return scores_list

def to_list_from_dict_format(scores_dict=None, score_key='score'):
    (_k, ) = set(scores_dict[0].keys()) - set([score_key])
    return sort_list([(_[_k], _[score_key]) for _ in scores_dict])

def softmax(x):
    x = np.asarray(map(float,x))
    return np.exp(x) / np.sum(np.exp(x), axis=0)


def to_cafa_scorlist_format(scores_list=None, sample_name="T00203", n=5):
    sl = np.asarray(scores_list)
    soft_maxs = softmax(sl[:,1])
    soft_maxs /= max(soft_maxs)
    _sl = scores_list[:n]
    return "\n".join(["%s\t%s\t%.2f" % (sample_name, _n, soft_maxs[_i])
                      for _i, (_n, _s) in zip(xrange(len(_sl)),_sl) if "%.2f"%soft_maxs[_i] != "0.00"])
        
def main():
    ontologies = ['C','F']
    n = 1500/len(ontologies)
    for filename, path in (iter_directory_files()):
        taxonID = filename[len('target.'):-(len('.fasta'))]
        print 'taxonID', taxonID
        _filename = 'Corraleads_1_'+taxonID+'.txt'
        if os.path.isfile(_filename): 
            continue
        with open(_filename,'w') as fout:
            fout.write("%s" % file_header_str)
            for target_name, seq in fasta_iter(path):
                _tn = target_name.split()[0]
                logging.info("\n\nProcessing: %s" % _tn)
                data = {'seq':seq}
                for onto in ontologies:
                    _ep = endpoint+onto+'/?n='+str(n)
                    r = requests.post(_ep, data=json.dumps(data), headers=headers).json()
                    scores_dict = r["predictions"][0]["prediction"]
                    scores_list = to_list_from_dict_format(scores_dict)
                    cafa_out = to_cafa_scorlist_format(scores_list=scores_list,sample_name=_tn,n=n)
                    fout.write("%s\n" % cafa_out )
            fout.write("%s" % "END")


if __name__ == '__main__':
    main()
