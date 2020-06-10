import os
import numpy as np
import pandas as pd
from collections import defaultdict

def load_interaction_sequences(path):
    f = open(path, "r")
    flines = f.readlines()
    pairs = []
    for line in flines:
        line = line.replace("\n","")
        vals = line.split("\t")
        seq1, seq2 = vals[0], vals[1]
        pairs.append((seq1, seq2))
    print("- Load seq pair file finished (size: {})".format(len(pairs)))
    return pairs


def load_sequence_id_dict(path, pssm_dir="pssms/train/{}.pssm"):
    f = open(path, "r")
    flines = f.readlines()
    seq_dict = {}
    for i in range(0,len(flines),2):
        tag = flines[i].replace("\n","")
        tag = tag.replace(">","")
        tag = tag.replace(" ","-")
        seq = flines[i+1].replace("\n","")
        seq_dict[seq] = (tag, pssm_dir.format(tag))
        #print("[{}]".format(i) + "-"*100)
        #print(tag, pssm_dir.format(tag))
        #print(seq)
        
    print("- Load seq_dict finished (total: {})".format(len(seq_dict)))
    return seq_dict

def min_max_norm(x):
    min_x = np.min(x)
    max_x = np.max(x)
    if max_x == 0:
        return x
    return (x-min_x)/(max_x-min_x)

def aa_onehot_dict():
    aa = "ARNDCQEGHILKMFPSTWYVU" # 21 amino-acids (+ selenocysteine U) but the PSSM only contains 20 amino-acids without 'U'
    def onehot(index):
        o = np.zeros(len(aa))
        o[index] = 1
        return o
    aa_onehot = {a: onehot(i) for i,a in enumerate(aa)}
    return aa_onehot

def pssm_feature(target_seq, pssm_path):
    try:
        f = open(pssm_path, "r")
    except:
        print("> {} not exists(not hit in PSIBLAST in PDB), so we do not use this data".format(pssm_path))
        return []
    flines = f.readlines()
    lines = []
    aa_onehot = aa_onehot_dict()
    
    LINE_HEAD  = 3
    LINE_START = LINE_HEAD
    LINE_END   = LINE_HEAD + len(target_seq)
    
    for i, line in enumerate(flines):
        line = line.replace("\n","")
        lines.append(line)
        #print("[line {}]".format(i), line)
    pssm_lines = lines[LINE_START:LINE_END]
    
    COL_INDEX          = 0
    COL_CHAR           = 1
    COL_RAW_START      = 2
    COL_RAW_END        = 2+19
    COL_WEIGHTED_START = 2+19+1
    COL_WEIGHTED_END   = 2+19+1+19
    COL_INFO_PER_POS   = 2+19+1+19+1
    COL_PSEUDO         = 2+19+1+19+1+1
    
    pssm_r = []
    pssm_w = []
    seq_onehot = []
    info_per_poses = []
    for line in pssm_lines:
        vals         = line.split()
        seq_index    = vals[COL_INDEX]
        seq_char     = vals[COL_CHAR]
        seq_raw_pssm = np.array(vals[COL_RAW_START:COL_RAW_END]).astype(float)
        seq_raw_pssm = min_max_norm(seq_raw_pssm)
        seq_weighted = np.array(vals[COL_WEIGHTED_START:COL_WEIGHTED_END]).astype(float) / 100
        seq_info_per_pos = float(vals[COL_INFO_PER_POS])
        seq_pseudo_count = float(vals[COL_PSEUDO])
        seq_onehot.append(aa_onehot[seq_char])
        pssm_r.append(seq_raw_pssm)
        pssm_w.append(seq_weighted)
        info_per_poses.append([seq_info_per_pos])
        
        """
        print(line)
        print(seq_index, seq_char)
        print(seq_raw_pssm)
        print(seq_weighted)
        print(seq_info_per_pos)
        print(seq_pseudo_count)
        """
    
    pssm_r = np.array(pssm_r)
    pssm_w = np.array(pssm_w)
    seq_onehot = np.array(seq_onehot)
    info_per_poses = np.array(info_per_poses)
    pssm_x = np.concatenate([seq_onehot, pssm_r, info_per_poses], axis=-1)
    #print(pssm_x)
    
    return pssm_x



def init_pssm_dict(seq_dict):
    pssm_dict = defaultdict(lambda: {"available": False, "pssm": None})
    num_items = len(seq_dict.items())
    for i, (seq, (seq_id, seq_pssm_path)) in enumerate(seq_dict.items()):
        pssm = pssm_feature(seq, seq_pssm_path)
        if len(pssm) > 0:
            pssm_dict[seq]["available"] = True
            pssm_dict[seq]["pssm"] = pssm
        if i % 1000 == 0:
            print("> pssm conversion processing .... {} / {}".format(i, num_items))
    print("- initialize pssm dict from seq dict (total: {}, remains: {})".format(len(seq_dict), len(pssm_dict)))
    return dict(pssm_dict)


# GLOBAL VALIABLES
TRAIN_NPZ_PATH = "npz/pp-pathway-train-lang.small.k0.pssm.npz"
TEST_NPZ_PATH  = "npz/pp-pathway-test-lang.small.k0.pssm.npz"
TRAIN_PAIRS = "PP-Pathways_ppi.train.txt"
TEST_PAIRS = "PP-Pathways_ppi.test.txt"
TRAIN_NODE_SEQUENCES = "PP-Pathways_ppi.train.txt.seq"
TEST_NODE_SEQUENCES = "PP-Pathways_ppi.test.txt.seq"
IS_TRAINSET = True

if __name__ == "__main__":
    train_pairs = load_interaction_sequences(path=TRAIN_PAIRS)
    test_pairs  = load_interaction_sequences(path=TEST_PAIRS)
    train_seq_dict = load_sequence_id_dict(path=TRAIN_NODE_SEQUENCES, pssm_dir="pssms/train/{}.pssm")
    test_seq_dict = load_sequence_id_dict(path=TEST_NODE_SEQUENCES, pssm_dir="pssms/test/{}.pssm")
    
    if IS_TRAINSET:
        train_pssm_dict = init_pssm_dict(train_seq_dict)
        np.savez(TRAIN_NPZ_PATH, pairs=np.array(train_pairs), seq_dict=train_seq_dict, pssm_dict=train_pssm_dict)
    else:
        test_pssm_dict = init_pssm_dict(test_seq_dict)
        np.savez(TEST_NPZ_PATH, pairs=np.array(test_pairs), seq_dict=test_seq_dict, pssm_dict=test_pssm_dict)

    # EXAMPLE - Load dataset
    """
    d = np.load(TRAIN_NPZ_PATH)
    _train_pairs = d["pairs"]
    _train_seq_dict = d["seq_dict"][()]
    _train_pssm_dict = d["pssm_dict"][()]
    print(_train_pairs.shape)
    print(len(_train_seq_dict))
    print(len(_train_pssm_dict))
    """