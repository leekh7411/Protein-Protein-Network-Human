from collections import defaultdict
import pandas as pd
import random

def uniprot_parser(key):
    kvals = key.split("|")
    gene_id = kvals[1]
    return gene_id

def fasta_parser(path, outpath):
    f = open(path, "r")
    flines = f.readlines()
    ret = {}
    keys, vals = [],[]
    for i, line in enumerate(flines):
        line = line.replace("\n","")
        if line[0] == ">":
            key = line[1:]
            keys.append(uniprot_parser(key))
            if i == 0 :
                val = ""
            else:
                vals.append(val)
                val = ""
        else:
            val += line
    vals.append(val)
    #print("num keys : {} / num vals : {}".format(len(keys), len(vals)))
    ret = {k:v for k,v in zip(keys, vals)}
    f.close()
    
    o = open(outpath,"w")
    for k,v in zip(keys, vals):
        o.write("{}\t{}\n".format(k,v))
    o.close()
    return ret

def load_entrez2uniprot_dict(path):
    f = open(path,"r")
    flines = f.readlines()
    sdict = {}
    for line in flines:
        line = line.replace("\n","")
        vals = line.split("\t")
        entrez_id, uniprot_id = vals[0], vals[1]
        sdict[entrez_id] = uniprot_id
    
    print("- Load Entrez-ID to Uniprot-ID converter finish (size : {})".format(len(sdict)))
    return sdict

def load_interaction_pairs(path):
    f = open(path, "r")
    flines = f.readlines()
    pairs = []
    ids = []
    for line in flines:
        line = line.replace("\n","")
        vals = line.split(",")
        ent1, ent2 = vals[0], vals[1]
        pairs.append((ent1, ent2))
        ids.append(ent1)
        ids.append(ent2)
    ids = list(set(ids))
    print("- Load {} pairs complete (total {} genes)".format(len(pairs), len(ids)))
    return pairs, ids

def convert_interactions_entrez_to_uniprot(interactios, E2U):
    uniprot_interactions = []
    for ent1, ent2 in interactions:
        uni1 = E2U[ent1]
        uni2 = E2U[ent2]
        uniprot_interactions.append((uni1, uni2))
    print("- Converting interaction IDs finished (size : {})".format(len(uniprot_interactions)))
    return uniprot_interactions

def link_interactions_and_uniprot_sequences(interactions, entz_interactions, uniprot_seqdict):
    num_failed = 0
    interactions_with_seq = []
    for (uni1, uni2), (ent1, ent2) in zip(interactions, entz_interactions):
        try:
            seq1 = uniprot_seqdict[uni1]
            seq2 = uniprot_seqdict[uni2]
            interactions_with_seq.append(((ent1, uni1, seq1), (ent2, uni2, seq2)))
        except:
            num_failed += 1
    print("- Linking interaction IDs to Uniprot sequences finished (num failed : {} and success : {})".format(num_failed, len(interactions_with_seq)))
    return interactions_with_seq

def save_interactions_as_csv(interactions_with_seq, path):
    df_dict = defaultdict(lambda: [])
    for (e1, u1, seq1), (e2, u2, seq2) in interactions_with_seq:
        df_dict["Entrez-ID 1"].append(e1)
        df_dict["Entrez-ID 2"].append(e2)
        df_dict["Uniprot-ID 1"].append(u1)
        df_dict["Uniprot-ID 2"].append(u2)
        df_dict["Uniprot Protein 1"].append(seq1)
        df_dict["Uniprot Protein 2"].append(seq2)
    df = pd.DataFrame.from_dict(df_dict)
    df.to_csv(path, index=False)
    print("- Preprocessed interaction finished (size : {})".format(len(df_dict)))
    return df

def k_mer_split(seq, k):
    ret = ""
    for i in range(0, len(seq)-k+1, k):
        ret += seq[i:i+k]
        if i == (len(seq)-k):
            ret += "."
        else:
            ret += " "
    return ret

def save_interactions_as_lang(interactions_with_seq, train_path, test_path, train_test_split=0.8, total_num=-1, len_constraint=-1):
    random.seed(0)
    random.shuffle(interactions_with_seq)
    if total_num > 0:
        interactions_with_seq = interactions_with_seq[:total_num]
    else:
        interactions_with_seq = interactions_with_seq
    
    if len_constraint > 0:
        _interactions_with_seq = []
        for (e1, u1, seq1), (e2, u2, seq2) in interactions_with_seq:
            if len(seq1) < len_constraint and len(seq2) < len_constraint:
                _interactions_with_seq.append(((e1, u1, seq1), (e2, u2, seq2)))
            interactions_with_seq = _interactions_with_seq
        print("- Total {} remains (max length : {})".format(len(interactions_with_seq), len_constraint))
            
    num_train = int(len(interactions_with_seq) * train_test_split)
    
    train_inters = interactions_with_seq[:num_train]
    test_inters  = interactions_with_seq[num_train:]
    train_seqs = {}
    test_seqs = {}
    
    f = open(train_path, "w")
    
    
    for (e1, u1, seq1), (e2, u2, seq2) in train_inters:
        train_seqs[seq1] = 1
        train_seqs[seq2] = 1
        #kseq1 = k_mer_split(seq1, k=3)
        #kseq2 = k_mer_split(seq2, k=3)
        kseq1 = seq1
        kseq2 = seq2
        f.write("{}\t{}\n".format(kseq1, kseq2))
    f.close()
    print("- Lang version trainset saved (size : {})".format(len(train_inters)))
    
    fseq = open(train_path + ".seq", "w")
    train_seqs = list(train_seqs.keys())
    for i, seq in enumerate(train_seqs):
        fseq.write(">train {}\n".format(i+1))
        fseq.write("{}\n".format(seq))
    fseq.close()
    print("- Trainset sequences saved (size : {})".format(len(train_seqs)))
    
    f = open(test_path, "w")
    for (e1, u1, seq1), (e2, u2, seq2) in test_inters:
        test_seqs[seq1] = 1
        test_seqs[seq2] = 1
        #kseq1 = k_mer_split(seq1, k=3)
        #kseq2 = k_mer_split(seq2, k=3)
        kseq1 = seq1
        kseq2 = seq2
        f.write("{}\t{}\n".format(kseq1, kseq2))
    f.close()
    print("- Lang version testset saved (size : {})".format(len(test_inters)))
    
    
    fseq = open(test_path + ".seq", "w")
    test_seqs = list(test_seqs.keys())
    for i, seq in enumerate(test_seqs):
        fseq.write(">test {}\n".format(i+1))
        fseq.write("{}\n".format(seq))
    fseq.close()
    print("- Testset sequences saved (size : {})".format(len(test_seqs)))
    

if __name__ == "__main__":
    
    uniprot_seq_dict = fasta_parser("uniprot_sprot.fasta", "uniprot_spot.pair.txt")
    print("- Load uniprot protein sequence database finished (number of seqs : {})".format(len(uniprot_seq_dict)))
    e2u_idict = load_entrez2uniprot_dict("PP-Pathways_ppi.entrez2uniprot.txt")
    interactions, gene_ids = load_interaction_pairs("PP-Pathways_ppi.csv")
    uniprot_interactions = convert_interactions_entrez_to_uniprot(interactions, e2u_idict)
    inters_with_seq = link_interactions_and_uniprot_sequences(uniprot_interactions, interactions, uniprot_seq_dict)
    interaction_df = save_interactions_as_csv(inters_with_seq, "PP-Pathways_ppi.Uniprot.seq.csv")
    save_interactions_as_lang(inters_with_seq, "PP-Pathways_ppi.train.txt", "PP-Pathways_ppi.test.txt", 0.8, 30000, 1000)
