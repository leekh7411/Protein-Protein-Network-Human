# Load sequences and make input file for the psi-blast program
import os
from multiprocessing import Process, Manager

def load_sequences(path):
    f = open(path, "r")
    flines = f.readlines()
    vals = []
    for i in range(0,len(flines),2):
        tag = flines[i].replace("\n","")
        tag = tag.replace(">","")
        seq = flines[i+1].replace("\n","")
        vals.append((tag, seq))
    print("- Load seq file finished (total: {})".format(len(vals)))
    return vals

# initialize pssm folder
def check_and_mkdir(path, isPrint=True):
    if not os.path.exists(path):
        os.makedirs(path)
        if isPrint: print("- Query path {} created".format(path))
    else:
        if isPrint: print("- Query path {} already exists".format(path))
            
        
def make_psiblast_querys(seqs, num_iters, out_path, base_dir, base_query):
    queries = []
    for tag, seq in seqs:
        tag = tag.replace(" ","-")
        out_pssm_path = base_dir + "/{}.pssm".format(tag)
        input_path = ".p{}.inp"
        if os.path.exists(out_pssm_path): continue
        target_query = base_query.format(input_path, num_iters, out_path, out_pssm_path)
        queries.append((target_query, tag, seq))
    print("- Initialize base dir {} queries finished (total : {})".format(base_dir, len(queries)))
    return queries

def excute_psiblast(queries, pid):
    for i, (q, tag, seq) in enumerate(queries):
        f = open(".p{}.inp".format(pid),"w")
        f.write(">{}\n{}\n".format(tag, seq))
        f.close()
        q = q.format(pid)
        pssm_output_path = q.split("out_ascii_pssm")[-1]
        pssm_output_path = pssm_output_path.replace(" ","")
        print("[P{}][{}/{}] {}".format(pid, i+1, len(queries), q))
        
        os.system(q)
        
    print("- Excution {} jobs finished".format(len(queries)))
    
def excute_psiblask_multiprocess(queries, n_jobs=4):
    batch_size = len(queries) // n_jobs
    batch_queries = []
    print("> n_jobs : {}".format(n_jobs))
    for qi in range(0, len(queries), batch_size):
        if qi+batch_size < len(queries):
            batch_queries.append(queries[qi:qi+batch_size])
            #print("batch size : {}".format(len(queries[qi:qi+batch_size])))
        else:
            batch_queries.append(queries[qi:])
            #print("batch size : {}".format(len(queries[qi:])))
    
    with Manager() as manager:
        processes = []
        for i in range(n_jobs):
            p = Process(target=excute_psiblast, args=(batch_queries[i], i))
            p.start()
            processes.append(p)
        for p in processes:
            p.join()
            
# GLOBAL VALIABLES
TRAIN_PSSM_DIR = "pssms/train"
TEST_PSSM_DIR = "pssms/test"
PSIBLAST_NUM_ITER = 0
PSIBLAST_OUT = ".null_out"
#PSIBLAST_QUERY = "psiblast -query {} -db /ncbi/db/pdbaa -num_iterations {} -out_pssm {} -out_ascii_pssm {} -out .null -comp_based_stats 1"
PSIBLAST_QUERY = "psiblast -query {} -db /ncbi/db/pdbaa -num_iterations {} -out_pssm {} -out_ascii_pssm {} -out .null"
NUM_PROCESS = 20
TRAIN_NODE_SEQUENCES = "PP-Pathways_ppi.train.txt.seq"
TEST_NODE_SEQUENCES = "PP-Pathways_ppi.test.txt.seq"

def main():
    train_seqs = load_sequences(path=TRAIN_NODE_SEQUENCES)
    test_seqs  = load_sequences(path=TEST_NODE_SEQUENCES)
    check_and_mkdir(TRAIN_PSSM_DIR)
    check_and_mkdir(TEST_PSSM_DIR)
    
    train_psiblast_queries = make_psiblast_querys(train_seqs, PSIBLAST_NUM_ITER, PSIBLAST_OUT, TRAIN_PSSM_DIR, PSIBLAST_QUERY)
    print(train_psiblast_queries[0])
    
    test_psiblast_queries = make_psiblast_querys(test_seqs, PSIBLAST_NUM_ITER, PSIBLAST_OUT, TEST_PSSM_DIR, PSIBLAST_QUERY)
    print(test_psiblast_queries[0])
    
    excute_psiblask_multiprocess(train_psiblast_queries, n_jobs=NUM_PROCESS)
    #excute_psiblask_multiprocess(test_psiblast_queries, n_jobs=NUM_PROCESS)
            
if __name__ == "__main__":
    main()