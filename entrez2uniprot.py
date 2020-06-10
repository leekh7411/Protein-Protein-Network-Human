from PyEntrezId import Conversion

# include your email address
Id = Conversion('__@__.__.__')

def read_pp_pathway_entrez_ids():
    f = open("PP-Pathways_ppi.geneids.txt", "r")
    flines = f.readlines()
    gene_ids = [gid.replace("\n","") for gid in flines]
    return gene_ids
gene_ids = read_pp_pathway_entrez_ids() 

def save_entrez_ids_to_uniprot(gene_ids, outpath):
    f = open(outpath, "w")
    for i, gid in enumerate(gene_ids):
        #print(gid)
        try:
            uid = Id.convert_entrez_to_uniprot(gid)
        except:
            uid = "NULL"
            
        f.write("{}\t{}\n".format(gid, uid))
        print("[{}/{}] {} -> {}".format(i+1, len(gene_ids), gid, uid))
    f.close()

if __name__ == "__main__":
    save_entrez_ids_to_uniprot(gene_ids, "PP-Pathways_ppi.entrez2uniprot.txt")