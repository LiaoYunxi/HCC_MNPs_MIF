import argparse
import sys
import scvelo as scv

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="RNAdata.loom", required=True)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="scvelo.gene.txt", required=True)
    base_group.add_argument("--paga", "--paga", type=str, dest="paga", metavar="scVelo.paga.png", required=True)
    base_group.add_argument("-p", type=int, dest="PCs", metavar="PCs", required=True)
    base_group.add_argument("-n", type=int, dest="neighbors", metavar="neighbors", required=True)
    base_group.add_argument("-t", type=int, dest="n_top_genes", metavar="n_top_genes", required=True)
    base_group.add_argument("--group", type=str, dest="group", metavar="group", required=True)
    return parser.parse_args(args)

def main(args):
    args = parse_args(args)
    f_in = args.input
    f_pc = args.PCs
    f_n = args.neighbors
    f_t = args.n_top_genes
    group_name = args.group

    uniq_prot_dict = defaultdict(list)
    adata = scv.read(f_in, cache=True)
    scv.pp.filter_and_normalize(adata, min_shared_counts=5, n_top_genes=f_t)
    scv.pp.moments(adata, n_pcs=f_pc, n_neighbors=f_n)
    scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)
    adata.uns['neighbors']['distances'] = adata.obsp['distances']
    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
    scv.tl.paga(adata, groups=group_name)   
    top_genes = adata.var['fit_likelihood'].sort_values(ascending=False)
    scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
                min_edge_width=2, node_size_scale=1.5,save=args.paga)
    top_genes[top_genes.fit_likelihood>0]. to_csv(args.output, sep="\t",index=None)

def run():
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
