import sys
import json
import matplotlib.pyplot as plt
from os import path

from pymoi.moi import *
from pymoi.parse import *
from pymoi.simdist import SimDist

def main(vcf_file, sample_id):

    fit_fig_path = f"./{sample_id}.bafscore.png"
    cluster_fig_path = f"./{sample_id}.bafplot.png"

    lines = parse_vcf(vcf_file)
    print("constructing BAF matrix")
    baf_matrix, bac_matrix = make_baf_matrix(lines, sample_id, filter_homo=True)

    gmm = finite_mixture_model(baf_matrix, max_k=5, responsibility_upper=0.1, plot_cluster_fit=True)
    plt.savefig(fit_fig_path, dpi='figure', format='png')
    plot_baf(baf_matrix, gmm)
    plt.savefig(cluster_fig_path, dpi='figure', format='png')

    resultsdict = { "sample_id" : sample_id,
                    "k" : gmm.k,
                    "fit_score" : gmm.score,
                    "boundary_prob" : gmm.prob,
                    "fit_fig_path" : fit_fig_path,
                    "cluster_fig_path" : cluster_fig_path }

    with open(f"{sample_id}.heterozygosity.json", "w") as fout:
        json.dump(resultsdict, fout, indent=4)

if __name__=="__main__":
    vcf_file = sys.argv[1]
    sample_id = sys.argv[2]
    main(vcf_file, sample_id)
