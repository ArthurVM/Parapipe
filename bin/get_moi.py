import sys
import json
import matplotlib.pyplot as plt
import pandas as pd
from os import path

from pymoi.VCF import *
from pymoi.plot import *


def gen_empty_plot(save_path):
    # Create an empty plot
    fig, ax = plt.subplots()

    # Add the text to the plot
    text = "too few heterogeneous alleles detected"
    ax.text(0.5, 0.5, text, ha='center', va='center', fontsize=16, color='red')

    # Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Save the plot as a PNG image
    plt.savefig(save_path, bbox_inches='tight', pad_inches=0)


def main(vcf_file, fws_csv):

    fws_df = pd.read_csv(fws_csv, index_col=0)

    vcf_obj = VCF(vcf_file)

    for sample_id in vcf_obj.getSampleIds():

        allele_df = vcf_obj.getAlleleCounts(sample_id)

        ## get only heterogeneous alleles
        het_allele_df = allele_df[allele_df["alt_freq"]<0.95]

        ## check that there are a sufficient number of heterogeneous alleles
        if len(het_allele_df) >= 20:
            gmm = finite_mixture_model(het_allele_df, max_k=5, responsibility_upper=0.1, plot_cluster_fit=True)
            plt.savefig(f"./{sample_id}.bafscore.png", dpi='figure', format='png')
            plot_baf(het_allele_df, sample_id, gmm)
            k = gmm.k
            fit_score = gmm.score
            boundary_prob = gmm.prob

        ## if not then produce an empty JSON
        else:
            print("Fewer than 20 heterogeneous alleles present. Exiting...")
            k = None
            fit_score = None
            boundary_prob = None
            gen_empty_plot(f"./{sample_id}.bafscore.png")
            gen_empty_plot(f"./{sample_id}_allelefreq.png")

        resultsdict = { "sample_id" : sample_id,
                        "k" : k,
                        "fit_score" : fit_score,
                        "boundary_prob" : boundary_prob,
                        "fws" : fws_df["fws"][sample_id],
                        "fit_fig_path" : f"./{sample_id}.bafscore.png",
                        "cluster_fig_path" : f"./{sample_id}_allelefreq.png" }

        with open(f"{sample_id}.moi.json", "w") as fout:
            json.dump(resultsdict, fout, indent=4)


if __name__=="__main__":
    vcf_file = sys.argv[1]
    fws_csv = sys.argv[2]
    main(vcf_file, fws_csv)
