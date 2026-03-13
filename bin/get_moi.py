import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from VCF import VCF  # noqa: E402
from plotting_funcs import load_chrom_lengths_from_fasta, plot_maf_across_genome  # noqa: E402


def gen_empty_plot(save_path: Path, text: str = "MAF plot placeholder") -> None:
    """
    Produce an empty plot placeholder to satisfy downstream report wiring.
    """
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=16, color="red")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.savefig(save_path, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def main(vcf_file: str, fws_csv_out: str | None = None, reference_fasta: str | None = None) -> None:
    """
    Combined Fws/MOI placeholder:
      - computes Fws per sample directly from the VCF
      - writes fws.csv (optional path)
      - emits per-sample JSON with sample_id, fws, maf_plot
      - generates per-sample MAF scatter plots across the genome (optionally scaled by reference lengths)
    """
    
    chrom_lengths = load_chrom_lengths_from_fasta(reference_fasta) if reference_fasta else None

    print("Reading vcf file... ", end="", flush=True)
    vcf_obj = VCF(vcf_file)
    print(f"Found {len(vcf_obj.sample_ids)} samples.", flush=True)
    print("Computing Fws... ", end="", flush=True)
    fws_values, hs_values = vcf_obj.fws()
    print("Done.", flush=True)

    # Write Fws table for downstream use
    fws_df = pd.DataFrame({"fws": fws_values, "hs": hs_values})
    fws_df.index.name = "sample_id"
    fws_out_path = Path(fws_csv_out) if fws_csv_out else Path("fws.csv")
    fws_df.to_csv(fws_out_path)

    for sample_id in vcf_obj.getSampleIds():
        maf_plot_path = Path(f"{sample_id}_maf_plot.png")
        plot_maf_across_genome(vcf_obj, sample_id, maf_plot_path, chrom_lengths=chrom_lengths)

        fws_value = fws_values.get(sample_id)
        hs_value = hs_values.get(sample_id)
        if pd.isna(fws_value):
            fws_value = None
        if pd.isna(hs_value):
            hs_value = None

        resultsdict = {
            "sample_id": sample_id,
            "fws": fws_value,
            "hs": hs_value,
            "maf_plot": str(maf_plot_path),
        }

        with open(f"{sample_id}.moi.json", "w") as fout:
            json.dump(resultsdict, fout, indent=4)


if __name__ == "__main__":
    vcf_file = sys.argv[1]
    fws_csv_out = sys.argv[2] if len(sys.argv) > 2 else None
    reference_fasta = sys.argv[3] if len(sys.argv) > 3 else None
    main(vcf_file, fws_csv_out, reference_fasta)
