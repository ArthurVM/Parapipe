
import numpy as np
import pandas as pd
from scipy.stats import linregress
from Variant import Variant


class VCF:
    """
    Represents a VCF file with methods for parsing and accessing sample-level allele depth and frequency data.
    """

    def __init__(self, vcf_file: str):
        self.vcf_file = vcf_file
        self.info = []
        self.var_lines = []
        self.sample_ids = []
        self._parse_vcf()
        
        print(f"{len(self.sample_ids)} samples detected.")
        print(self.sample_ids)

    def _parse_vcf(self):
        with open(self.vcf_file, "r") as fin:
            for line in fin:
                if line.startswith("#CHROM"):
                    self.sample_ids = line.strip().split("\t")[9:]
                elif line.startswith("#"):
                    self.info.append(line)
                else:
                    variant = Variant()
                    variant.parseLine(line, self.sample_ids)
                    self.var_lines.append(variant)

    def getVars(self):
        return self.var_lines

    def getSampleIds(self):
        return self.sample_ids

    def getAlleleCounts(self, sample_id: str, filter_homo: bool = False) -> pd.DataFrame:
        """
        Outputs ref/alt depths and allele frequencies for the given sample.
        """
        allele_list = []

        for variant in self.var_lines:
            if sample_id not in variant.samples:
                raise ValueError(f"{sample_id} not found in variant.")

            ad = variant.samples[sample_id].get("AD", None)
            if ad is None or ad[0] in (".", None):
                continue

            try:
                r_depth = int(ad[0])
                a_depth = sum(int(x) for x in ad[1:] if x not in (".", None))
            except (ValueError, TypeError):
                continue

            total = r_depth + a_depth
            if total == 0:
                continue

            r_freq = r_depth / total
            a_freq = a_depth / total

            if r_freq < 1.0:
                if filter_homo and r_freq == 0.0:
                    continue
                allele_list.append({
                    "chromosome": variant.chrom,
                    "position": variant.pos,
                    "ref_freq": r_freq,
                    "alt_freq": a_freq,
                    "ref_depth": r_depth,
                    "alt_depth": a_depth
                })

        return pd.DataFrame.from_records(allele_list)

    def getMAF(self, filter_homo: bool = False) -> pd.DataFrame:
        """
        Returns sample-wise minor allele frequency matrix.
        """
        maf_matrix = {}

        for sample_id in self.sample_ids:
            mafs = []
            for variant in self.var_lines:
                ad = variant.samples.get(sample_id, {}).get("AD", None)
                if ad is None or ad[0] in (".", None):
                    mafs.append(np.nan)
                    continue

                try:
                    r_depth = int(ad[0])
                    a_depth = sum(int(x) for x in ad[1:] if x not in (".", None))
                except (ValueError, TypeError):
                    mafs.append(np.nan)
                    continue

                total = r_depth + a_depth
                if total == 0:
                    mafs.append(np.nan)
                else:
                    maf = min(r_depth / total, a_depth / total)
                    mafs.append(maf)

            maf_matrix[sample_id] = mafs

        return pd.DataFrame(maf_matrix)

    def getPopMAF(self, filter_homo: bool = False) -> np.ndarray:
        """
        Returns a population-level MAF vector across all samples.
        """
        maf_list = []

        for variant in self.var_lines:
            total = 0
            mac = 0
            for sample_id in self.sample_ids:
                ad = variant.samples.get(sample_id, {}).get("AD", None)
                if ad is None or ad[0] in (".", None):
                    continue

                try:
                    r_depth = int(ad[0])
                    a_depth = sum(int(x) for x in ad[1:] if x not in (".", None))
                except (ValueError, TypeError):
                    continue

                total += r_depth + a_depth
                mac += min(r_depth, a_depth)

            if total > 0:
                maf_list.append(mac / total)
            else:
                maf_list.append(np.nan)

        return np.array(maf_list)

    def fws(self, maf_matrix: np.ndarray, n_maf_bins: int = 10) -> dict:
        """
        Computes Fws statistic per sample.
        """
        if maf_matrix.ndim != 2:
            raise ValueError("Input to fws() must be 2D. Use getMAF(), not getPopMAF().")

        fws_values = {}

        h_matrix = 2 * maf_matrix * (1 - maf_matrix)
        h_s = np.nanmean(h_matrix, axis=0)

        maf_bins = np.linspace(0, 0.5, n_maf_bins + 1)

        for i, sample_id in enumerate(self.sample_ids):
            binwise_fws = []

            for j in range(n_maf_bins):
                lower, upper = maf_bins[j], maf_bins[j + 1]
                idx = np.where((maf_matrix[i] >= lower) & (maf_matrix[i] < upper))[0]

                if len(idx) > 0:
                    h_w = np.nanmean(h_matrix[i, idx])
                    h_s_bin = np.nanmean(h_s[idx])
                    if h_s_bin > 0:
                        ratio = h_w / h_s_bin
                        ratio = min(ratio, 1.0)  # Clamp to avoid negative Fws
                        binwise_fws.append(1 - ratio)

            fws_values[sample_id] = np.nanmean(binwise_fws) if binwise_fws else np.nan

        return fws_values