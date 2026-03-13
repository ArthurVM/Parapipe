import os, sys
import numpy as np
import pandas as pd
import subprocess
import json
import math
import re
from Bio.Seq import Seq
from collections import defaultdict, Counter

FAMILY_FILTERS = {
    "IIa": {"min_cov": 3, "min_af": 0.05, "l3": 0.15, "l6": 0.075},
    "IIb": {"min_cov": 3, "min_af": 0.05, "l3": 0.15, "l6": 0.075},
    "IIc": {"min_cov": 3, "min_af": 0.05, "l3": 0.15, "l6": 0.075},
    "IId": {"min_cov": 3, "min_af": 0.05, "l3": 0.15, "l6": 0.075},
    "IIr": {"min_cov": 3, "min_af": 0.05, "l3": 0.15, "l6": 0.075},
}

def get_family_filters(family):
    ## individual filtering parameters arent really necessary as far as I know
    ## but keep as a placeholder
    # return FAMILY_FILTERS.get(family, FAMILY_FILTERS["IIa"])
    return {"min_cov": 3, "min_af": 0.05, "l3": 0.15, "l6": 0.075}

def parseReport(report):

    prefix = report.split("/")[-1].split(".")[0]

    seqtups = []
    fragtups = []
    with open(report, 'r') as fin:
        lines = fin.read().split("Sequence variants:")[-1].split("\nLength variants:")[0]
        for line in lines.split("\n"):
            if line == "":
                continue
            seqtup, count = line.split("\t")
            _, seq, _, family, _ = seqtup.split("'")
            count = int(count)

            fullseq = str(Seq(seq).reverse_complement())
            
            subtype = inferGP60Subtype(fullseq, family)

            seqtups.append((fullseq, subtype, count))
    return sorted(seqtups, key=lambda x: x[1], reverse=True)


def filterAlleles(sid, seqtups, min_cov=5, min_af=0.05, stutter_filter=True, l3=0.16, l6=0.08):
    if not seqtups:
        return []
    
    filtered_seqtups = []
    s_alleles = sorted(seqtups, key=lambda x: x[-1], reverse=True)
    major_allele = s_alleles[0]
    
    print(s_alleles, major_allele)

    ## calculate total coverage to determine allele frequencies
    total_cov = sum(count for _, _, count in seqtups)
    if total_cov == 0:
        return []

    primary_cov = major_allele[-1]
    lmaj = len(major_allele[0])
    
    maj_fam = major_allele[1][-3:]

    ## the major allele is always kept, assuming it represents the primary pop
    filtered_seqtups.append(major_allele)

    ## if minor alleles exist, filter them
    if len(s_alleles) > 1:
        minor_alleles = s_alleles[1:]

        for seq, subtype, count in minor_alleles:
            poly_family = subtype[-3:] != maj_fam
            if poly_family:
                print(f"poly family: {sid} ({maj_fam}) {subtype[-3:]} {count}")

            ## GENERAL FILTERS
            ##  filter by minimum coverage
            if count < min_cov:
                continue

            ## filter by minimum allele frequency
            af = count / total_cov
            if af < min_af:
                continue

            ## STUTTER FILTERS
            if stutter_filter:
                lmin = len(seq)
                is_stutter_3 = (lmin == lmaj - 3)
                is_stutter_6 = (lmin == lmaj - 6)

                if is_stutter_3 or is_stutter_6:
                    ## calculate frequency relative to the major allele for stutter check
                    primary_af = count / primary_cov
                    if is_stutter_3 and primary_af < l3:
                        continue  ## filtered as stutter artifact
                    if is_stutter_6 and primary_af < l6:
                        continue  ## filtered as stutter artifact

            ## if the allele has passed all applicable filters then add it to the list
            filtered_seqtups.append((seq, subtype, count))

    return filtered_seqtups


def inferGP60Subtype(seq, family):
    """
    Infer gp60 subtype string using CryptoGenotyper repeat/motif rules.
    """
    if not seq:
        return family

    seq = seq.upper().replace(" ", "")
    family_core = family.split("|")[0].split("(")[0]

    ## legacy probe quirk for IId: drop the trailing base to keep frame.
    if family_core == "IId" and len(seq) > 1:
        seq = seq[:-1]

    repeat_match = re.match(r"(?:TCA|TCG|TCT)+", seq)
    repeat_region = repeat_match.group(0) if repeat_match else ""
    after_repeat = seq[len(repeat_region):] if repeat_region else ""
    if not repeat_region:
        repeat_region, after_repeat = seq, ""

    codons = [
        repeat_region[i : i + 3]
        for i in range(0, len(repeat_region) - len(repeat_region) % 3, 3)
    ]
    codon_counts = Counter(codons)

    subtype_parts = []
    good_repeat = True

    num_tca = codon_counts.get("TCA", 0)
    num_tcg = codon_counts.get("TCG", 0)
    num_tct = codon_counts.get("TCT", 0)

    if num_tca:
        subtype_parts.append(f"A{num_tca}")
    if num_tcg:
        subtype_parts.append(f"G{num_tcg}")

    allows_t = (
        family_core.startswith("Ib")
        or family_core.startswith("Ie")
        or family_core.startswith("IV")
        or family_core == "VIIa"
        or family_core == "XIa"
        or "XIV" in family_core
    )
    if num_tct:
        if allows_t:
            subtype_parts.append(f"T{num_tct}")
        elif family_core.startswith("XVIIa") and re.search(
            r"GGTGTTACCACTGCTCCTGTGGCA", after_repeat
        ):
            subtype_parts.append("R1")
        else:
            good_repeat = False

    num_acatca = after_repeat.count("ACATCA")
    if num_acatca and (
        family_core in {"IIa", "IIl", "IIt"} or family_core.startswith("XIII")
    ):
        subtype_parts.append(f"R{num_acatca}")

    if family_core.startswith("Ia"):
        hom_r = sum(
            after_repeat.count(p)
            for p in (
                "AAAACGGTGGTAAGG",
                "AGAACGGTGGTAAGG",
                "AAGACGGTGGTAAGG",
                "AGGACGGTGGTAAGG",
            )
        )
        if hom_r >= 2:
            hom_r += sum(
                after_repeat.count(p)
                for p in (
                    "AAAACGGTGAAGG",
                    "AAGACGGTGAAGG",
                    "AGAACGGTGAAGG",
                    "AGGACGGTGAAGG",
                )
            )
        if hom_r:
            subtype_parts.append(f"R{hom_r}")

    if family_core.startswith("If"):
        homf_r = 0
        homf_r += len(re.findall(r"AAGAAGGCAAGAGAAG", after_repeat))
        homf_r += len(re.findall(r".AGA.GGCA..GAAG", after_repeat))
        if homf_r:
            subtype_parts.append(f"R{homf_r}")

    if family_core.startswith("Vb"):
        m = re.search(r"(ACA)+", after_repeat)
        if m:
            subtype_parts.append(f"R{m.group(0).count('ACA')}")

    if family_core.startswith("XXV"):
        num_suis = len(re.findall(r"GGTG.{1}TCAAG.{1}GAATGC.{1}CAG", seq))
        subtype_parts = [f"R{num_suis}"] if num_suis else []
    elif family_core.startswith("XXIV"):
        subtype_parts = []

    if not good_repeat or re.match(
        r"XX[a-z]{1}|XXI[a-z]{1}|XXII[a-z]{1}|XXIII[a-z]{1}|XXVI[a-z]{1}", family_core
    ):
        subtype_parts = []

    return family_core + "".join(subtype_parts)
    

def filter_only_mixed(moi_dict):

    filtered_dict = defaultdict(dict)

    for sample, allele_dict in moi_dict.items():
        if len(allele_dict) > 1:
            filtered_dict[sample] = allele_dict
    
    return filtered_dict


def load_polyfamily_json(poly_json):
    """Allow both {family: [[seq, count], ...]} and {sample: {family: ...}} layouts."""
    with open(poly_json, "r") as fin:
        jdata = json.load(fin)

    if not jdata:
        return {}

    # Multi-sample shape: {sample: {family: [[seq, count], ...]}}
    first_val = next(iter(jdata.values()))
    if isinstance(first_val, dict):
        return jdata

    # Single-sample shape: {family: [[seq, count], ...]}
    sample_id = os.path.basename(poly_json).split(".")[0]
    return {sample_id: jdata}




def filter_family_groups(family_dict, min_reads=3, min_af=0.05):
    if not family_dict:
        return {}

    if len(family_dict) <= 1:
        return family_dict

    totals = {}
    for family, entries in family_dict.items():
        total = 0
        if entries:
            for entry in entries:
                if not entry or len(entry) < 2:
                    continue
                try:
                    total += int(entry[1])
                except (TypeError, ValueError):
                    continue
        totals[family] = total

    overall = sum(totals.values())
    if overall == 0:
        return {}

    filtered = {}
    for family, entries in family_dict.items():
        fam_total = totals.get(family, 0)
        if fam_total < min_reads:
            continue
        if (fam_total / overall) < min_af:
            continue
        filtered[family] = entries

    return filtered


def process_family(sample_id, family, entries):
    """Name gp60 alleles for a single family and apply family-specific artifact filters."""
    seqtups = []
    for entry in entries:
        if not entry or len(entry) < 2:
            continue
        seq, count = entry[0], entry[1]
        fullseq = str(Seq(seq))
        subtype = inferGP60Subtype(fullseq, family) or f"UNK-{family}"
        seqtups.append((fullseq, subtype, count))

    filt_cfg = get_family_filters(family)
    return filterAlleles(
        f"{sample_id}:{family}",
        seqtups,
        min_af=filt_cfg["min_af"],
        min_cov=filt_cfg["min_cov"],
        stutter_filter=True,
        l3=filt_cfg["l3"],
        l6=filt_cfg["l6"],
    )


def run(poly_json, outdir=None):

    data = load_polyfamily_json(poly_json)
    results = {}

    for sid, family_dict in data.items():
        family_dict = filter_family_groups(family_dict)
        fam_results = {}
        for family, entries in family_dict.items():
            filtered = process_family(sid, family, entries)
            if filtered:
                fam_results[family] = filtered
        if fam_results:
            results[sid] = fam_results

    outdir = outdir or os.path.dirname(poly_json) or "."
    prefix = os.path.basename(poly_json).rsplit(".json", 1)[0]
    out_path = os.path.join(outdir, f"{prefix}.filtered.json")

    with open(out_path, "w") as f:
        json.dump(results, f, indent=4)

    print(f"Wrote filtered gp60 calls with subtypes to {out_path}")


def main(poly_json, outdir=None):
    run(poly_json, outdir=outdir)


if __name__ == "__main__":
    poly_json = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else None
    main(poly_json, outdir=outdir)
