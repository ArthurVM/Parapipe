import argparse
import numpy as np
import pandas as pd
import allel
import statsmodels.api as sm
import os

# --- Argument Parsing ---
def parse_args():
    """
    Parses command-line arguments for the FWS calculation script.

    Returns:
        argparse.Namespace: An object containing the parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Compute FWS (within-host F statistics of sample heterozygosity) from a VCF file."
    )
    parser.add_argument(
        "-v", "--vcf", type=str, required=True,
        help="Path to the merged VCF file."
    )
    args = parser.parse_args()

    if not os.path.exists(args.vcf):
        raise FileNotFoundError(f"VCF file not found at: {args.vcf}")

    return args

# --- Core FWS Calculation Functions ---

def get_maf(callset):
    """
    Computes Minor Allele Frequency (MAF) for each variant.

    Args:
        callset (allel.model.CallSet): An allel CallSet object containing
                                        genomic data, specifically the 'AD' (Allelic Depth) field.

    Returns:
        numpy.ndarray: A 1D array of MAF values for each variant.
    """
    # Check if 'AD' (Allelic Depth) field is present in the callset
    if 'variants/AD' not in callset:
        raise ValueError("CallSet must contain 'variants/AD' (Allelic Depth) field to compute MAF.")

    # Access the allelic depths (AD) array
    # The AD array typically has dimensions (n_variants, n_samples, n_alleles)
    ad = callset['variants/AD']

    # Sum reference allele depths across all samples for each variant
    ref_depth = ad[:, :, 0].sum(axis=1)
    # Sum alternate allele depths across all samples for each variant
    alt_depth = ad[:, :, 1].sum(axis=1)

    # Calculate total depth for each variant
    total_depth = ref_depth + alt_depth

    # Calculate reference and alternate allele frequencies
    # Use np.divide with a 'where' condition to avoid division by zero
    # If total_depth is 0, the frequency will be NaN
    with np.errstate(divide='ignore', invalid='ignore'):
        ref_freq = np.divide(ref_depth, total_depth, where=total_depth != 0)
        alt_freq = np.divide(alt_depth, total_depth, where=total_depth != 0)

    # MAF is the minimum of reference and alternate allele frequencies
    maf = np.minimum(ref_freq, alt_freq)

    return maf

def get_population_heterozygosity(maf_array):
    """
    Computes population heterozygosity based on Minor Allele Frequency (MAF).
    Formula: H_pop = 1 - (p^2 + q^2) where p=MAF and q=1-MAF.

    Args:
        maf_array (numpy.ndarray): A 1D array of MAF values.

    Returns:
        numpy.ndarray: A 1D array of population heterozygosity values.
    """
    # Calculate heterozygosity using the formula 1 - (p^2 + q^2)
    # where p is MAF and q is (1 - MAF)
    population_het = 1 - (maf_array**2 + (1 - maf_array)**2)
    return population_het

def get_sample_heterozygosity(callset):
    """
    Computes heterozygosity for each variant within each sample.

    Args:
        callset (allel.model.CallSet): An allel CallSet object containing
                                        genomic data, specifically the 'AD' (Allelic Depth) field.

    Returns:
        numpy.ndarray: A 2D array (n_samples, n_variants) of heterozygosity values,
                       with NaN where depth is zero.
    """
    # Check if 'AD' (Allelic Depth) field is present in the callset
    if 'variants/AD' not in callset:
        raise ValueError("CallSet must contain 'variants/AD' (Allelic Depth) field to compute sample heterozygosity.")

    # Access the allelic depths (AD) array
    # Dimensions: (n_variants, n_samples, n_alleles)
    ad = callset['variants/AD']

    # Transpose AD to (n_samples, n_variants, n_alleles) for easier per-sample calculation
    sample_ad = ad.transpose(1, 0, 2)

    # Calculate total depth per sample per variant
    depth_per_sample_variant = sample_ad[:, :, 0] + sample_ad[:, :, 1]

    # Calculate allele frequencies (p and q) for each sample at each variant
    # Use np.divide with a 'where' condition to avoid division by zero, resulting in NaN
    with np.errstate(divide='ignore', invalid='ignore'):
        p_per_sample_variant = np.divide(sample_ad[:, :, 0], depth_per_sample_variant,
                                         where=depth_per_sample_variant != 0)
        q_per_sample_variant = np.divide(sample_ad[:, :, 1], depth_per_sample_variant,
                                         where=depth_per_sample_variant != 0)

    # Calculate heterozygosity for each sample at each variant
    # H_sample = 1 - (p^2 + q^2)
    het_per_sample_variant = 1 - (p_per_sample_variant**2 + q_per_sample_variant**2)

    return het_per_sample_variant

def get_fws(callset):
    """
    Computes FWS (within-host F statistics of sample heterozygosity).
    This involves:
    1. Calculating MAF and population/sample heterozygosity.
    2. Binning variants by MAF.
    3. Calculating mean heterozygosity within each MAF bin for both population and samples.
    4. Performing a linear regression of sample heterozygosity against population heterozygosity
       (without intercept) for each sample.
    5. FWS is then 1 - beta, where beta is the regression coefficient.

    Args:
        callset (allel.model.CallSet): An allel CallSet object.

    Returns:
        dict: A dictionary where keys are sample IDs and values are their FWS statistics.
    """
    sample_ids = callset['samples']
    variant_ids = callset['variants/ID'].astype(str) # Ensure variant IDs are strings

    # 1. Get MAF, population heterozygosity, and sample heterozygosity
    maf = get_maf(callset)
    population_het = get_population_heterozygosity(maf)
    sample_het = get_sample_heterozygosity(callset) # (n_samples, n_variants)

    # Convert to pandas DataFrame for easier indexing and operations, matching R output structure
    # Samples as index, variants as columns
    sample_het_df = pd.DataFrame(sample_het, index=sample_ids, columns=variant_ids)
    population_het_series = pd.Series(population_het, index=variant_ids)

    # 2. Bin variants by MAF
    # Replicate R's seq(0, 0.5, length.out = 11) for bin edges
    # This creates 10 bins: [0, 0.05), [0.05, 0.1), ..., [0.45, 0.5]
    bins = np.linspace(0, 0.5, num=11)
    
    # Use pd.cut to assign each MAF value to a bin interval
    # The 'right=True' argument is default for pd.cut and matches R's findInterval behavior
    # For MAF 0.5, it falls into the last bin, which is correct
    maf_bins = pd.cut(maf, bins=bins, include_lowest=True, labels=False, right=True)
    # The R code uses labels=FALSE, so we want integer labels (0 to 9) for bins
    # Ensure variants with MAF > 0.5 (which is rare for true MAF but good to handle)
    # or NaN MAF are treated as NaN for binning, so they don't affect means.
    maf_bins = pd.Series(maf_bins, index=variant_ids)


    # 3. Calculate mean heterozygosity within each MAF bin
    # Mean population heterozygosity per bin
    # Create a DataFrame for population_het and its corresponding bins
    pop_het_df = pd.DataFrame({'heterozygosity': population_het_series, 'maf_bin': maf_bins})
    mu_population_het = pop_het_df.groupby('maf_bin')['heterozygosity'].mean().dropna()

    # Mean sample heterozygosity per bin for each sample
    mu_sample_het = {}
    for sample_id in sample_ids:
        # Create a DataFrame for current sample's heterozygosity and corresponding bins
        sample_data = pd.DataFrame({
            'heterozygosity': sample_het_df.loc[sample_id],
            'maf_bin': maf_bins
        })
        # Group by MAF bin and calculate mean, then drop NaNs (bins with no data for this sample)
        mu_sample_het[sample_id] = sample_data.groupby('maf_bin')['heterozygosity'].mean().dropna()

    # 4. Perform linear regression and calculate FWS for each sample
    fws_results = {}
    for sample_id in sample_ids:
        # Align the sample's mean heterozygosity with the population's mean heterozygosity by bin
        # Only consider bins where both population and sample have data
        aligned_pop_het = mu_population_het.align(mu_sample_het[sample_id], join='inner')[0]
        aligned_sample_het = mu_population_het.align(mu_sample_het[sample_id], join='inner')[1]

        # Ensure there's enough data for regression (at least 2 points)
        if len(aligned_pop_het) < 2:
            print(f"Warning: Not enough data points for regression for sample {sample_id}. Setting FWS to NaN.")
            fws_results[sample_id] = np.nan
            continue
        
        # Prepare data for OLS regression
        # X is the independent variable (mu.population.het)
        # y is the dependent variable (mu.sample.het for the current sample)
        X = aligned_pop_het.values
        y = aligned_sample_het.values

        # Perform Ordinary Least Squares (OLS) regression without an intercept (as in R's lm(y ~ x - 1))
        # The coefficient of X will be our beta (coeff)
        try:
            model = sm.OLS(y, X) # y is dependent, X is independent
            results = model.fit()
            beta = results.params[0] # Get the single coefficient
            
            # FWS is 1 - beta
            fws = 1 - beta
            fws_results[sample_id] = fws
        except Exception as e:
            print(f"Error during regression for sample {sample_id}: {e}. Setting FWS to NaN.")
            fws_results[sample_id] = np.nan

    return fws_results

# --- Main execution block ---
def main():
    """
    Main function to run the FWS calculation script.
    It parses arguments, reads the VCF, computes FWS, and saves the results to a CSV file.
    """
    opt = parse_args()
    vcf_path = opt.vcf

    print(f"Reading VCF file: {vcf_path}")
    # Read the VCF file using allel.read_vcf.
    # We explicitly request 'AD' (Allelic Depth) field.
    # The `fields='*'` will load all fields, but specifically requesting 'AD' ensures it's available.
    try:
        callset = allel.read_vcf(vcf_path, fields=['samples', 'variants/ID', 'variants/AD'])
        if callset is None:
            raise ValueError(f"Could not read VCF file or it is empty: {vcf_path}")
        if 'variants/AD' not in callset or callset['variants/AD'] is None:
            raise ValueError(f"The VCF file '{vcf_path}' does not contain 'AD' (Allelic Depth) information in the 'FORMAT' field, which is required.")
        if 'samples' not in callset or callset['samples'] is None:
             raise ValueError(f"The VCF file '{vcf_path}' does not contain 'samples' information, which is required.")

    except Exception as e:
        print(f"Error reading VCF file with allel: {e}")
        print("Please ensure the VCF file is valid and contains 'AD' (Allelic Depth) information in the FORMAT field.")
        return

    print("Computing FWS...")
    fws_values = get_fws(callset)

    # Convert the FWS results to a pandas DataFrame and save to CSV
    fws_df = pd.DataFrame.from_dict(fws_values, orient='index', columns=['FWS'])
    fws_df.index.name = 'Sample_ID'

    output_filename = "fws.csv"
    fws_df.to_csv(output_filename)
    print(f"FWS results saved to '{output_filename}'")
    print("\nCalculation complete.")

if __name__ == "__main__":
    main()
