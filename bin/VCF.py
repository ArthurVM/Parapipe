import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

try:
    import zarr
    from numcodecs import Blosc, VLenUTF8
except ImportError:  # pragma: no cover - optional dependency
    zarr = None  # type: ignore
    Blosc = None  # type: ignore
    VLenUTF8 = None  # type: ignore

from Variant import Variant


class VCF:
    """
    Memory-aware VCF reader with optional Zarr backing store for allele depth data.
    """

    def __init__(
        self,
        vcf_file: Optional[str] = None,
        *,
        zarr_store: Optional[str] = None,
        load_into_memory: bool = False,
    ):
        if vcf_file is None and zarr_store is None:
            raise ValueError("Provide either a VCF file path or a Zarr store.")

        self.vcf_file = vcf_file
        self.zarr_store = zarr_store
        self.load_into_memory = load_into_memory

        self.info: List[str] = []
        self.sample_ids: List[str] = []
        self.var_lines: Optional[List[Variant]] = [] if load_into_memory else None
        self._zarr_group = None

        if zarr_store and Path(zarr_store).exists():
            self._open_zarr(zarr_store)

        if vcf_file:
            self._parse_header()
            if load_into_memory:
                self._load_variants_into_memory()
        elif not self.sample_ids:
            raise ValueError("Unable to determine sample IDs from the provided inputs.")

    @classmethod
    def from_zarr(cls, store_path: str) -> "VCF":
        """
        Construct a VCF instance backed by an existing Zarr store.
        """
        return cls(vcf_file=None, zarr_store=store_path, load_into_memory=False)

    @staticmethod
    def _coerce_int(value) -> Optional[int]:
        if value in (None, ".", ""):
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            try:
                return int(float(value))
            except (TypeError, ValueError):
                return None

    def _open_zarr(self, store_path: str, mode: str = "r") -> None:
        if zarr is None:
            raise ImportError("Install zarr to read or write Zarr stores.")

        self._zarr_group = zarr.open_group(store_path, mode=mode)
        self.zarr_store = store_path

        sample_ids = self._zarr_group.attrs.get("sample_ids")
        if sample_ids is None:
            raise ValueError("Zarr store is missing required 'sample_ids' attribute.")

        self.sample_ids = list(sample_ids)
        info_lines = self._zarr_group.attrs.get("info")
        if info_lines:
            self.info = list(info_lines)

    def _parse_header(self) -> None:
        if not self.vcf_file:
            return

        with open(self.vcf_file, "r") as fin:
            for line in fin:
                if line.startswith("#CHROM"):
                    self.sample_ids = line.strip().split("\t")[9:]
                    break
                if line.startswith("#"):
                    self.info.append(line)

        if not self.sample_ids:
            raise ValueError("VCF header missing '#CHROM' line with sample identifiers.")

    def _load_variants_into_memory(self) -> None:
        if self.var_lines is None:
            self.var_lines = []

        for variant in self._iter_vcf_records():
            self.var_lines.append(variant)

    def _iter_vcf_records(self) -> Iterator[Variant]:
        if not self.vcf_file:
            raise ValueError("No VCF file available for iteration.")

        with open(self.vcf_file, "r") as fin:
            for line in fin:
                if line.startswith("#"):
                    continue
                variant = Variant()
                variant.parseLine(line, self.sample_ids)
                yield variant

    def _extract_depths(self, variant: Variant) -> Tuple[np.ndarray, np.ndarray]:
        n_samples = len(self.sample_ids)
        ref_depths = np.full(n_samples, np.nan, dtype=float)
        alt_depths = np.full(n_samples, np.nan, dtype=float)

        for idx, sample_id in enumerate(self.sample_ids):
            sample_data = variant.samples.get(sample_id)
            if not sample_data:
                continue

            ad = sample_data.get("AD")
            if ad is None:
                continue

            if isinstance(ad, (list, tuple)):
                ref = self._coerce_int(ad[0]) if len(ad) > 0 else None
                alt_values = [self._coerce_int(x) for x in ad[1:]]
            else:
                ref = self._coerce_int(ad)
                alt_values: List[Optional[int]] = []

            if ref is None and not alt_values:
                continue

            ref_val = float(ref) if ref is not None else 0.0
            alt_val = float(sum(x for x in alt_values if x is not None))

            if ref_val == 0.0 and alt_val == 0.0:
                continue

            ref_depths[idx] = ref_val
            alt_depths[idx] = alt_val

        return ref_depths, alt_depths

    def _iter_vcf_depths(self) -> Iterator[Tuple[str, int, np.ndarray, np.ndarray]]:
        for variant in self._iter_vcf_records():
            ref_depths, alt_depths = self._extract_depths(variant)
            yield variant.chrom, variant.pos, ref_depths, alt_depths

    def _iter_memory_depths(self) -> Iterator[Tuple[str, int, np.ndarray, np.ndarray]]:
        if not self.var_lines:
            return

        for variant in self.var_lines:
            ref_depths, alt_depths = self._extract_depths(variant)
            yield variant.chrom, variant.pos, ref_depths, alt_depths

    def _iter_zarr_depths(self) -> Iterator[Tuple[str, int, np.ndarray, np.ndarray]]:
        if self._zarr_group is None:
            raise RuntimeError("Zarr group is not initialised.")

        ref_arr = self._zarr_group["ref_depth"]
        alt_arr = self._zarr_group["alt_depth"]
        chrom_arr = self._zarr_group.get("chrom")
        pos_arr = self._zarr_group.get("pos")

        n_variants = ref_arr.shape[0]
        for idx in range(n_variants):
            ref_depths = np.asarray(ref_arr[idx, :], dtype=float)
            alt_depths = np.asarray(alt_arr[idx, :], dtype=float)

            chrom = chrom_arr[idx] if chrom_arr is not None else ""
            if isinstance(chrom, bytes):
                chrom = chrom.decode("utf-8")
            pos = int(pos_arr[idx]) if pos_arr is not None else idx

            yield chrom, pos, ref_depths, alt_depths

    def _iter_allele_depths(self) -> Iterator[Tuple[str, int, np.ndarray, np.ndarray]]:
        if self._zarr_group is not None:
            yield from self._iter_zarr_depths()
            return

        if self.var_lines:
            yield from self._iter_memory_depths()
            return

        yield from self._iter_vcf_depths()

    def getVars(self):
        if self.var_lines is None:
            raise RuntimeError("Variants not loaded into memory; set load_into_memory=True.")
        return self.var_lines

    def getSampleIds(self):
        return self.sample_ids

    def getAlleleCounts(self, sample_id: str, filter_homo: bool = False) -> pd.DataFrame:
        """
        Stream ref/alt depths and allele frequencies for the given sample.
        """
        if sample_id not in self.sample_ids:
            raise ValueError(f"{sample_id} not present in VCF samples.")

        sample_idx = self.sample_ids.index(sample_id)
        allele_list = []

        for chrom, pos, ref_depths, alt_depths in self._iter_allele_depths():
            ref = ref_depths[sample_idx]
            alt = alt_depths[sample_idx]

            if np.isnan(ref) or np.isnan(alt):
                continue

            total = ref + alt
            if total <= 0:
                continue

            r_freq = ref / total
            a_freq = alt / total

            if r_freq < 1.0:
                if filter_homo and r_freq == 0.0:
                    continue
                allele_list.append(
                    {
                        "chromosome": chrom,
                        "position": pos,
                        "ref_freq": r_freq,
                        "alt_freq": a_freq,
                        "ref_depth": ref,
                        "alt_depth": alt,
                    }
                )

        return pd.DataFrame.from_records(allele_list)

    def getMAF(self, filter_homo: bool = False) -> pd.DataFrame:
        """
        Returns sample-wise minor allele frequency matrix.
        Note: constructing this matrix requires holding all values in memory.
        """
        maf_matrix: Dict[str, List[float]] = {sid: [] for sid in self.sample_ids}

        for _, _, ref_depths, alt_depths in self._iter_allele_depths():
            totals = ref_depths + alt_depths
            minors = np.minimum(ref_depths, alt_depths)

            valid = (~np.isnan(totals)) & (~np.isnan(minors)) & (totals > 0)
            maf_row = np.full(len(self.sample_ids), np.nan, dtype=float)
            maf_row[valid] = minors[valid] / totals[valid]

            for idx, sample_id in enumerate(self.sample_ids):
                maf_matrix[sample_id].append(float(maf_row[idx]))

        return pd.DataFrame(maf_matrix)

    def getPopMAF(self, filter_homo: bool = False) -> np.ndarray:
        """
        Returns a population-level MAF vector across all samples using streaming aggregation.
        """
        maf_list: List[float] = []

        for _, _, ref_depths, alt_depths in self._iter_allele_depths():
            total_ref = np.nansum(ref_depths)
            total_alt = np.nansum(alt_depths)
            total_depth = total_ref + total_alt

            if total_depth > 0:
                maf_list.append(min(total_ref, total_alt) / total_depth)
            else:
                maf_list.append(np.nan)

        return np.array(maf_list, dtype=float)

    def fws(
        self,
        maf_matrix: Optional[np.ndarray | pd.DataFrame] = None,
        n_maf_bins: int = 10,
    ) -> tuple[dict, dict]:
        """
        Computes the Fws statistic per sample. Defaults to a streaming implementation
        that avoids materialising the full variant-by-sample matrix in memory.
        Returns a tuple of (fws_per_sample, hs_per_sample).
        """
        if maf_matrix is not None:
            return self._fws_from_matrix(maf_matrix, n_maf_bins)

        return self._fws_streaming(n_maf_bins)

    def _fws_streaming(self, n_maf_bins: int) -> dict:
        n_samples = len(self.sample_ids)
        maf_bin_edges = np.linspace(0.0, 0.5, n_maf_bins + 1)

        pop_sum = np.zeros(n_maf_bins, dtype=float)
        pop_count = np.zeros(n_maf_bins, dtype=float)
        sample_sum = np.zeros((n_maf_bins, n_samples), dtype=float)
        sample_count = np.zeros((n_maf_bins, n_samples), dtype=float)

        for _, _, ref_depths, alt_depths in self._iter_allele_depths():
            if np.all(np.isnan(ref_depths)) and np.all(np.isnan(alt_depths)):
                continue

            ref_depths = np.nan_to_num(ref_depths, nan=0.0)
            alt_depths = np.nan_to_num(alt_depths, nan=0.0)

            total_ref = ref_depths.sum()
            total_alt = alt_depths.sum()
            total_depth = total_ref + total_alt

            if total_depth <= 0:
                continue

            maf = min(total_ref, total_alt) / total_depth
            bin_idx = np.searchsorted(maf_bin_edges, maf, side="right") - 1
            if bin_idx < 0 or bin_idx >= n_maf_bins:
                continue

            pop_p = total_ref / total_depth
            pop_q = total_alt / total_depth
            pop_het = 1.0 - (pop_p**2 + pop_q**2)

            pop_sum[bin_idx] += pop_het
            pop_count[bin_idx] += 1

            sample_totals = ref_depths + alt_depths
            valid_samples = sample_totals > 0
            if not np.any(valid_samples):
                continue

            p = np.divide(ref_depths, sample_totals, out=np.zeros_like(ref_depths), where=valid_samples)
            q = np.divide(alt_depths, sample_totals, out=np.zeros_like(alt_depths), where=valid_samples)
            sample_het = 1.0 - (p**2 + q**2)

            sample_sum[bin_idx, valid_samples] += sample_het[valid_samples]
            sample_count[bin_idx, valid_samples] += 1

        mu_pop_het = np.full(n_maf_bins, np.nan, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            valid = pop_count > 0
            mu_pop_het[valid] = pop_sum[valid] / pop_count[valid]

        mu_sample_het = np.full((n_maf_bins, n_samples), np.nan, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            valid_samples = sample_count > 0
            mu_sample_het[valid_samples] = sample_sum[valid_samples] / sample_count[valid_samples]

        fws_values = {}
        hs_values = {}
        with np.errstate(divide="ignore", invalid="ignore"):
            total_sample_counts = sample_count.sum(axis=0)
            total_sample_sum = sample_sum.sum(axis=0)
            hs_means = np.divide(
                total_sample_sum,
                total_sample_counts,
                out=np.full(n_samples, np.nan, dtype=float),
                where=total_sample_counts > 0,
            )

        for sample_idx, sample_id in enumerate(self.sample_ids):
            sample_bin_means = mu_sample_het[:, sample_idx]
            valid_bins = (~np.isnan(mu_pop_het)) & (~np.isnan(sample_bin_means))

            if not np.any(valid_bins):
                fws_values[sample_id] = np.nan
                hs_values[sample_id] = hs_means[sample_idx]
                continue

            x = mu_pop_het[valid_bins]
            y = sample_bin_means[valid_bins]
            denominator = np.dot(x, x)

            if np.isclose(denominator, 0.0):
                fws_values[sample_id] = np.nan
                hs_values[sample_id] = hs_means[sample_idx]
                continue

            slope = np.dot(x, y) / denominator
            fws_values[sample_id] = 1.0 - slope
            hs_values[sample_id] = hs_means[sample_idx]

        return fws_values, hs_values

    def _fws_from_matrix(
        self,
        maf_matrix: np.ndarray | pd.DataFrame,
        n_maf_bins: int,
    ) -> tuple[dict, dict]:
        if isinstance(maf_matrix, pd.DataFrame):
            maf_df = maf_matrix.copy()
        else:
            maf_array = np.asarray(maf_matrix, dtype=float)
            if maf_array.ndim != 2:
                raise ValueError("MAF matrix must be two-dimensional.")
            if maf_array.shape[1] == len(self.sample_ids):
                maf_df = pd.DataFrame(maf_array, columns=self.sample_ids)
            elif maf_array.shape[0] == len(self.sample_ids):
                maf_df = pd.DataFrame(maf_array.T, columns=self.sample_ids)
            else:
                raise ValueError("MAF matrix dimensions do not align with sample IDs.")

        missing_samples = [sid for sid in self.sample_ids if sid not in maf_df.columns]
        if missing_samples:
            raise ValueError(f"MAF matrix is missing samples: {missing_samples}")

        maf_df = maf_df[self.sample_ids]
        maf_values = maf_df.to_numpy(dtype=float)

        pop_maf = self.getPopMAF()
        if pop_maf.shape[0] != maf_values.shape[0]:
            raise ValueError("Population MAF vector length does not match number of variants.")

        pop_het = 1 - (pop_maf**2 + (1 - pop_maf) ** 2)
        sample_het = 1 - (maf_values**2 + (1 - maf_values) ** 2)

        maf_bin_edges = np.linspace(0.0, 0.5, n_maf_bins + 1)
        bin_indices = np.full(pop_maf.shape, -1, dtype=int)
        valid_pop = ~np.isnan(pop_maf)
        bin_indices[valid_pop] = np.searchsorted(maf_bin_edges, pop_maf[valid_pop], side="right") - 1
        bin_indices = np.clip(bin_indices, -1, n_maf_bins - 1)

        mu_pop_het = np.full(n_maf_bins, np.nan, dtype=float)
        mu_sample_het = np.full((n_maf_bins, maf_values.shape[1]), np.nan, dtype=float)

        for b in range(n_maf_bins):
            mask = bin_indices == b
            if not np.any(mask):
                continue

            pop_vals = pop_het[mask]
            pop_vals = pop_vals[~np.isnan(pop_vals)]
            if pop_vals.size:
                mu_pop_het[b] = pop_vals.mean()

            sample_vals = sample_het[mask, :]
            if sample_vals.size:
                with np.errstate(invalid="ignore"):
                    means = np.nanmean(sample_vals, axis=0)
                mu_sample_het[b, :] = means

        fws_values = {}
        hs_values = {}
        with np.errstate(invalid="ignore"):
            variant_mask = bin_indices >= 0
            hs_means = np.nanmean(sample_het[variant_mask], axis=0) if np.any(variant_mask) else np.full(maf_values.shape[1], np.nan)

        for sample_idx, sample_id in enumerate(self.sample_ids):
            sample_bin_means = mu_sample_het[:, sample_idx]
            valid_bins = (~np.isnan(mu_pop_het)) & (~np.isnan(sample_bin_means))
            x = mu_pop_het[valid_bins]
            y = sample_bin_means[valid_bins]

            if x.size == 0 or np.isclose(np.dot(x, x), 0.0):
                fws_values[sample_id] = np.nan
                hs_values[sample_id] = hs_means[sample_idx]
                continue

            slope = np.dot(x, y) / np.dot(x, x)
            fws_values[sample_id] = 1 - slope
            hs_values[sample_id] = hs_means[sample_idx]

        return fws_values, hs_values

    def write_zarr(
        self,
        store_path: str,
        *,
        chunk_size: int = 1024,
        compressor=None,
    ) -> None:
        """
        Convert the underlying VCF into a Zarr store with chunked allele depth arrays.
        """
        if self.vcf_file is None:
            raise ValueError("A source VCF file is required to write a Zarr store.")
        if zarr is None:
            raise ImportError("Install zarr to enable writing to Zarr stores.")
        if VLenUTF8 is None or Blosc is None:
            raise ImportError("Install numcodecs to enable compression and UTF-8 storage.")

        store_path = str(store_path)
        group = zarr.open_group(store_path, mode="w")
        group.attrs["sample_ids"] = list(self.sample_ids)
        if self.info:
            group.attrs["info"] = list(self.info)

        n_samples = len(self.sample_ids)
        if compressor is None:
            compressor = Blosc(cname="zstd", clevel=5, shuffle=Blosc.BITSHUFFLE)

        ref_arr = group.create_dataset(
            "ref_depth",
            shape=(0, n_samples),
            chunks=(chunk_size, n_samples),
            maxshape=(None, n_samples),
            dtype="f4",
            compressor=compressor,
        )
        alt_arr = group.create_dataset(
            "alt_depth",
            shape=(0, n_samples),
            chunks=(chunk_size, n_samples),
            maxshape=(None, n_samples),
            dtype="f4",
            compressor=compressor,
        )
        chrom_arr = group.create_dataset(
            "chrom",
            shape=(0,),
            chunks=(chunk_size,),
            maxshape=(None,),
            dtype=object,
            object_codec=VLenUTF8(),
        )
        pos_arr = group.create_dataset(
            "pos",
            shape=(0,),
            chunks=(chunk_size,),
            maxshape=(None,),
            dtype="i8",
        )

        buffer_ref: List[np.ndarray] = []
        buffer_alt: List[np.ndarray] = []
        buffer_chr: List[str] = []
        buffer_pos: List[int] = []

        def flush_buffers() -> None:
            if not buffer_ref:
                return

            chunk_len = len(buffer_ref)
            start = ref_arr.shape[0]
            stop = start + chunk_len

            ref_arr.resize(stop, axis=0)
            alt_arr.resize(stop, axis=0)
            chrom_arr.resize(stop, axis=0)
            pos_arr.resize(stop, axis=0)

            ref_arr[start:stop, :] = np.vstack(buffer_ref).astype("float32", copy=False)
            alt_arr[start:stop, :] = np.vstack(buffer_alt).astype("float32", copy=False)
            chrom_arr[start:stop] = np.array(buffer_chr, dtype=object)
            pos_arr[start:stop] = np.array(buffer_pos, dtype="int64")

            buffer_ref.clear()
            buffer_alt.clear()
            buffer_chr.clear()
            buffer_pos.clear()

        for chrom, pos, ref_depths, alt_depths in self._iter_vcf_depths():
            buffer_ref.append(ref_depths.astype("float32", copy=False))
            buffer_alt.append(alt_depths.astype("float32", copy=False))
            buffer_chr.append(chrom)
            buffer_pos.append(int(pos))

            if len(buffer_ref) >= chunk_size:
                flush_buffers()

        flush_buffers()

        self._open_zarr(store_path, mode="r")
