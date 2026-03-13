#!/bin/bash

# Build all per-process containers and report success/failure for each.

set -u
set -o pipefail

# Where to store build logs (can override with LOG_FILE env)
LOG_FILE="${LOG_FILE:-singularity/build.log}"
mkdir -p "$(dirname "${LOG_FILE}")"
echo "=== Singularity build run: $(date) ===" > "${LOG_FILE}"

# Allow overriding the singularity command (e.g., SINGULARITY_CMD="sudo apptainer")
SINGULARITY_CMD="${SINGULARITY_CMD:-sudo singularity}"

container_list=(
  "captureenv"
  "prepref"
  "indexgp60db"
  "fqchecks"
  "fastp"
  "fastqc"
  "multiqc"
  "trimgalore"
  "summarise"
  "snippy"
  "moi"
  "moltyping"
  "wgsnv"
  "reporting"
  "bloomine"
)

success=()
failed=()

for item in "${container_list[@]}"; do
  def_file="singularity.${item}"
  sif_file="${item}.sif"

  if [[ ! -f "${def_file}" ]]; then
    echo "[SKIP] ${item}: definition file ${def_file} not found"
    failed+=("${item} (missing def)")
    continue
  fi

  echo "[BUILD] ${item} -> ${sif_file}"
  {
    echo
    echo "[BUILD] ${item} -> ${sif_file}"
  } >> "${LOG_FILE}"

  if ${SINGULARITY_CMD} build "${sif_file}" "${def_file}" >> "${LOG_FILE}" 2>&1; then
    echo "[OK] ${item}"
    success+=("${item}")
  else
    echo "[FAIL] ${item}"
    failed+=("${item}")
  fi
  echo ""
done

echo "========== Build summary =========="
echo "Succeeded: ${#success[@]} -> ${success[*]:-none}"
echo "Failed:    ${#failed[@]} -> ${failed[*]:-none}"
{
  echo "========== Build summary =========="
  echo "Succeeded: ${#success[@]} -> ${success[*]:-none}"
  echo "Failed:    ${#failed[@]} -> ${failed[*]:-none}"
} >> "${LOG_FILE}"

if [[ ${#failed[@]} -gt 0 ]]; then
  exit 1
fi
