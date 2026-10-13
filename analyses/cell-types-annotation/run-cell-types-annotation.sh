#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")"

CONFIG_FILE="../../project_parameters.Config.yaml"
if [[ ! -f "${CONFIG_FILE}" ]]; then
  echo "Error: configuration file not found: ${CONFIG_FILE}" >&2
  exit 1
fi

# Which annotation pipeline steps to run (must match run-merge-cell-types-annotations-all.R)
method=$(grep '^method_annotation_module:' "${CONFIG_FILE}" | awk '{print $2}' | tr -d '"')

run_singler=false
run_gene_markers=false
run_reference=false

case "${method}" in
  all)
    run_singler=true
    run_gene_markers=true
    run_reference=true
    ;;
  singler)
    run_singler=true
    ;;
  gene_markers)
    run_gene_markers=true
    ;;
  reference)
    run_reference=true
    ;;
  singler\&gene_markers)
    run_singler=true
    run_gene_markers=true
    ;;
  singler\&reference)
    run_singler=true
    run_reference=true
    ;;
  gene_markers\&reference)
    run_gene_markers=true
    run_reference=true
    ;;
  *)
    echo "Error: unsupported method_annotation_module: '${method}'" >&2
    echo "Valid options: all, singler, gene_markers, reference, singler&gene_markers, singler&reference, gene_markers&reference" >&2
    exit 1
    ;;
esac

################################################################################################################
# Run module (steps selected via method_annotation_module in project_parameters.Config.yaml)
if [[ "${run_singler}" == true ]]; then
  Rscript --vanilla run-cell-types-annotation-SingleR.R
fi
if [[ "${run_gene_markers}" == true ]]; then
  Rscript --vanilla run-cell-types-annotation-gene-markers.R
fi
if [[ "${run_reference}" == true ]]; then
  Rscript --vanilla run-cell-types-annotation-reference.R
fi
Rscript --vanilla run-merge-cell-types-annotations-all.R
################################################################################################################
