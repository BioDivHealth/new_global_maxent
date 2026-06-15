#!/bin/zsh
# Deprecated exploratory/manual workflow. The active broad-taxa entrypoint is
# 05_build_broad_taxa_support.R, which calls the R metadata stage only when
# --refresh-ncbi-metadata is explicitly requested. This shell version is kept
# for provenance around older resolution/TSV/slim review outputs and should not
# be treated as part of the current rebuild contract.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

WHO_DISEASES_OLD_DIR="${REPO_ROOT}/pathogen_association_data/WHO/who_diseases"
WHO_DISEASES_STAGED_BROAD_TAXA_DIR="${REPO_ROOT}/pathogen_association_data/staged/who_diseases/broad_taxa"
WHO_DISEASES_MANUAL_BROAD_TAXA_DIR="${REPO_ROOT}/pathogen_association_data/manual/who_diseases/broad_taxa"

prefer_existing_path() {
  local primary="$1"
  local fallback="$2"
  if [[ -e "${primary}" || ! -e "${fallback}" ]]; then
    printf "%s" "${primary}"
  else
    printf "%s" "${fallback}"
  fi
}

INPUT_CSV="$(prefer_existing_path "${WHO_DISEASES_STAGED_BROAD_TAXA_DIR}/who_broad_taxa_candidate_strains.csv" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_strains.csv")"
RAW_JSONL="$(prefer_existing_path "${WHO_DISEASES_STAGED_BROAD_TAXA_DIR}/who_broad_taxa_candidate_strains_ncbi_raw.jsonl" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_strains_ncbi_raw.jsonl")"
RESOLUTION_CSV="$(prefer_existing_path "${WHO_DISEASES_STAGED_BROAD_TAXA_DIR}/who_broad_taxa_candidate_strains_ncbi_resolution.csv" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_strains_ncbi_resolution.csv")"
METADATA_TSV="$(prefer_existing_path "${WHO_DISEASES_STAGED_BROAD_TAXA_DIR}/who_broad_taxa_candidate_strains_ncbi_metadata.tsv" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_strains_ncbi_metadata.tsv")"
METADATA_CSV="$(prefer_existing_path "${WHO_DISEASES_STAGED_BROAD_TAXA_DIR}/who_broad_taxa_candidate_strains_ncbi_metadata.csv" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_strains_ncbi_metadata.csv")"
ENRICHED_CSV="$(prefer_existing_path "${WHO_DISEASES_STAGED_BROAD_TAXA_DIR}/who_broad_taxa_candidate_strains_ncbi_enriched.csv" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_strains_ncbi_enriched.csv")"
ENRICHED_SLIM_CSV="$(prefer_existing_path "${WHO_DISEASES_STAGED_BROAD_TAXA_DIR}/who_broad_taxa_candidate_strains_ncbi_enriched_slim.csv" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_strains_ncbi_enriched_slim.csv")"
HOST_OVERRIDE_CSV="$(prefer_existing_path "${WHO_DISEASES_MANUAL_BROAD_TAXA_DIR}/who_broad_taxa_candidate_host_overrides.csv" "${WHO_DISEASES_OLD_DIR}/who_broad_taxa_candidate_host_overrides.csv")"

DATASETS_BIN="${REPO_ROOT}/ncbi/datasets"
DATAFORMAT_BIN="${REPO_ROOT}/ncbi/dataformat"

if [[ ! -f "${INPUT_CSV}" ]]; then
  echo "Missing input CSV: ${INPUT_CSV}" >&2
  exit 1
fi

if [[ ! -x "${DATASETS_BIN}" ]]; then
  echo "Missing executable datasets binary: ${DATASETS_BIN}" >&2
  exit 1
fi

if [[ ! -x "${DATAFORMAT_BIN}" ]]; then
  echo "Missing executable dataformat binary: ${DATAFORMAT_BIN}" >&2
  exit 1
fi

TMP_ACCESSIONS="$(mktemp)"
trap 'rm -f "${TMP_ACCESSIONS}"' EXIT

Rscript -e '
input_csv <- commandArgs(trailingOnly = TRUE)[1]
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
acc <- unique(df$accession)
acc <- sub("[.][0-9]+$", "", acc)
acc <- acc[!is.na(acc) & nzchar(acc)]
writeLines(acc, con = stdout())
' "${INPUT_CSV}" > "${TMP_ACCESSIONS}"

: > "${RAW_JSONL}"
printf "accession_base,resolved_accession,ncbi_lookup_status,output_source\n" > "${RESOLUTION_CSV}"

SUCCESS_COUNT=0
FAIL_COUNT=0

while IFS= read -r accession_base; do
  [[ -z "${accession_base}" ]] && continue

  resolved_accession=""
  json_payload=""
  output_source=""

  for suffix in "" ".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8"; do
    accession_try="${accession_base}${suffix}"

    tmp_stdout="$(mktemp)"
    tmp_stderr="$(mktemp)"

    "${DATASETS_BIN}" summary virus genome accession "${accession_try}" --as-json-lines > "${tmp_stdout}" 2> "${tmp_stderr}" < /dev/null || true

    output_stdout="$(cat "${tmp_stdout}")"
    output_stderr="$(cat "${tmp_stderr}")"

    rm -f "${tmp_stdout}" "${tmp_stderr}"

    if [[ "${output_stdout}" == \{* && "${output_stdout}" == *'"accession"'* ]]; then
      resolved_accession="${accession_try}"
      json_payload="${output_stdout}"
      output_source="stdout"
      break
    fi

    if [[ "${output_stderr}" == \{* && "${output_stderr}" == *'"accession"'* ]]; then
      resolved_accession="${accession_try}"
      json_payload="${output_stderr}"
      output_source="stderr"
      break
    fi
  done

  if [[ -n "${resolved_accession}" ]]; then
    printf "%s\n" "${json_payload}" >> "${RAW_JSONL}"
    printf "%s,%s,ok,%s\n" "${accession_base}" "${resolved_accession}" "${output_source}" >> "${RESOLUTION_CSV}"
    SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
  else
    printf "%s,,not_found,\n" "${accession_base}" >> "${RESOLUTION_CSV}"
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
done < "${TMP_ACCESSIONS}"

if [[ -s "${RAW_JSONL}" ]]; then
  "${DATAFORMAT_BIN}" tsv virus-genome \
    --inputfile "${RAW_JSONL}" \
    --fields accession,completeness,geo-location,geo-region,host-common-name,host-name,host-tax-id,is-annotated,is-lab-host,is-vaccine-strain,isolate-collection-date,length,protein-count,release-date,sourcedb,submitter-affiliation,submitter-country,submitter-names,update-date,virus-common-name,virus-infraspecific-isolate,virus-name,virus-tax-id \
    > "${METADATA_TSV}"
else
  printf "accession\tcompleteness\tgeo-location\tgeo-region\thost-common-name\thost-name\thost-tax-id\tis-annotated\tis-lab-host\tis-vaccine-strain\tisolate-collection-date\tlength\tprotein-count\trelease-date\tsourcedb\tsubmitter-affiliation\tsubmitter-country\tsubmitter-names\tupdate-date\tvirus-common-name\tvirus-infraspecific-isolate\tvirus-name\tvirus-tax-id\n" > "${METADATA_TSV}"
fi

Rscript -e '
args <- commandArgs(trailingOnly = TRUE)
candidate_csv <- args[1]
resolution_csv <- args[2]
metadata_tsv <- args[3]
metadata_csv <- args[4]
enriched_csv <- args[5]
enriched_slim_csv <- args[6]
host_override_csv <- args[7]

candidate <- read.csv(candidate_csv, stringsAsFactors = FALSE, check.names = FALSE)
resolution <- read.csv(resolution_csv, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.delim(metadata_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
host_overrides <- if (file.exists(host_override_csv)) {
  read.csv(host_override_csv, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(accession_base = character(), stringsAsFactors = FALSE, check.names = FALSE)
}

if ("Accession" %in% names(metadata) && !"accession" %in% names(metadata)) {
  names(metadata)[names(metadata) == "Accession"] <- "accession"
}

if (!"accession" %in% names(metadata)) {
  metadata <- data.frame(accession = character(), stringsAsFactors = FALSE, check.names = FALSE)
}

metadata_output <- merge(
  resolution,
  metadata,
  by.x = "resolved_accession",
  by.y = "accession",
  all.x = TRUE,
  sort = TRUE
)

candidate$accession_base <- sub("[.][0-9]+$", "", candidate$accession)

enriched_output <- merge(
  candidate,
  metadata_output,
  by.x = "accession_base",
  by.y = "accession_base",
  all.x = TRUE,
  sort = FALSE
)

enriched_output <- merge(
  enriched_output,
  host_overrides,
  by = "accession_base",
  all.x = TRUE,
  sort = FALSE
)

if ("Host Name" %in% names(enriched_output)) {
  enriched_output[["host_name_datasets"]] <- enriched_output[["Host Name"]]
} else {
  enriched_output[["host_name_datasets"]] <- NA_character_
}

enriched_output[["host_name_record_or_override"]] <- ifelse(
  !is.na(enriched_output[["host_name_record"]]) & nzchar(enriched_output[["host_name_record"]]),
  enriched_output[["host_name_record"]],
  enriched_output[["host_name_datasets"]]
)

enriched_output[["host_name_current_curated"]] <- ifelse(
  !is.na(enriched_output[["host_name_current_override"]]) & nzchar(enriched_output[["host_name_current_override"]]),
  enriched_output[["host_name_current_override"]],
  enriched_output[["host_name_record_or_override"]]
)

slim_columns <- c(
  "accession_base",
  "broad_group",
  "proposed_active_unit",
  "virus_name",
  "isolate",
  "abbrev",
  "accession",
  "resolved_accession",
  "ncbi_lookup_status",
  "output_source",
  "candidate_role",
  "ictv_genus",
  "ictv_subgenus",
  "ictv_species",
  "available_sequence",
  "Geographic Location",
  "Geographic Region",
  "Host Name",
  "Host Taxonomic ID",
  "host_name_datasets",
  "host_name_record",
  "host_name_current_override",
  "host_name_record_or_override",
  "host_name_current_curated",
  "host_override_note",
  "Virus Name",
  "Virus Taxonomic ID"
)

enriched_output_slim <- enriched_output[, intersect(slim_columns, names(enriched_output)), drop = FALSE]

write.csv(metadata_output, metadata_csv, row.names = FALSE, na = "")
write.csv(enriched_output, enriched_csv, row.names = FALSE, na = "")
write.csv(enriched_output_slim, enriched_slim_csv, row.names = FALSE, na = "")

cat("Candidate rows read:", nrow(candidate), "\n")
cat("Distinct accession bases queried:", nrow(resolution), "\n")
cat("Successful NCBI metadata lookups:", sum(resolution$ncbi_lookup_status == "ok"), "\n")
cat("Failed NCBI metadata lookups:", sum(resolution$ncbi_lookup_status != "ok"), "\n")
cat("Wrote resolution table to", resolution_csv, "\n")
cat("Wrote metadata TSV to", metadata_tsv, "\n")
cat("Wrote metadata CSV to", metadata_csv, "\n")
cat("Wrote enriched candidate table to", enriched_csv, "\n")
cat("Wrote slim enriched candidate table to", enriched_slim_csv, "\n")
' "${INPUT_CSV}" "${RESOLUTION_CSV}" "${METADATA_TSV}" "${METADATA_CSV}" "${ENRICHED_CSV}" "${ENRICHED_SLIM_CSV}" "${HOST_OVERRIDE_CSV}"

echo "Wrote raw JSONL to ${RAW_JSONL}"
