# Host-Vector Evidence

Integrated host-vector evidence surface generated from the staged VectorMap and
MapVEu source-family outputs.

## Files

- `vector_host_links_analysis_ready.csv`: Combined record-level host-vector
  evidence table built from the VectorMap and MapVEu analysis-ready outputs,
  with source provenance preserved.
- `vector_host_links_analysis_summary.csv`: Deduplicated host-vector summary
  table aggregated from the combined record-level evidence.
- `vector_host_links_join_ready.csv`: Canonical host-vector join table used by
  downstream disease-host-vector and pathogen-host-vector joins.
- `vector_host_links_join_blocked.csv`: Non-joinable or cautionary records kept
  for QA and review rather than dropped silently.

## Usage

Active scripts should use `vector_host_outputs_dir` from
`scripts/associations/working_inputs.R` rather than hard-coding this path.
