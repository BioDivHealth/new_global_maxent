# VecTraits exploratory probes

These scripts are optional API-backed probes for VectorByte VecTraits data. They
are not part of the routine vector-screening pipeline and are not used as
disease-vector inclusion evidence.

Outputs are written under:

- `pathogen_association_data/staged/vector_screening/vectraits/`

That output folder is ignored because these files are local exploratory
snapshots of an external API surface.

Scripts:

- `probe_vectraits_access.R`: small access/schema probe for selected searches.
- `summarize_vectraits_traits.R`: trait-oriented summary across selected
  species and keywords.
- `build_vectraits_species_tables.R`: row-level exploratory trait tables for
  selected vector species.

Run these from the repository root. Set `VECTRAITS_USE_QA=true` to query the
QA VecTraits API host where supported by the script.
