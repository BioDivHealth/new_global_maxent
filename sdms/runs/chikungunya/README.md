# Chikungunya SDM Pilot

This folder tracks the disease-level SDM run queue for the Chikungunya pilot.
Final species-level SDM objects should stay under `sdms/models/{species_name}/`
so they can be reused by other diseases.

## Files

- `sdm_target_manifest.csv` is the current run manifest for Chikungunya species.
  It is derived from the generated pilot package under
  `pathogen_association_data/readiness/disease_modelling_pilot_package/`.

## Current Scope

- Disease: `Chikungunya`
- Manifest rows: one row per host or vector species in the pilot SDM table.
- Existing host SDMs are marked `already_available`.
- Missing host and vector SDMs are marked `not_run`.
- Occurrence-count fields are placeholders until vector and host occurrence
  inputs are assembled.

## Storage Boundary

Use this folder for Chikungunya-specific manifests, logs, and run notes. Do not
store final model objects here. Store final reusable species SDMs in `sdms/models/`
and update the catalog outputs under `sdms/outputs/catalog/`.
