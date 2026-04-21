# GenBank Metadata Source Rules

This note defines when a WHO pathogen should be handled through `nuccore` versus `biosample` for country-level metadata extraction.

## Initial Routing

Use `nuccore` when:

- the pathogen is a virus with `virion_tax_ids`
- the pathogen has clear viral taxonomy support, such as an ICTV/MSL viral name
- the pathogen family is obviously viral
- the pathogen is a virus-like target where accession source qualifiers are usually informative

Why:

- viral accessions often represent one genome or fragment per sampled isolate
- viral `nuccore` records more often carry `country`, `geo_loc_name`, `collection_date`, `host`, and `isolate` directly on the accession

Use `biosample` when:

- the pathogen is a CLOVER-linked bacterial pathogen with `clover_tax_ids`
- the pathogen is non-viral or taxonomically ambiguous
- the biology suggests accession records will be dominated by assemblies, plasmids, contigs, or repeat submissions rather than one-sample-one-record style metadata

Why:

- bacteria often generate many sequence records per biological sample
- country metadata are often better represented at the sample level than on every nucleotide accession

## Repo Rule

In this repo, the default initial routing is:

- `virion_tax_ids` or strong viral signal -> `nuccore`
- `clover_tax_ids` -> `biosample`
- everything else non-viral or ambiguous -> `biosample`

These decisions are now written into the GenBank manifests as:

- `preferred_metadata_source`
- `metadata_source_reason`

## When To Escalate A `nuccore` Pathogen To `biosample`

Even if a pathogen starts in `nuccore`, switch to `biosample` review when:

- `records_found >= 100000` and `countries_observed <= 10`
- or a very large search plateaus early and still yields weak geographic coverage

This is meant to catch cases like large bacterial taxa where the accession universe is huge but country yield stays poor.

The `5_7` second-pass table now writes:

- `followup_metadata_source`
- `followup_metadata_reason`

So we can distinguish:

- `nuccore` second pass: worth fetching more accessions
- `biosample_review`: accession metadata look too redundant or too sparse geographically

## Practical Rule Of Thumb

- viruses: start with `nuccore`
- bacteria: start with `biosample`
- if `nuccore` gives huge accession counts but weak country growth, stop treating it as the primary geography source and move to `biosample`
