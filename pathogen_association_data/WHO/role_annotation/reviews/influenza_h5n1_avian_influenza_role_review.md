# Influenza H5N1 Avian Influenza Role Review

Phase: `Phase N`
Started: `2026-05-08`
Last updated: `2026-05-08`

## Disease Scope And Local Candidate Counts

- Disease name: Influenza (H5N1 avian influenza)
- Source pathogen or analysis unit: Alphainfluenzavirus influenzae (H5N1)
- Host candidate rows: 538
- Vector candidate rows: 0
- Competence-linked vector rows: 0
- Local candidate snapshot: large host-only candidate set dominated by avian hosts with livestock and mammal spillover candidates.

## Current Host Candidate Highlights

- Wild aquatic birds are supported at group level as natural hosts/reservoirs for avian influenza A viruses.
- `Gallus gallus` is present as a domestic poultry candidate and is supported as a poultry outbreak/amplification host rather than natural reservoir.
- `Bos taurus` is present as a local candidate and is supported as a recent susceptible mammalian host in A(H5N1) dairy cattle outbreaks.
- `Homo sapiens` is present as a local candidate and is treated as an incidental zoonotic host because sustained human-to-human transmission is not demonstrated.

## Current Vector Candidate Highlights

- `not_applicable_non_vectored_scope`: no current disease-vector rows are present for this Phase N disease.

## Sources Searched

| Source | Type | URL or local path | Used for rows? | Notes |
|---|---|---|---|---|
| CDC Bird Flu Current Situation in Wild Birds | official public health guidance | https://www.cdc.gov/bird-flu/situation-summary/wildbirds.html | yes | Used for wild aquatic bird natural-host group row. |
| WHO Influenza avian and other zoonotic fact sheet | official public health factsheet | https://www.who.int/news-room/fact-sheets/detail/influenza-%28avian-and-other-zoonotic%29 | yes | Used for poultry epizootic cattle susceptibility and human spillover rows. |
| CDC Avian Influenza Type A | official public health background | https://www.cdc.gov/bird-flu/about/avian-influenza-type-a.html | supporting context | Used to cross-check host groups and human infection framing. |
| Local host role candidates | local candidate table | `pathogen_association_data/WHO/role_annotation/host_role_candidates.csv` | yes | Used to verify local candidate presence and tax_id for `Gallus gallus` `Bos taurus` and `Homo sapiens`. |

## Source-Backed Host Role Findings

| Host or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| wild aquatic birds | `reservoir_host_group` | supports | high | yes | CDC identifies wild aquatic birds including gulls terns shorebirds waterfowl ducks geese and swans as natural hosts. |
| `Gallus gallus` | `amplifying_host` | supports | medium | yes | WHO describes recurrent poultry epizootics and hundreds of millions of poultry infections. |
| `Bos taurus` | `susceptible_host_only` | supports | medium | yes | WHO reports A(H5N1) outbreaks in United States dairy cattle in 2024 and infections among exposed workers. |
| `Homo sapiens` | `incidental_host` | supports | high | no | WHO states zoonotic influenza human infections are sporadic and sustained human-to-human transmission has not been demonstrated. |

## Source-Backed Vector Role Findings

| Vector or group | Role claim | Evidence direction | Confidence | Manual review? | Evidence note |
|---|---|---|---|---|---|
| not applicable | `not_applicable_non_vectored_scope` | not applicable | not applicable | no | Phase N host-only review; no vector rows added. |

## Rows Added To Evidence CSVs

- Host evidence rows: 4
- Vector evidence rows: 0

## Draft Assignments Added

- Host assignments: 4
- Vector assignments: 0

## Deferred Candidates And Why

- Avian species-level assignments: deferred because the strongest reviewed source supports wild aquatic birds at group level.
- Other mammals: deferred unless source-backed evidence supports a role beyond infection or host presence.
- Swine and other livestock: deferred pending H5N1-specific role evidence that distinguishes susceptible host from amplifier or reassortment concern.

## Open Questions For Collaborator Review

- Whether domestic poultry should remain `amplifying_host` or move to a poultry-specific vocabulary value in a future role schema.
- Whether dairy cattle should remain evidence-only or receive a stronger assignment after a dedicated 2024-2026 H5N1 cattle literature review.
