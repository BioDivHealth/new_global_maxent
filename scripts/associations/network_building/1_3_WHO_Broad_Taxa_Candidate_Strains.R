# ------------------------------------------------------------------------------
# 1_3_WHO_Broad_Taxa_Candidate_Strains.R
# ------------------------------------------------------------------------------
# Purpose: Build a curation inventory of ICTV-supported candidate strains and
#          exemplar viruses for the broad taxa currently under review.
#
# Input  : ICTV-derived member tables for Sarbecovirus, Merbecovirus, and
#          Vesiculovirus, plus the current zoonotic WHO analysis-unit shortlist.
# Output : pathogen_association_data/WHO/who_diseases/
#          who_broad_taxa_candidate_strains.csv
# ------------------------------------------------------------------------------

library(pacman)
p_load(dplyr, here, readr, stringr, tibble)

clean_text <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "NA", "NaN", "No data", "null", "Null")] <- NA_character_
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_replace_all(x, "[\r\n\t]+", " ")
  x <- stringr::str_squish(x)
  x[x == ""] <- NA_character_
  x
}

who_dir <- here("pathogen_association_data", "WHO", "who_diseases")
output_path <- file.path(who_dir, "who_broad_taxa_candidate_strains.csv")
analysis_units_keep_path <- file.path(who_dir, "who_pathogen_analysis_units_keep.csv")

analysis_units_keep <- read_csv(
  analysis_units_keep_path,
  show_col_types = FALSE,
  na = c("", "NA")
) %>%
  mutate(across(where(is.character), clean_text))

candidate_strains <- tibble::tribble(
  ~broad_group, ~ictv_genus, ~ictv_subgenus, ~ictv_species, ~virus_name, ~isolate, ~accession, ~available_sequence, ~abbrev, ~candidate_role, ~proposed_active_unit, ~decision, ~decision_reason, ~ictv_taxon_anchor, ~ictv_taxon_source, ~notes,
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", "Betacoronavirus pandemicum", "severe acute respiratory syndrome coronavirus", "Tor2", "AY274119", "Complete genome", "SARS-CoV", "active_unit_exemplar", "SARS-CoV-1", "keep_active", "ICTV exemplar for SARS-CoV-1", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Formal ICTV exemplar row for the first SARS outbreak virus.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", "Betacoronavirus pandemicum", "SARS coronavirus", "PC4-227", "AY613950", "Complete genome", "SARS-CoV", "supporting_exemplar", "SARS-CoV-1", "keep_supporting_example", "Second ICTV exemplar for SARS-CoV-1", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Use as supporting evidence for the SARS-CoV-1 analytic unit.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", "Betacoronavirus pandemicum", "severe acute respiratory syndrome coronavirus 2", "Wuhan-Hu-1", "MN908947", "Complete genome", "SARS-CoV-2", "active_unit_exemplar", "SARS-CoV-2", "keep_active", "ICTV exemplar for SARS-CoV-2", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Formal ICTV exemplar row for the COVID-19 virus.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", "Betacoronavirus pandemicum", "severe acute respiratory syndrome-related coronavirus", "BtKY72", "KY352407", "Complete genome", "SARSr-CoV", "wildlife_group_exemplar", "SARS-like bat sarbecoviruses", "keep_active", "ICTV exemplar for a bat sarbecovirus candidate group", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Keep as the bat-sarbecovirus group exemplar while collecting additional NCBI names and host data.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat coronavirus RaTG13", "RaTG13", "MN996532", "Complete genome", "RaTG13", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed bat sarbecovirus example with an NCBI accession available for host/location enrichment.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat coronavirus RmYN02", "RmYN02", NA, "Partial NCBI gene accessions only", "RmYN02", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed bat sarbecovirus example. No single clean whole-genome NCBI nucleotide accession pinned in this pass; current NCBI search returns partial gene accessions including MW201981.1 and MW201982.1.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS-like coronavirus WIV1", "WIV1", "KF367457", "Complete genome", "WIV1", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed bat sarbecovirus example with reported ACE2-usage evidence.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS-like coronavirus SHC014", "SHC014", "KC881005", "Complete genome", "SHC014", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed SHC014 example. Uses the same NCBI accession as the RsSHC014 naming variant and should be treated as an alias-level supporting row for now.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS-like coronavirus Rs3367", "Rs3367", "KC881006", "Complete genome", "Rs3367", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed bat sarbecovirus example from Rhinolophus sinicus.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS-like coronavirus RsSHC014", "RsSHC014", "KC881005", "Complete genome", "RsSHC014", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed bat sarbecovirus example from Rhinolophus sinicus. Shares an accession with the SHC014 naming variant.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS-like coronavirus Rs4028", "Rs4028", NA, "Name unresolved in current NCBI pass", "Rs4028", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed bat sarbecovirus example. Exact NCBI accession unresolved in this pass, so left blank rather than inferring from similarly named strains.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat coronavirus isolate PREDICT/PDF-2370/OTBA35RSV", "PDF-2370", "MT726044", "Complete genome", "PDF-2370", "supporting_named_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Named sarbecovirus example highlighted for zoonotic-potential review", "Subgenus Sarbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506129&taxon_name=Sarbecovirus", "Dave-listed bat sarbecovirus example. Restored with the NCBI accession MT726044.1, but keep in mind the record metadata should still be manually reconciled against Dave's shorthand location description.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS coronavirus HKU3-1", "HKU3-1", "DQ022305", "Complete genome", "HKU3", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=694009&lvl=", "Added as a representative HKU3 lineage example from the Sarbecovirus taxonomy tree using the commonly cited complete-genome accession DQ022305.2.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Rhinolophus affinis coronavirus isolate LYRa11", "LYRa11", "KF569996", "Complete genome", "LYRa11", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=694009&lvl=", "Added as a host-linked sarbecovirus example from the NCBI taxonomy tree with a clear Rhinolophus affinis label.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "BtRs-BetaCoV/GX2013", "GX2013", "KJ473815", "Complete genome", "GX2013", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=694009&lvl=", "Added as a 2013 bat sarbecovirus genome from the NCBI taxonomy tree with a direct GenBank accession.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "BtRs-BetaCoV/HuB2013", "HuB2013", "KJ473814", "Complete genome", "HuB2013", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=694009&lvl=", "Added as a 2013 bat sarbecovirus genome from the NCBI taxonomy tree. Accession KJ473814 is widely cited in sarbecovirus comparison papers.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "BtRs-BetaCoV/YN2013", "YN2013", "KJ473816", "Complete genome", "YN2013", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=694009&lvl=", "Added as a 2013 bat sarbecovirus genome from the NCBI taxonomy tree. Accession KJ473816 is widely cited in sarbecovirus comparison papers.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS CoV Rp3/2004", "Rp3/2004", "DQ071615", "Complete genome", "Rp3", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=694009&lvl=", "Added as a classic bat SARS-related coronavirus genome commonly used in comparative sarbecovirus analyses.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "SARSr-Rf-BatCoV YNLF_31C", "YNLF_31C", "KP886808", "GenBank accession provided", "YNLF_31C", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509511&lvl=3", "Added from the Sarbecovirus taxonomy tree and comparative literature. Manual nuccore lookup pinned accession KP886808.1.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "SARSr-Rf-BatCoV YNLF_34C", "YNLF_34C", "KP886809", "GenBank accession provided", "YNLF_34C", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509511&lvl=3", "Added from the Sarbecovirus taxonomy tree and comparative literature. Manual nuccore lookup pinned accession KP886809.1.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "BtRf-BetaCoV/HeB2013", "BtRf-BetaCoV/HeB2013", "KJ473812", "GenBank accession provided", "HeB2013", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509511&lvl=3", "Added from the Sarbecovirus taxonomy tree. Manual nuccore lookup pinned accession KJ473812.1.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "BtRf-BetaCoV/JL2012", "BtRf-BetaCoV/JL2012", "KJ473811", "GenBank accession provided", "JL2012", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509511&lvl=3", "Added from the Sarbecovirus taxonomy tree. Manual nuccore lookup pinned accession KJ473811.1.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "BtRf-BetaCoV/SX2013", "BtRf-BetaCoV/SX2013", "KJ473813", "GenBank accession provided", "SX2013", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509511&lvl=3", "Added from the Sarbecovirus taxonomy tree. Manual nuccore lookup pinned accession KJ473813.1.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS-related coronavirus Rf1/2004", "Rf1/2004", "DQ412042", "GenBank accession provided", "Rf1", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509511&lvl=3", "Added from the Sarbecovirus taxonomy tree. Manual nuccore lookup pinned accession DQ412042.1.",
  "Sarbecovirus", "Betacoronavirus", "Sarbecovirus", NA, "Bat SARS-related coronavirus Rm1/2004", "Rm1/2004", "DQ412043", "GenBank accession provided", "Rm1", "supporting_tree_example", "SARS-like bat sarbecoviruses", "keep_supporting_example", "Wildlife sarbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Sarbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509511&lvl=3", "Added from the Sarbecovirus taxonomy tree. Manual nuccore lookup pinned accession DQ412043.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", "Betacoronavirus cameli", "Middle East respiratory syndrome-related coronavirus", "HCoV-EMC", "JX869059", "Complete genome", "MERS-CoV", "active_unit_exemplar", "MERS-CoV", "keep_active", "ICTV exemplar for MERS-CoV", "Subgenus Merbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506125&taxon_name=Merbecovirus", "Formal ICTV exemplar row for the MERS virus.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", "Betacoronavirus erinacei", "hedgehog coronavirus 1; Erinaceus europaeus coronavirus 1", "2012-174/GER/2012", "KC545383", "Complete genome", "EriCoV1", "supporting_example", "Merbecovirus review set", "keep_supporting_example", "Extra ICTV merbecovirus exemplar, not part of the initial active set", "Subgenus Merbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506125&taxon_name=Merbecovirus", "Useful as a reference candidate but not one of the initial active host-distribution targets.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", "Betacoronavirus pipistrelli", "Pipistrellus bat coronavirus HKU5", "LMH03f", "EF065509", "Complete genome", "Pi-BatCoV_HKU5", "wildlife_group_exemplar", "MERS-like bat merbecoviruses", "keep_active", "ICTV exemplar for a bat merbecovirus candidate group", "Subgenus Merbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506125&taxon_name=Merbecovirus", "Use as a bat-merbecovirus group exemplar while collecting host and zoonotic-potential evidence.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", "Betacoronavirus tylonycteridis", "Tylonycteris bat coronavirus HKU4", "B04f", "EF065505", "Complete genome", "Ty-BatCoV_HKU4", "wildlife_group_exemplar", "MERS-like bat merbecoviruses", "keep_active", "ICTV exemplar for a bat merbecovirus candidate group", "Subgenus Merbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506125&taxon_name=Merbecovirus", "Use as a bat-merbecovirus group exemplar while collecting host and zoonotic-potential evidence.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus isolate PREDICT/PDF-2180", "PDF-2180", "KX574227", "Complete genome", "PDF-2180", "supporting_named_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Named merbecovirus example highlighted for zoonotic-potential review", "Subgenus Merbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506125&taxon_name=Merbecovirus", "Dave-listed bat merbecovirus example. Restored with the GenBank accession KX574227.1; NC_034440.1 is the matching RefSeq companion if needed for later comparison.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "BtVs-BetaCoV/SC2013", "BtVs-BetaCoV/SC2013", "KJ473821", "Complete genome", "BtVs-BetaCoV/SC2013", "supporting_named_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Named merbecovirus example highlighted for zoonotic-potential review", "Subgenus Merbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506125&taxon_name=Merbecovirus", "Dave-listed bat merbecovirus example from Vespertilio superans with an NCBI accession available for enrichment.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "NeoCoV", "NeoCoV", "KC869678", "Complete genome", "NeoCoV", "supporting_named_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Named merbecovirus example highlighted for zoonotic-potential review", "Subgenus Merbecovirus", "https://ictv.global/taxonomy/taxondetails?taxnode_id=202506125&taxon_name=Merbecovirus", "Dave-listed bat merbecovirus example with published receptor-binding concern and an NCBI accession available for enrichment.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "BtPa-BetaCoV/GD2013", "BtPa-BetaCoV/GD2013", "KJ473820", "GenBank accession provided", "GD2013", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as an HKU5-lineage supporting example. Manual nuccore lookup pinned accession KJ473820.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "BtTp-BetaCoV/GX2012", "BtTp-BetaCoV/GX2012", "KJ473822", "GenBank accession provided", "GX2012", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus example promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as an HKU4-lineage supporting example. Manual nuccore lookup pinned accession KJ473822.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU5-1", "HKU5-1", NA, "Accession unresolved in current pass", "HKU5-1", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU5 isolate variant. Exact NCBI nucleotide accession not pinned in this pass.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU5-2", "HKU5-2", "EF065510", "GenBank accession provided", "HKU5-2", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU5 isolate variant. Manual nuccore lookup pinned accession EF065510.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU5-3", "HKU5-3", "EF065511", "GenBank accession provided", "HKU5-3", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU5 isolate variant. Manual nuccore lookup pinned accession EF065511.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU5-5", "HKU5-5", "EF065512", "GenBank accession provided", "HKU5-5", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU5 isolate variant. Manual nuccore lookup pinned accession EF065512.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU4-1", "HKU4-1", NA, "Accession unresolved in current pass", "HKU4-1", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU4 isolate variant. Exact NCBI nucleotide accession not pinned in this pass.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU4-2", "HKU4-2", "EF065506", "GenBank accession provided", "HKU4-2", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU4 isolate variant. Manual nuccore lookup pinned accession EF065506.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU4-3", "HKU4-3", "EF065507", "GenBank accession provided", "HKU4-3", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU4 isolate variant. Manual nuccore lookup pinned accession EF065507.1.",
  "Merbecovirus", "Betacoronavirus", "Merbecovirus", NA, "Bat coronavirus HKU4-4", "HKU4-4", "EF065508", "GenBank accession provided", "HKU4-4", "supporting_tree_example", "MERS-like bat merbecoviruses", "keep_supporting_example", "Wildlife merbecovirus isolate variant promoted from the NCBI taxonomy tree", "Subgenus Merbecovirus", "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?command=show&mode=tree&id=2509494&lvl=3", "Added from the Merbecovirus taxonomy tree as a named HKU4 isolate variant. Manual nuccore lookup pinned accession EF065508.1.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus alagoas", "vesicular stomatitis Alagoas virus", "Indiana 3", "EU373658", "Complete genome", "VSAV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus bogdanovac", "Yug Bogdanovac virus", "YU4-76", "JF911700", "Coding-complete genome", "YBV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus carajas", "Carajás virus", "BeAr411391", "KM205015", "Coding-complete genome", "CARV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus chandipura", "Chandipura virus", "CIN 0451", "GU212856", "Coding-complete genome", "CHPV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus cocal", "Cocal virus", "TRVL40233", "EU373657", "Complete genome", "COCV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus indiana", "vesicular stomatitis Indiana virus", "98COE", "AF473864", "Complete genome", "VSIV", "active_unit_exemplar", "Vesicular stomatitis Indiana virus", "keep_active", "ICTV exemplar for the active Vesiculovirus unit VSIV", NA, NA, "Active unit to carry forward for host/vector/amplifier lookups.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus isfahan", "Isfahan virus", "91026-167", "AJ810084", "Complete genome", "ISFV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus jurona", "Jurona virus", "BeAr40578", "KM204996", "Coding-complete genome", "JURV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus malpais", "Malpais Spring virus", "85-488NM", "KC412247", "Complete genome", "MSPV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus maraba", "Maraba virus", "BeAr 411459", "HQ660076", "Complete genome", "MARAV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus mejal", "Mejal virus", "JAL10", "MW798173", "Coding-complete genome", "MEJV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus morreton", "Morreton virus", "CoAr191048", "KM205007", "Complete genome", "MORV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus newjersey", "vesicular stomatitis New Jersey virus", "NJ1184HDB", "JX121109", "Complete genome", "VSNJV", "active_unit_exemplar", "Vesicular stomatitis New Jersey virus", "keep_active", "ICTV exemplar for the active Vesiculovirus unit VSNJV", NA, NA, "Active unit to carry forward for host/vector/amplifier lookups.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus perinet", "Perinet virus", "DAkAr Mg802", "HM566195", "Complete genome", "PERV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus piry", "Piry virus", "BeAn2423", "KU178986", "Coding-complete genome", "PIRYV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review.",
  "Vesiculovirus", "Vesiculovirus", NA, "Vesiculovirus radi", "Radi virus", "ISS Ph1-166", "KM205024", "Complete genome", "RADV", "supporting_example", "Vesiculovirus review set", "keep_supporting_example", "ICTV vesiculovirus species candidate, not in the initial active vesicular-stomatitis pair", NA, NA, "Useful supporting species for later Vesiculovirus review."
)

candidate_strains <- candidate_strains %>%
  mutate(across(where(is.character), clean_text)) %>%
  mutate(
    source_rationale = decision_reason,
    source_notes = notes
  ) %>%
  select(-notes) %>%
  left_join(
    analysis_units_keep %>%
      transmute(
        analysis_unit,
        keep_reason,
        keep_rationale = rationale,
        keep_notes = notes,
        do_not_split_further_yet
      ),
    by = c("proposed_active_unit" = "analysis_unit")
  ) %>%
  mutate(
    keep_reason = dplyr::coalesce(keep_reason, decision_reason),
    keep_rationale = dplyr::coalesce(keep_rationale, source_rationale),
    keep_notes = dplyr::coalesce(keep_notes, source_notes),
    support_status = dplyr::case_when(
      decision == "keep_active" ~ "active_unit",
      decision == "keep_supporting_example" ~ "supporting_example",
      TRUE ~ "review_only"
    ),
    ncbi_followup_hint = dplyr::case_when(
      decision == "keep_active" ~ paste0(virus_name, " | ", isolate, " | ", accession),
      TRUE ~ paste0(virus_name, " | ", isolate)
    )
  ) %>%
  select(
    broad_group,
    ictv_genus,
    ictv_subgenus,
    ictv_species,
    virus_name,
    isolate,
    accession,
    available_sequence,
    abbrev,
    candidate_role,
    proposed_active_unit,
    decision,
    decision_reason,
    support_status,
    ncbi_followup_hint,
    ictv_taxon_anchor,
    ictv_taxon_source,
    keep_reason,
    keep_rationale,
    keep_notes,
    source_rationale,
    source_notes
  ) %>%
  arrange(broad_group, desc(decision == "keep_active"), ictv_species, virus_name, isolate)

write_csv(candidate_strains, output_path, na = "")

cat("Candidate strain rows written:", nrow(candidate_strains), "\n")
cat("Active-unit rows:", sum(candidate_strains$decision == "keep_active"), "\n")
cat("Supporting-example rows:", sum(candidate_strains$decision == "keep_supporting_example"), "\n")
cat("Review-only rows:", sum(candidate_strains$decision == "review_only"), "\n")
cat("Wrote candidate strain table to", output_path, "\n")
