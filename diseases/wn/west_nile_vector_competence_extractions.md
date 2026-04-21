# West Nile Vector Competence Extractions By PDF

## Bellini et al. - 2014 - A review of the vector management methods to prevent and control outbreaks of West Nile virus infect.pdf

```csv
West Nile virus,Culex pipiens s.l.,mosquito,competent,narrative_review,yes,yes,Europe,Bellini et al. 2014,major vector role in Europe
West Nile virus,Culex modestus,mosquito,unclear,narrative_review,,,Europe,Bellini et al. 2014,regional role; field-collected species
```

Extraction notes

- Broad Europe-focused management review; partial support only.
- `Culex pipiens s.l.` is the strongest row because the review explicitly says the major vector role in outbreaks seems to be covered by this taxon.
- `Culex modestus` is retained conservatively because the review only gives it a regional role.

## Engler et al. - 2013 - European Surveillance for West Nile Virus in Mosquito Populations.pdf

```csv
West Nile virus,Culex pipiens,mosquito,unclear,narrative_review,,,Europe,Engler et al. 2013,important vector species; WNV-positive pools
West Nile virus,Culex modestus,mosquito,unclear,narrative_review,,,Europe,Engler et al. 2013,important vector species; WNV-positive pools
West Nile virus,Culex perexiguus,mosquito,unclear,narrative_review,,,Spain,Engler et al. 2013,WNV-positive pools
West Nile virus,Ochlerotatus caspius,mosquito,unclear,narrative_review,,,Italy,Engler et al. 2013,WNV-positive pools
```

Extraction notes

- Broad surveillance review; these rows are surveillance-based rather than direct competence experiments.
- I kept the taxa conservative as `unclear` because the source primarily reports positive mosquito pools and vector importance language.
- `Culex pipiens` and `Culex modestus` are the most informative rows here because the review calls them important vector species.

## Ferraguti - 2024 - Mosquito species identity matters unraveling the complex interplay in vector-borne diseases.pdf

```csv
West Nile virus,Aedes caspius,mosquito,unclear,narrative_review,,,Europe,Ferraguti 2024,natural infection; competence debated
West Nile virus,Culex modestus,mosquito,competent,narrative_review,yes,,Europe,Ferraguti 2024,primary WNV vector in Europe
West Nile virus,Culex perexiguus,mosquito,competent,narrative_review,yes,,Europe,Ferraguti 2024,primary WNV vector in Europe
West Nile virus,Culex pipiens,mosquito,competent,narrative_review,yes,,Europe,Ferraguti 2024,primary WNV vector in Europe
```

Extraction notes

- Broad narrative review focused on vector-pathogen specificity.
- `Aedes caspius` is the only clearly uncertain row because the review says its competence remains debated after natural infection detection.
- The three `Culex` rows are treated as competent because the review explicitly calls them primary WNV vectors in Europe.

## Goddard et al. - 2002 - Vector Competence of California Mosquitoes for West Nile virus.pdf

```csv
West Nile virus,Culex tarsalis,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,most efficient laboratory vector
West Nile virus,Cx. p. pipiens,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,highly efficient laboratory vector
West Nile virus,Cx. p. quinquefasciatus,mosquito,mixed,lab_experiment,mixed,,California USA,Goddard et al. 2002,geographic variation; Bakersfield higher than Coachella Valley and Orange County
West Nile virus,Cx. stigmatosoma,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,efficient laboratory vector
West Nile virus,Cx. erythrothorax,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,bridge potential and moderate transmission
West Nile virus,Ochlerotatus dorsalis,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,low to moderate efficiency
West Nile virus,Oc. melanimon,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,low to moderate efficiency
West Nile virus,Oc. sierrensis,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,poor vector; transmitted at low levels
West Nile virus,Aedes vexans,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,secondary role possible
West Nile virus,Culiseta inornata,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,minor role; winter mosquito
```

Extraction notes

- Narrow laboratory competence study of 10 California mosquito species.
- The source is strong for competence because all 10 species became infected and transmitted WNV at some level.
- `Cx. p. quinquefasciatus` is the only mixed row because the paper explicitly shows strong geographic variation in competence.

## Gray and Webb - 2014 - A review of the epidemiological and clinical aspects of West Nile virus.pdf

```csv
West Nile virus,Culex spp.,mosquito,competent,narrative_review,yes,,global,Gray and Webb 2014,primary vectors
West Nile virus,Aedes spp.,mosquito,unclear,narrative_review,,,global,Gray and Webb 2014,lab only; ecological barriers likely limit outbreak importance
West Nile virus,Culex tarsalis,mosquito,competent,narrative_review,yes,,North America,Gray and Webb 2014,major vector
West Nile virus,Culex modestus,mosquito,competent,narrative_review,yes,,Europe,Gray and Webb 2014,locally important vector
```

Extraction notes

- Broad epidemiological review; useful for role language but not a complete competence list.
- I kept the genus-level `Culex spp.` row because the review explicitly frames Culex mosquitoes as the primary vectors.
- `Aedes spp.` remains uncertain because the review only suggests possible efficiency and immediately qualifies the field relevance.

## Hernández-Triana et al. - 2014 - Emergence of West Nile Virus Lineage 2 in Europe A Review on the Introduction and Spread of a Mosqu.pdf

```csv
West Nile virus,Cx. pipiens s.l.,mosquito,unclear,narrative_review,,,Italy,Hernández-Triana et al. 2014,WNV detected in mosquito pools; bridge-vector risk
West Nile virus,Aedes albopictus,mosquito,competent,narrative_review,yes,,Europe,Hernández-Triana et al. 2014,laboratory studies show competent vector
```

Extraction notes

- Broad lineage-2 review; partial support only.
- `Cx. pipiens s.l.` is kept as an uncertainty row because the review mentions WNV-positive pools but does not itself demonstrate competence.
- `Aedes albopictus` is the strongest row because the review explicitly states laboratory infectious studies showed competence for WNV.

## Jansen et al. - 2013 - The Role of Australian Mosquito Species in the Transmission of Endemic and Exotic West Nile Virus St.pdf

```csv
West Nile virus,Cx. annulirostris,mosquito,competent,narrative_review,yes,,Australia,Jansen et al. 2013,primary vector; most competent for exotic WNVNY-99
West Nile virus,Cx. quinquefasciatus,mosquito,competent,narrative_review,yes,,Australia,Jansen et al. 2013,competent vector; primary enzootic vector candidate
West Nile virus,Cx. gelidus,mosquito,competent,narrative_review,yes,,northern Australia,Jansen et al. 2013,highly efficient laboratory vector
West Nile virus,Cx. sitiens,mosquito,unclear,narrative_review,yes,,coastal Australia,Jansen et al. 2013,poor laboratory vector; regional candidate
West Nile virus,Cx. squamosus,mosquito,unclear,narrative_review,,,northern Queensland,Jansen et al. 2013,WNVKUN isolate; regional candidate
West Nile virus,Cx. australicus,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate; potential vector
West Nile virus,Aedes notoscriptus,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes vigilax,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes alternans,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate; secondary vector
West Nile virus,Aedes normanensis,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate; secondary vector
West Nile virus,Anopheles amictus,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate
```

Extraction notes

- Broad Australia-focused review; partial support only.
- The three `Culex` rows at the top are the strongest because the review explicitly calls them competent or highly efficient vectors.
- The remaining rows are conservative candidate or field-isolate rows; several are tied to WNVKUN rather than exotic WNVNY-99, so I kept them `unclear`.

## Martinet et al. - 2019 - Mosquitoes of North-Western Europe as Potential Vectors of Arboviruses A Review.pdf

```csv
West Nile virus,Aedes caspius,mosquito,not_competent,systematic_review,no,,France,Martinet et al. 2019,susceptible to infection but not able to transmit
West Nile virus,Aedes detritus,mosquito,competent,systematic_review,yes,,United Kingdom,Martinet et al. 2019,competence demonstrated
West Nile virus,Aedes geniculatus,mosquito,competent,systematic_review,yes,,Europe,Martinet et al. 2019,tested competent for WNV
West Nile virus,Aedes japonicus japonicus,mosquito,not_competent,systematic_review,no,,Germany,Martinet et al. 2019,could not be infected nor transmit WNV
West Nile virus,Anopheles plumbeus,mosquito,competent,systematic_review,yes,,Europe,Martinet et al. 2019,tested competent; role uncertain
West Nile virus,Culex modestus,mosquito,competent,systematic_review,yes,,France,Martinet et al. 2019,vector incriminated; high transmission rates
West Nile virus,Culex pipiens biotype pipiens,mosquito,competent,systematic_review,yes,,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex pipiens biotype molestus,mosquito,competent,systematic_review,yes,,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex torrentium,mosquito,competent,systematic_review,yes,,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex pipiens s.l.,mosquito,mixed,systematic_review,mixed,,Switzerland,Martinet et al. 2019,susceptible to infection but not competent for WNV lineage 1 FIN Italy
```

Extraction notes

- Systematic review with explicit positive and negative experimental evidence.
- `Aedes caspius`, `Aedes japonicus japonicus`, and `Culex pipiens s.l.` are the key negative or mixed rows.
- I kept `Culex pipiens s.l.` as mixed because the review contains both positive competence evidence for other Culex studies and a negative Swiss population result for WNV lineage 1 FIN Italy.

## Martinet et al. - 2023 - Assessing vector competence of mosquitoes from northeastern France to West Nile virus and Usutu viru.pdf

```csv
West Nile virus,Cx. pipiens,mosquito,mixed,lab_experiment,mixed,,northeastern France,Martinet et al. 2023,some populations transmitted WNV; Maine and Verzy did not
West Nile virus,Ae. rusticus,mosquito,competent,lab_experiment,yes,,northeastern France,Martinet et al. 2023,new putative vector; transmitted from 7 dpi
West Nile virus,Ae. albopictus,mosquito,competent,lab_experiment,yes,,northeastern France,Martinet et al. 2023,transmits WNV and to a lesser extent USUV
West Nile virus,An. plumbeus,mosquito,not_competent,lab_experiment,no,,northeastern France,Martinet et al. 2023,infected but did not transmit WNV
```

Extraction notes

- Narrow laboratory study from northeastern France.
- `Cx. pipiens` is the only mixed row because transmission depended on population: Machault transmitted, Maine and Verzy did not.
- `An. plumbeus` is the clear negative row because the paper explicitly says it became infected but did not transmit WNV.

## Vogels et al. - 2017 - Vector competence of European mosquitoes for West Nile virus.pdf

```csv
West Nile virus,Aedes albopictus,mosquito,competent,narrative_review,yes,,Spain; Italy,Vogels et al. 2017,competent; field relevance likely low in Europe
West Nile virus,Aedes caspius,mosquito,not_competent,narrative_review,no,,France,Vogels et al. 2017,not competent
West Nile virus,Aedes detritus,mosquito,competent,narrative_review,yes,,United Kingdom,Vogels et al. 2017,competent
West Nile virus,Aedes japonicus japonicus,mosquito,not_competent,narrative_review,no,,Germany,Vogels et al. 2017,not competent
West Nile virus,Culex modestus,mosquito,competent,narrative_review,yes,,France,Vogels et al. 2017,efficient vector
West Nile virus,Culex pipiens s.l.,mosquito,competent,narrative_review,yes,,Europe,Vogels et al. 2017,most important vector; transmission 0 to 60 percent
West Nile virus,Culex pipiens pipiens,mosquito,mixed,narrative_review,mixed,,The Netherlands; Italy,Vogels et al. 2017,temperature dependent; transmission varies
West Nile virus,Culex pipiens molestus,mosquito,mixed,narrative_review,mixed,,The Netherlands,Vogels et al. 2017,temperature dependent; transmission varies
West Nile virus,Culex pipiens Hybrid (pipiens x molestus),mosquito,mixed,narrative_review,mixed,,The Netherlands,Vogels et al. 2017,temperature dependent; transmission varies
West Nile virus,Culex torrentium,mosquito,unclear,narrative_review,,,Germany,Vogels et al. 2017,infection and dissemination only; no transmission rate reported
```

Extraction notes

- Broad narrative review of European WNV competence studies.
- The negative rows are `Aedes caspius` and `Aedes japonicus japonicus`; the mixed rows are the `Culex pipiens` biotypes and hybrid because competence varies with temperature.
- `Culex torrentium` stays `unclear` because the review reports infection and dissemination only in the summarized studies and does not give a transmission conclusion for the table row.

## ALL TOGETHER

```csv
West Nile virus,Culex pipiens s.l.,mosquito,competent,narrative_review,yes,yes,Europe,Bellini et al. 2014,major vector role in Europe
West Nile virus,Culex modestus,mosquito,unclear,narrative_review,,,Europe,Bellini et al. 2014,regional role; field-collected species
West Nile virus,Culex pipiens,mosquito,unclear,narrative_review,,,Europe,Engler et al. 2013,important vector species; WNV-positive pools
West Nile virus,Culex modestus,mosquito,unclear,narrative_review,,,Europe,Engler et al. 2013,important vector species; WNV-positive pools
West Nile virus,Culex perexiguus,mosquito,unclear,narrative_review,,,Spain,Engler et al. 2013,WNV-positive pools
West Nile virus,Ochlerotatus caspius,mosquito,unclear,narrative_review,,,Italy,Engler et al. 2013,WNV-positive pools
West Nile virus,Aedes caspius,mosquito,unclear,narrative_review,,,Europe,Ferraguti 2024,natural infection; competence debated
West Nile virus,Culex modestus,mosquito,competent,narrative_review,yes,,Europe,Ferraguti 2024,primary WNV vector in Europe
West Nile virus,Culex perexiguus,mosquito,competent,narrative_review,yes,,Europe,Ferraguti 2024,primary WNV vector in Europe
West Nile virus,Culex pipiens,mosquito,competent,narrative_review,yes,,Europe,Ferraguti 2024,primary WNV vector in Europe
West Nile virus,Culex tarsalis,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,most efficient laboratory vector
West Nile virus,Cx. p. pipiens,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,highly efficient laboratory vector
West Nile virus,Cx. p. quinquefasciatus,mosquito,mixed,lab_experiment,mixed,,California USA,Goddard et al. 2002,geographic variation; Bakersfield higher than Coachella Valley and Orange County
West Nile virus,Cx. stigmatosoma,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,efficient laboratory vector
West Nile virus,Cx. erythrothorax,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,bridge potential and moderate transmission
West Nile virus,Ochlerotatus dorsalis,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,low to moderate efficiency
West Nile virus,Oc. melanimon,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,low to moderate efficiency
West Nile virus,Oc. sierrensis,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,poor vector; transmitted at low levels
West Nile virus,Aedes vexans,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,secondary role possible
West Nile virus,Culiseta inornata,mosquito,competent,lab_experiment,yes,,California USA,Goddard et al. 2002,minor role; winter mosquito
West Nile virus,Culex spp.,mosquito,competent,narrative_review,yes,,global,Gray and Webb 2014,primary vectors
West Nile virus,Aedes spp.,mosquito,unclear,narrative_review,,,global,Gray and Webb 2014,lab only; ecological barriers likely limit outbreak importance
West Nile virus,Culex tarsalis,mosquito,competent,narrative_review,yes,,North America,Gray and Webb 2014,major vector
West Nile virus,Culex modestus,mosquito,competent,narrative_review,yes,,Europe,Gray and Webb 2014,locally important vector
West Nile virus,Cx. pipiens s.l.,mosquito,unclear,narrative_review,,,Italy,Hernández-Triana et al. 2014,WNV detected in mosquito pools; bridge-vector risk
West Nile virus,Aedes albopictus,mosquito,competent,narrative_review,yes,,Europe,Hernández-Triana et al. 2014,laboratory studies show competent vector
West Nile virus,Cx. annulirostris,mosquito,competent,narrative_review,yes,,Australia,Jansen et al. 2013,primary vector; most competent for exotic WNVNY-99
West Nile virus,Cx. quinquefasciatus,mosquito,competent,narrative_review,yes,,Australia,Jansen et al. 2013,competent vector; primary enzootic vector candidate
West Nile virus,Cx. gelidus,mosquito,competent,narrative_review,yes,,northern Australia,Jansen et al. 2013,highly efficient laboratory vector
West Nile virus,Cx. sitiens,mosquito,unclear,narrative_review,yes,,coastal Australia,Jansen et al. 2013,poor laboratory vector; regional candidate
West Nile virus,Cx. squamosus,mosquito,unclear,narrative_review,,,northern Queensland,Jansen et al. 2013,WNVKUN isolate; regional candidate
West Nile virus,Cx. australicus,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate; potential vector
West Nile virus,Aedes notoscriptus,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes vigilax,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes alternans,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate; secondary vector
West Nile virus,Aedes normanensis,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate; secondary vector
West Nile virus,Anopheles amictus,mosquito,unclear,narrative_review,,,Australia,Jansen et al. 2013,WNVKUN isolate
West Nile virus,Aedes caspius,mosquito,not_competent,systematic_review,no,,France,Martinet et al. 2019,susceptible to infection but not able to transmit
West Nile virus,Aedes detritus,mosquito,competent,systematic_review,yes,,United Kingdom,Martinet et al. 2019,competence demonstrated
West Nile virus,Aedes geniculatus,mosquito,competent,systematic_review,yes,,Europe,Martinet et al. 2019,tested competent for WNV
West Nile virus,Aedes japonicus japonicus,mosquito,not_competent,systematic_review,no,,Germany,Martinet et al. 2019,could not be infected nor transmit WNV
West Nile virus,Anopheles plumbeus,mosquito,competent,systematic_review,yes,,Europe,Martinet et al. 2019,tested competent; role uncertain
West Nile virus,Culex modestus,mosquito,competent,systematic_review,yes,,France,Martinet et al. 2019,vector incriminated; high transmission rates
West Nile virus,Culex pipiens biotype pipiens,mosquito,competent,systematic_review,yes,,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex pipiens biotype molestus,mosquito,competent,systematic_review,yes,,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex torrentium,mosquito,competent,systematic_review,yes,,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex pipiens s.l.,mosquito,mixed,systematic_review,mixed,,Switzerland,Martinet et al. 2019,susceptible to infection but not competent for WNV lineage 1 FIN Italy
West Nile virus,Cx. pipiens,mosquito,mixed,lab_experiment,mixed,,northeastern France,Martinet et al. 2023,some populations transmitted WNV; Maine and Verzy did not
West Nile virus,Ae. rusticus,mosquito,competent,lab_experiment,yes,,northeastern France,Martinet et al. 2023,new putative vector; transmitted from 7 dpi
West Nile virus,Ae. albopictus,mosquito,competent,lab_experiment,yes,,northeastern France,Martinet et al. 2023,transmits WNV and to a lesser extent USUV
West Nile virus,An. plumbeus,mosquito,not_competent,lab_experiment,no,,northeastern France,Martinet et al. 2023,infected but did not transmit WNV
West Nile virus,Aedes albopictus,mosquito,competent,narrative_review,yes,,Spain; Italy,Vogels et al. 2017,competent; field relevance likely low in Europe
West Nile virus,Aedes caspius,mosquito,not_competent,narrative_review,no,,France,Vogels et al. 2017,not competent
West Nile virus,Aedes detritus,mosquito,competent,narrative_review,yes,,United Kingdom,Vogels et al. 2017,competent
West Nile virus,Aedes japonicus japonicus,mosquito,not_competent,narrative_review,no,,Germany,Vogels et al. 2017,not competent
West Nile virus,Culex modestus,mosquito,competent,narrative_review,yes,,France,Vogels et al. 2017,efficient vector
West Nile virus,Culex pipiens s.l.,mosquito,competent,narrative_review,yes,,Europe,Vogels et al. 2017,most important vector; transmission 0 to 60 percent
West Nile virus,Culex pipiens pipiens,mosquito,mixed,narrative_review,mixed,,The Netherlands; Italy,Vogels et al. 2017,temperature dependent; transmission varies
West Nile virus,Culex pipiens molestus,mosquito,mixed,narrative_review,mixed,,The Netherlands,Vogels et al. 2017,temperature dependent; transmission varies
West Nile virus,Culex pipiens Hybrid (pipiens x molestus),mosquito,mixed,narrative_review,mixed,,The Netherlands,Vogels et al. 2017,temperature dependent; transmission varies
West Nile virus,Culex torrentium,mosquito,unclear,narrative_review,,,Germany,Vogels et al. 2017,infection and dissemination only; no transmission rate reported
```

