# West Nile Vector Extractions By PDF

## Bellini et al. - 2014 - A review of the vector management methods to prevent and control outbreaks of West Nile virus infect.pdf

```csv
West Nile virus,Culex pipiens s.l.,mosquito,confirmed,review,Europe,Bellini et al. 2014,major vector role
West Nile virus,Culex modestus,mosquito,probable,review,Europe,Bellini et al. 2014,regional role
West Nile virus,Culex univittatus,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Coquillettidia richiardii,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes cantans,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes caspius,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes excrucians,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes vexans,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Anopheles maculipennis s.s.,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Anopheles atroparvus,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
```

Extraction notes

- Broad Europe-focused review; partial support only, not a complete species list.
- `Culex pipiens s.l.` is the strongest row because the paper explicitly says it covers the major vector role in European outbreaks.
- The remaining species were retained only as conservative candidate rows because the source frames them as field-collected WNV-positive mosquitoes rather than established vectors.

## Engler et al. - 2013 - European Surveillance for West Nile Virus in Mosquito Populations.pdf

```csv
West Nile virus,Culex pipiens,mosquito,probable,review,Europe,Engler et al. 2013,important vector species
West Nile virus,Culex modestus,mosquito,probable,review,Europe,Engler et al. 2013,important vector species
West Nile virus,Culex perexiguus,mosquito,candidate,field,Spain,Engler et al. 2013,WNV-positive pools
West Nile virus,Ochlerotatus caspius,mosquito,weak,field,Italy,Engler et al. 2013,WNV-positive pools
```

Extraction notes

- Broad surveillance review; partial support only, not a complete species list.
- I kept `Culex pipiens` and `Culex modestus` at `probable` because the paper explicitly calls them important vector species.
- `Culex perexiguus` and `Ochlerotatus caspius` were kept conservative because this source supports them mainly through positive-pool surveillance.

## Ferraguti - 2024 - Mosquito species identity matters unraveling the complex interplay in vector-borne diseases.pdf

```csv
West Nile virus,Aedes caspius,mosquito,weak,review,Europe,Ferraguti 2024,natural infection; debated
West Nile virus,Culex modestus,mosquito,confirmed,review,Europe,Ferraguti 2024,primary vector in Europe
West Nile virus,Culex perexiguus,mosquito,confirmed,review,Europe,Ferraguti 2024,primary vector in Europe
West Nile virus,Culex pipiens,mosquito,confirmed,review,Europe,Ferraguti 2024,primary vector in Europe
```

Extraction notes

- Broad review on mosquito-pathogen specificity; partial support only, not a complete species list.
- I trimmed the initial agent output to the taxa the paper explicitly discusses for West Nile virus in the main text.
- `Aedes caspius` is the only clearly uncertain row because the paper says its competence remains debated.

## Goddard et al. - 2002 - Vector Competence of California Mosquitoes for West Nile virus.pdf

```csv
West Nile virus,Culex tarsalis,mosquito,probable,lab,California USA,Goddard et al. 2002,primary role
West Nile virus,Culex pipiens pipiens,mosquito,probable,lab,California USA,Goddard et al. 2002,enzootic vector candidate
West Nile virus,Culex pipiens quinquefasciatus,mosquito,candidate,lab,California USA,Goddard et al. 2002,variable populations
West Nile virus,Culex stigmatosoma,mosquito,probable,lab,California USA,Goddard et al. 2002,urban settings
West Nile virus,Culex erythrothorax,mosquito,probable,lab,California USA,Goddard et al. 2002,bridge potential
West Nile virus,Ochlerotatus dorsalis,mosquito,candidate,lab,California USA,Goddard et al. 2002,secondary role questioned
West Nile virus,Ochlerotatus melanimon,mosquito,candidate,lab,California USA,Goddard et al. 2002,secondary role questioned
West Nile virus,Ochlerotatus sierrensis,mosquito,weak,lab,California USA,Goddard et al. 2002,probably not vector
West Nile virus,Aedes vexans,mosquito,candidate,lab,California USA,Goddard et al. 2002,mammal-feeding
West Nile virus,Culiseta inornata,mosquito,candidate,lab,California USA,Goddard et al. 2002,minor role
```

Extraction notes

- Narrow laboratory competence study of 10 California mosquito species; not a complete species list.
- All rows come from direct experimental infection and transmission results, but I kept evidence levels below `confirmed` except where the paper clearly leans toward a major role.
- The most uncertain rows are the non-`Culex` species and `Ochlerotatus sierrensis`, which the paper explicitly downplays.

## Gray and Webb - 2014 - A review of the epidemiological and clinical aspects of West Nile virus.pdf

```csv
West Nile virus,Culex spp.,mosquito,confirmed,review,global,Gray and Webb 2014,primary vectors
West Nile virus,Aedes spp.,mosquito,candidate,review,,Gray and Webb 2014,lab only
West Nile virus,Culex tarsalis,mosquito,confirmed,review,North America,Gray and Webb 2014,major vector
West Nile virus,Culex modestus,mosquito,probable,review,Europe,Gray and Webb 2014,locally important vector
```

Extraction notes

- Broad epidemiological review; partial support only, not a complete species list.
- I kept the genus-level rows because this paper states the vector roles at that taxonomic level in the main text.
- `Aedes spp.` is the weakest row because the paper immediately cautions that ecological barriers likely limit outbreak importance in the field.

## Hernández-Triana et al. - 2014 - Emergence of West Nile Virus Lineage 2 in Europe A Review on the Introduction and Spread of a Mosqu.pdf

```csv
West Nile virus,Culex pipiens,mosquito,candidate,field,Italy; Greece,Hernandez-Triana et al. 2014,WNV-positive pools
West Nile virus,Culex modestus,mosquito,weak,field,Italy; Greece,Hernandez-Triana et al. 2014,to a lesser extent
West Nile virus,Aedes albopictus,mosquito,probable,lab,Europe,Hernandez-Triana et al. 2014,competent vector
```

Extraction notes

- Broad lineage-2 review; partial support only, not a complete species list.
- The two `Culex` rows are based on mosquito surveillance detections summarized in the paper, so I kept them conservative.
- `Aedes albopictus` is supported by explicit laboratory competence language and is the strongest species-level row here after the main `Culex` surveillance findings.

## Jansen et al. - 2013 - The Role of Australian Mosquito Species in the Transmission of Endemic and Exotic West Nile Virus St.pdf

```csv
West Nile virus,Culex annulirostris,mosquito,confirmed,field+lab,Australia,Jansen et al. 2013,principal vector
West Nile virus,Culex quinquefasciatus,mosquito,probable,field+lab,Australia,Jansen et al. 2013,primary enzootic candidate
West Nile virus,Culex gelidus,mosquito,candidate,lab,northern Australia,Jansen et al. 2013,introduced lab vector
West Nile virus,Culex squamosus,mosquito,candidate,field,northern Queensland,Jansen et al. 2013,regional candidate
West Nile virus,Culex sitiens,mosquito,weak,field+lab,coastal Australia,Jansen et al. 2013,poor lab vector; regional candidate
West Nile virus,Culex australicus,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Culex molestus,mosquito,candidate,review,Australia,Jansen et al. 2013,potential vector
West Nile virus,Aedes vigilax,mosquito,candidate,review,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes notoscriptus,mosquito,candidate,review,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes alternans,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Aedes normanensis,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Aedes tremulus,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Anopheles amictus,mosquito,weak,field,Australia,Jansen et al. 2013,field isolate
```

Extraction notes

- Broad Australia-focused review; partial support only, not a complete species list.
- `Culex annulirostris` is the clearest row because the paper explicitly treats it as the accepted primary WNVKUN vector and strongest candidate for exotic-strain transmission.
- The weakest rows are the occasional field-isolate and bridge-vector mentions, which I kept only at `candidate` or `weak`.

## Martinet et al. - 2019 - Mosquitoes of North-Western Europe as Potential Vectors of Arboviruses A Review.pdf

```csv
West Nile virus,Culex modestus,mosquito,confirmed,review,France,Martinet et al. 2019,historical vector in France
West Nile virus,Culex pipiens biotype pipiens,mosquito,probable,review,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex pipiens biotype molestus,mosquito,probable,review,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex torrentium,mosquito,probable,review,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Aedes detritus,mosquito,probable,review,United Kingdom,Martinet et al. 2019,competence demonstrated
West Nile virus,Anopheles plumbeus,mosquito,candidate,review,,Martinet et al. 2019,tested competent; role uncertain
West Nile virus,Aedes geniculatus,mosquito,candidate,review,,Martinet et al. 2019,tested competent; role uncertain
```

Extraction notes

- Broad north-western Europe review; partial support only, not a complete species list.
- I downgraded several review-summary rows from `confirmed` to `probable` where the support is competence-focused rather than a clearly established field role.
- `Anopheles plumbeus` and `Aedes geniculatus` are the most uncertain rows because the review says their role remains unknown.

## Martinet et al. - 2023 - Assessing vector competence of mosquitoes from northeastern France to West Nile virus and Usutu viru.pdf

```csv
West Nile virus,Culex pipiens,mosquito,probable,lab,northeastern France (Machault; Maine; Verzy),Martinet et al. 2023,main vector
West Nile virus,Aedes rusticus,mosquito,probable,lab,northeastern France (Berru),Martinet et al. 2023,new putative vector
West Nile virus,Aedes albopictus,mosquito,probable,lab,northeastern France (Strasbourg),Martinet et al. 2023,new putative vector
West Nile virus,Anopheles plumbeus,mosquito,weak,lab,northeastern France (Beaumont),Martinet et al. 2023,infected but did not transmit
```

Extraction notes

- Narrow experimental study from northeastern France; not a complete species list.
- `Aedes rusticus` and `Aedes albopictus` are retained as `probable` because the paper itself frames them as new putative vectors rather than settled established vectors.
- `Anopheles plumbeus` is included only as a weak negative-comparison row.

## Vogels et al. - 2017 - Vector competence of European mosquitoes for West Nile virus.pdf

```csv
West Nile virus,Aedes albopictus,mosquito,confirmed,review,Spain; Italy,Vogels et al. 2017,competent
West Nile virus,Aedes caspius,mosquito,weak,review,France,Vogels et al. 2017,not competent
West Nile virus,Aedes detritus,mosquito,confirmed,review,United Kingdom,Vogels et al. 2017,competent
West Nile virus,Aedes japonicus japonicus,mosquito,weak,review,Germany,Vogels et al. 2017,not competent
West Nile virus,Culex modestus,mosquito,confirmed,review,France,Vogels et al. 2017,efficient vector
West Nile virus,Culex torrentium,mosquito,candidate,review,Germany,Vogels et al. 2017,infection and dissemination only
West Nile virus,Culex pipiens,mosquito,confirmed,review,Europe,Vogels et al. 2017,most important vector
```

Extraction notes

- Broad review of European WNV competence studies; partial support only, not a complete species list.
- I kept the rows the paper explicitly summarizes in its conclusions and omitted more granular biotype rows from the initial agent output.
- `Aedes caspius`, `Aedes japonicus japonicus`, and `Culex torrentium` are the most uncertain rows because the paper emphasizes non-competence or incomplete competence evidence.

## ALL TOGETHER

```csv
West Nile virus,Culex pipiens s.l.,mosquito,confirmed,review,Europe,Bellini et al. 2014,major vector role
West Nile virus,Culex modestus,mosquito,probable,review,Europe,Bellini et al. 2014,regional role
West Nile virus,Culex univittatus,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Coquillettidia richiardii,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes cantans,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes caspius,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes excrucians,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Aedes vexans,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Anopheles maculipennis s.s.,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Anopheles atroparvus,mosquito,candidate,review,Europe,Bellini et al. 2014,field-collected example
West Nile virus,Culex pipiens,mosquito,probable,review,Europe,Engler et al. 2013,important vector species
West Nile virus,Culex modestus,mosquito,probable,review,Europe,Engler et al. 2013,important vector species
West Nile virus,Culex perexiguus,mosquito,candidate,field,Spain,Engler et al. 2013,WNV-positive pools
West Nile virus,Ochlerotatus caspius,mosquito,weak,field,Italy,Engler et al. 2013,WNV-positive pools
West Nile virus,Aedes caspius,mosquito,weak,review,Europe,Ferraguti 2024,natural infection; debated
West Nile virus,Culex modestus,mosquito,confirmed,review,Europe,Ferraguti 2024,primary vector in Europe
West Nile virus,Culex perexiguus,mosquito,confirmed,review,Europe,Ferraguti 2024,primary vector in Europe
West Nile virus,Culex pipiens,mosquito,confirmed,review,Europe,Ferraguti 2024,primary vector in Europe
West Nile virus,Culex tarsalis,mosquito,probable,lab,California USA,Goddard et al. 2002,primary role
West Nile virus,Culex pipiens pipiens,mosquito,probable,lab,California USA,Goddard et al. 2002,enzootic vector candidate
West Nile virus,Culex pipiens quinquefasciatus,mosquito,candidate,lab,California USA,Goddard et al. 2002,variable populations
West Nile virus,Culex stigmatosoma,mosquito,probable,lab,California USA,Goddard et al. 2002,urban settings
West Nile virus,Culex erythrothorax,mosquito,probable,lab,California USA,Goddard et al. 2002,bridge potential
West Nile virus,Ochlerotatus dorsalis,mosquito,candidate,lab,California USA,Goddard et al. 2002,secondary role questioned
West Nile virus,Ochlerotatus melanimon,mosquito,candidate,lab,California USA,Goddard et al. 2002,secondary role questioned
West Nile virus,Ochlerotatus sierrensis,mosquito,weak,lab,California USA,Goddard et al. 2002,probably not vector
West Nile virus,Aedes vexans,mosquito,candidate,lab,California USA,Goddard et al. 2002,mammal-feeding
West Nile virus,Culiseta inornata,mosquito,candidate,lab,California USA,Goddard et al. 2002,minor role
West Nile virus,Culex spp.,mosquito,confirmed,review,global,Gray and Webb 2014,primary vectors
West Nile virus,Aedes spp.,mosquito,candidate,review,,Gray and Webb 2014,lab only
West Nile virus,Culex tarsalis,mosquito,confirmed,review,North America,Gray and Webb 2014,major vector
West Nile virus,Culex modestus,mosquito,probable,review,Europe,Gray and Webb 2014,locally important vector
West Nile virus,Culex pipiens,mosquito,candidate,field,Italy; Greece,Hernandez-Triana et al. 2014,WNV-positive pools
West Nile virus,Culex modestus,mosquito,weak,field,Italy; Greece,Hernandez-Triana et al. 2014,to a lesser extent
West Nile virus,Aedes albopictus,mosquito,probable,lab,Europe,Hernandez-Triana et al. 2014,competent vector
West Nile virus,Culex annulirostris,mosquito,confirmed,field+lab,Australia,Jansen et al. 2013,principal vector
West Nile virus,Culex quinquefasciatus,mosquito,probable,field+lab,Australia,Jansen et al. 2013,primary enzootic candidate
West Nile virus,Culex gelidus,mosquito,candidate,lab,northern Australia,Jansen et al. 2013,introduced lab vector
West Nile virus,Culex squamosus,mosquito,candidate,field,northern Queensland,Jansen et al. 2013,regional candidate
West Nile virus,Culex sitiens,mosquito,weak,field+lab,coastal Australia,Jansen et al. 2013,poor lab vector; regional candidate
West Nile virus,Culex australicus,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Culex molestus,mosquito,candidate,review,Australia,Jansen et al. 2013,potential vector
West Nile virus,Aedes vigilax,mosquito,candidate,review,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes notoscriptus,mosquito,candidate,review,Australia,Jansen et al. 2013,bridge vector candidate
West Nile virus,Aedes alternans,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Aedes normanensis,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Aedes tremulus,mosquito,candidate,field,Australia,Jansen et al. 2013,field isolate; secondary vector
West Nile virus,Anopheles amictus,mosquito,weak,field,Australia,Jansen et al. 2013,field isolate
West Nile virus,Culex modestus,mosquito,confirmed,review,France,Martinet et al. 2019,historical vector in France
West Nile virus,Culex pipiens biotype pipiens,mosquito,probable,review,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex pipiens biotype molestus,mosquito,probable,review,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Culex torrentium,mosquito,probable,review,France; The Netherlands; Switzerland; Germany,Martinet et al. 2019,competent for lineages 1 and 2
West Nile virus,Aedes detritus,mosquito,probable,review,United Kingdom,Martinet et al. 2019,competence demonstrated
West Nile virus,Anopheles plumbeus,mosquito,candidate,review,,Martinet et al. 2019,tested competent; role uncertain
West Nile virus,Aedes geniculatus,mosquito,candidate,review,,Martinet et al. 2019,tested competent; role uncertain
West Nile virus,Culex pipiens,mosquito,probable,lab,northeastern France (Machault; Maine; Verzy),Martinet et al. 2023,main vector
West Nile virus,Aedes rusticus,mosquito,probable,lab,northeastern France (Berru),Martinet et al. 2023,new putative vector
West Nile virus,Aedes albopictus,mosquito,probable,lab,northeastern France (Strasbourg),Martinet et al. 2023,new putative vector
West Nile virus,Anopheles plumbeus,mosquito,weak,lab,northeastern France (Beaumont),Martinet et al. 2023,infected but did not transmit
West Nile virus,Aedes albopictus,mosquito,confirmed,review,Spain; Italy,Vogels et al. 2017,competent
West Nile virus,Aedes caspius,mosquito,weak,review,France,Vogels et al. 2017,not competent
West Nile virus,Aedes detritus,mosquito,confirmed,review,United Kingdom,Vogels et al. 2017,competent
West Nile virus,Aedes japonicus japonicus,mosquito,weak,review,Germany,Vogels et al. 2017,not competent
West Nile virus,Culex modestus,mosquito,confirmed,review,France,Vogels et al. 2017,efficient vector
West Nile virus,Culex torrentium,mosquito,candidate,review,Germany,Vogels et al. 2017,infection and dissemination only
West Nile virus,Culex pipiens,mosquito,confirmed,review,Europe,Vogels et al. 2017,most important vector
```
