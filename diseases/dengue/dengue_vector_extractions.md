# Dengue Vector Extractions By PDF

## Althouse et al. - 2012 - Synchrony of Sylvatic Dengue Isolations A Multi-Host, Multi-Vector SIR Model of Dengue Virus Transm.pdf

```csv
dengue,Aedes furcifer,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes taylori,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes luteocephalus,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes vittatus,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes aegypti,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
```

Extraction notes
- Broad modeling paper with a Senegal-specific vector summary; partial support only, not a complete dengue vector list.
- I kept only the five mosquito taxa explicitly named in the PDF as the main dengue vectors in Senegal.

## Barrera - 2025 - Surveillance and Control of Dengue Vectors in the United States and Territories.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,global,Barrera 2025,principal vector
dengue,Aedes albopictus,mosquito,probable,review,,Barrera 2025,secondary vector
dengue,Aedes polynesiensis,mosquito,probable,review,South Pacific Islands,Barrera 2025,transmits dengue
dengue,Aedes mediovittatus,mosquito,candidate,review,Puerto Rico and US Virgin Islands,Barrera 2025,potential vector
```

Extraction notes
- Broad review chapter on the United States and territories; partial support only, not a complete dengue vector list.
- `Aedes polynesiensis` was kept at `probable` rather than `confirmed` because the chapter says it transmits dengue but still frames `Aedes aegypti` as the main vector on the islands.
- `Aedes mediovittatus` is the most uncertain row because the chapter says it could potentially transmit DENV but has not been incriminated in outbreaks.

## Delrieu et al. - 2023 - Temperature and transmission of chikungunya, dengue, and Zika viruses A systematic review of experi.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,,Delrieu et al. 2023,major vector
dengue,Aedes albopictus,mosquito,confirmed,review,,Delrieu et al. 2023,
```

Extraction notes
- Broad multi-regional review, but taxonomically narrow; partial support only, not a complete dengue vector list.
- Only `Aedes aegypti` and `Aedes albopictus` are explicitly covered for dengue in this source.

## Doeurk et al. - 2024 - Review of dengue vectors in Cambodia distribution, bionomics, vector competence, control and insect.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,Cambodia,Doeurk et al. 2024,main vector
dengue,Aedes albopictus,mosquito,confirmed,review,Cambodia,Doeurk et al. 2024,main vector
dengue,Aedes malayensis,mosquito,candidate,review,Cambodia,Doeurk et al. 2024,potential vector
dengue,Aedes scutellaris,mosquito,candidate,review,Cambodia,Doeurk et al. 2024,potential vector
```

Extraction notes
- Broad Cambodia-focused review; partial support only, not a complete dengue vector list.
- `Aedes aegypti` and `Aedes albopictus` are explicitly described as the main vectors.
- `Aedes malayensis` and `Aedes scutellaris` are uncertain rows because the review treats them as potential vectors only.

## Islam et al. - 2021 - Production, Transmission, Pathogenesis, and Control of Dengue Virus A Literature-Based Undivided Pe.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,,Islam et al. 2021,principal vector
dengue,Aedes albopictus,mosquito,confirmed,review,,Islam et al. 2021,principal vector
dengue,Aedes polynesiensis,mosquito,probable,review,,Islam et al. 2021,secondary vector
dengue,Aedes niveus,mosquito,probable,review,,Islam et al. 2021,secondary vector
```

Extraction notes
- Broad review; partial support only, not a complete dengue vector list.
- I kept only the taxa explicitly tied to dengue in the main text and figure caption.
- `Aedes polynesiensis` and `Aedes niveus` are the uncertain rows because they are presented as secondary vectors in some regions rather than primary vectors.

## Jones et al. - 2020 - Arbovirus vectors of epidemiological concern in the Americas A scoping review of entomological stud.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Aedes albopictus,mosquito,confirmed,review,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Culex quinquefasciatus,mosquito,weak,review,Americas,Jones et al. 2020,field infection only not experimentally competent
```

Extraction notes
- Broad scoping review of the Americas; partial support only, not a complete dengue vector list.
- I omitted non-species entries such as `Aedes sp.`.
- `Culex quinquefasciatus` is the only uncertain row and is kept `weak` because the review reports field infection but explicitly says it was not experimentally competent.

## Kristan et al. - 2025 - Quantifying the potential relative roles of Aedes aegypti and Ae. albopictus in dengue.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,,Kristan et al. 2025,principal vector
dengue,Aedes albopictus,mosquito,probable,review,,Kristan et al. 2025,secondary vector
```

Extraction notes
- Broad multi-region systematic review/meta-analysis; partial support only, not a complete dengue vector list.
- Only `Aedes aegypti` and `Aedes albopictus` are explicitly named for dengue in this source.
- `Aedes albopictus` is the more uncertain row because the paper frames it as secondary rather than primary.

## Poole-Smith et al. - 2015 - Comparison of Vector Competence of Aedes mediovittatus and Aedes aegypti for Dengue Virus Implicati.pdf

```csv
dengue,Aedes mediovittatus,mosquito,probable,lab,Puerto Rico,Poole-Smith et al. 2015,secondary vector
dengue,Aedes aegypti,mosquito,confirmed,lab,Puerto Rico,Poole-Smith et al. 2015,primary vector
```

Extraction notes
- Narrow laboratory comparison from Puerto Rico; not a complete dengue vector list.
- `Aedes mediovittatus` is kept at `probable` because the paper supports competence and a possible secondary role, while `Aedes aegypti` is explicit as the primary vector.

## Rakotonirina et al. - 2022 - MALDI-TOF MS An effective tool for a global surveillance of dengue vector species.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,global,Rakotonirina et al. 2022,primary vector
dengue,Aedes albopictus,mosquito,confirmed,review,global,Rakotonirina et al. 2022,primary vector
dengue,Aedes polynesiensis,mosquito,confirmed,review,Pacific region,Rakotonirina et al. 2022,main vector
dengue,Aedes scutellaris,mosquito,candidate,review,Pacific region,Rakotonirina et al. 2022,highly suspected
dengue,Aedes pseudoscutellaris,mosquito,candidate,review,Pacific region,Rakotonirina et al. 2022,highly suspected
dengue,Aedes malayensis,mosquito,confirmed,review,,Rakotonirina et al. 2022,vector of DENV and CHIKV
```

Extraction notes
- Broad multi-regional paper; partial support only, not a complete dengue vector list.
- `Aedes scutellaris` and `Aedes pseudoscutellaris` are the most uncertain rows because the paper says these species are vectors, or highly suspected to be vectors, of DENV.
- `Aedes malayensis` is supported in a later limitation sentence rather than the main vector summary.

## Schaffner and Mathis - 2014 - Dengue and dengue vectors in the WHO European region past, present, and scenarios for the future.pdf

```csv
dengue,Aedes aegypti,mosquito,confirmed,review,WHO European region,Schaffner and Mathis 2014,principal vector
dengue,Aedes albopictus,mosquito,probable,review,WHO European region,Schaffner and Mathis 2014,main established vector
dengue,Aedes cretinus,mosquito,candidate,review,Greece and Turkey,Schaffner and Mathis 2014,potential vector
dengue,Aedes japonicus,mosquito,candidate,review,central Europe,Schaffner and Mathis 2014,potential vector
dengue,Aedes vittatus,mosquito,weak,review,India; western Mediterranean,Schaffner and Mathis 2014,occasionally mentioned; role likely restricted
```

Extraction notes
- Broad review of dengue vectors in the WHO European region; partial support only, not a complete dengue vector list.
- `Aedes albopictus` was downgraded to `probable` in the consolidated file because the paper calls it the main established vector in Europe but also describes it more broadly as the secondary vector globally.
- `Aedes cretinus`, `Aedes japonicus`, and `Aedes vittatus` are the uncertain rows because the paper frames them as potential or downplayed vectors.

## ALL TOGETHER

```csv
dengue,Aedes furcifer,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes taylori,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes luteocephalus,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes vittatus,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes aegypti,mosquito,confirmed,review,Senegal,Althouse et al. 2012,main vector in Senegal
dengue,Aedes aegypti,mosquito,confirmed,review,global,Barrera 2025,principal vector
dengue,Aedes albopictus,mosquito,probable,review,,Barrera 2025,secondary vector
dengue,Aedes polynesiensis,mosquito,probable,review,South Pacific Islands,Barrera 2025,transmits dengue
dengue,Aedes mediovittatus,mosquito,candidate,review,Puerto Rico and US Virgin Islands,Barrera 2025,potential vector
dengue,Aedes aegypti,mosquito,confirmed,review,,Delrieu et al. 2023,major vector
dengue,Aedes albopictus,mosquito,confirmed,review,,Delrieu et al. 2023,
dengue,Aedes aegypti,mosquito,confirmed,review,Cambodia,Doeurk et al. 2024,main vector
dengue,Aedes albopictus,mosquito,confirmed,review,Cambodia,Doeurk et al. 2024,main vector
dengue,Aedes malayensis,mosquito,candidate,review,Cambodia,Doeurk et al. 2024,potential vector
dengue,Aedes scutellaris,mosquito,candidate,review,Cambodia,Doeurk et al. 2024,potential vector
dengue,Aedes aegypti,mosquito,confirmed,review,,Islam et al. 2021,principal vector
dengue,Aedes albopictus,mosquito,confirmed,review,,Islam et al. 2021,principal vector
dengue,Aedes polynesiensis,mosquito,probable,review,,Islam et al. 2021,secondary vector
dengue,Aedes niveus,mosquito,probable,review,,Islam et al. 2021,secondary vector
dengue,Aedes aegypti,mosquito,confirmed,review,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Aedes albopictus,mosquito,confirmed,review,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Culex quinquefasciatus,mosquito,weak,review,Americas,Jones et al. 2020,field infection only not experimentally competent
dengue,Aedes aegypti,mosquito,confirmed,review,,Kristan et al. 2025,principal vector
dengue,Aedes albopictus,mosquito,probable,review,,Kristan et al. 2025,secondary vector
dengue,Aedes mediovittatus,mosquito,probable,lab,Puerto Rico,Poole-Smith et al. 2015,secondary vector
dengue,Aedes aegypti,mosquito,confirmed,lab,Puerto Rico,Poole-Smith et al. 2015,primary vector
dengue,Aedes aegypti,mosquito,confirmed,review,global,Rakotonirina et al. 2022,primary vector
dengue,Aedes albopictus,mosquito,confirmed,review,global,Rakotonirina et al. 2022,primary vector
dengue,Aedes polynesiensis,mosquito,confirmed,review,Pacific region,Rakotonirina et al. 2022,main vector
dengue,Aedes scutellaris,mosquito,candidate,review,Pacific region,Rakotonirina et al. 2022,highly suspected
dengue,Aedes pseudoscutellaris,mosquito,candidate,review,Pacific region,Rakotonirina et al. 2022,highly suspected
dengue,Aedes malayensis,mosquito,confirmed,review,,Rakotonirina et al. 2022,vector of DENV and CHIKV
dengue,Aedes aegypti,mosquito,confirmed,review,WHO European region,Schaffner and Mathis 2014,principal vector
dengue,Aedes albopictus,mosquito,probable,review,WHO European region,Schaffner and Mathis 2014,main established vector
dengue,Aedes cretinus,mosquito,candidate,review,Greece and Turkey,Schaffner and Mathis 2014,potential vector
dengue,Aedes japonicus,mosquito,candidate,review,central Europe,Schaffner and Mathis 2014,potential vector
dengue,Aedes vittatus,mosquito,weak,review,India; western Mediterranean,Schaffner and Mathis 2014,occasionally mentioned; role likely restricted
```
