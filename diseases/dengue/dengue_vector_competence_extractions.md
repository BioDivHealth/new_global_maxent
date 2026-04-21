# Dengue Vector Competence Extractions By PDF

## Althouse et al. - 2012 - Synchrony of Sylvatic Dengue Isolations A Multi-Host, Multi-Vector SIR Model of Dengue Virus Transm.pdf

```csv
dengue,Aedes furcifer,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes taylori,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes luteocephalus,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes vittatus,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes aegypti,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
```

Extraction notes
- Broad modeling paper with Senegal-focused background; partial support only.
- The source does not present new competence experiments; it summarizes the Senegal vector set and notes wide variation in competence estimates.

## Barrera - 2025 - Surveillance and Control of Dengue Vectors in the United States and Territories.pdf

```csv
dengue,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Barrera 2025,principal vector worldwide; high vector competence and vector capacity
dengue,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Barrera 2025,less efficient than Ae. aegypti; still involved in dengue transmission
dengue,Aedes polynesiensis,mosquito,competent,narrative_review,yes,,South Pacific Islands,Barrera 2025,can transmit dengue; main vector on the islands
dengue,Aedes mediovittatus,mosquito,unclear,narrative_review,,,Puerto Rico and US Virgin Islands,Barrera 2025,potential vector; not yet incriminated in actual outbreaks
```

Extraction notes
- Broad review chapter on the US and territories; partial support only.
- `Aedes mediovittatus` is the most uncertain row because the chapter describes it as a potential vector rather than a proven outbreak vector.

## Delrieu et al. - 2023 - Temperature and transmission of chikungunya, dengue, and Zika viruses A systematic review of experi.pdf

```csv
dengue,Aedes aegypti,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent competence; most studies show higher temperatures increase infection dissemination and transmission
dengue,Aedes albopictus,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent competence; most studies show higher temperatures increase infection dissemination and transmission
```

Extraction notes
- Broad multi-regional systematic review; taxonomically narrow for dengue because it only covers `Aedes aegypti` and `Aedes albopictus`.
- Mixed status reflects the review’s temperature-dependent results and the fact that several studies differed from the overall pattern.

## Doeurk et al. - 2024 - Review of dengue vectors in Cambodia distribution, bionomics, vector competence, control and insect.pdf

```csv
dengue,Aedes aegypti,mosquito,competent,field_plus_lab,yes,,Cambodia,Doeurk et al. 2024,main vector; 2023 Phnom Penh DENV-1 experiment showed infection dissemination and saliva positivity
dengue,Aedes albopictus,mosquito,competent,field_plus_lab,yes,,Cambodia,Doeurk et al. 2024,main vector; 2023 Phnom Penh DENV-1 experiment showed infection dissemination and saliva positivity
dengue,Aedes malayensis,mosquito,unclear,narrative_review,,,Cambodia,Doeurk et al. 2024,potential vector in Cambodia; cited as competent in Thailand Laos and Singapore
dengue,Aedes scutellaris,mosquito,unclear,narrative_review,,,Cambodia,Doeurk et al. 2024,potential vector in Cambodia; cited as competent in Thailand Laos and Singapore
```

Extraction notes
- Broad Cambodia-focused review with added unpublished competence data; partial support only.
- The two main vectors have direct saliva-positive lab evidence in the paper.
- `Aedes malayensis` and `Aedes scutellaris` remain uncertain because the review treats them as potential vectors rather than showing Cambodia-specific competence results.

## Islam et al. - 2021 - Production, Transmission, Pathogenesis, and Control of Dengue Virus A Literature-Based Undivided Pe.pdf

```csv
dengue,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Islam et al. 2021,primary vector
dengue,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Islam et al. 2021,primary vector
dengue,Aedes polynesiensis,mosquito,unclear,narrative_review,,,global,Islam et al. 2021,secondary vector in some regions
dengue,Aedes niveus,mosquito,unclear,narrative_review,,,global,Islam et al. 2021,secondary vector in some regions
```

Extraction notes
- Broad literature review; partial support only.
- The secondary-vector rows are kept uncertain because the source states role language rather than presenting direct competence experiments.

## Jones et al. - 2020 - Arbovirus vectors of epidemiological concern in the Americas A scoping review of entomological stud.pdf

```csv
dengue,Aedes aegypti,mosquito,competent,systematic_review,yes,yes,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Aedes albopictus,mosquito,competent,systematic_review,yes,yes,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Culex quinquefasciatus,mosquito,not_competent,systematic_review,,yes,Americas,Jones et al. 2020,field infection observed but not identified as experimentally competent
```

Extraction notes
- Broad scoping review of the Americas; partial support only.
- `Culex quinquefasciatus` is the negative row because the review reports field infection but says it was not experimentally competent.

## Kristan et al. - 2025 - Quantifying the potential relative roles of Aedes aegypti and Ae. albopictus in dengue.pdf

```csv
dengue,Aedes aegypti,mosquito,unclear,systematic_review,,yes,global,Kristan et al. 2025,field DENV prevalence meta-analysis; primary vector role but no direct competence test
dengue,Aedes albopictus,mosquito,mixed,systematic_review,,yes,global,Kristan et al. 2025,field DENV prevalence meta-analysis; overall lower prevalence than Ae. aegypti but post-2000 differences were not significant
```

Extraction notes
- Systematic review and meta-analysis of field DENV prevalence; this is proxy evidence rather than direct competence testing.
- `Aedes albopictus` is mixed because the paper’s overall result differs from the post-2000 subgroup and the authors stress context dependence.

## Poole-Smith et al. - 2015 - Comparison of Vector Competence of Aedes mediovittatus and Aedes aegypti for Dengue Virus Implicati.pdf

```csv
dengue,Aedes mediovittatus,mosquito,mixed,lab_experiment,yes,,Puerto Rico,Poole-Smith et al. 2015,competent for DENV-1 to -3; DENV-4 transmission lower than Ae. aegypti
dengue,Aedes aegypti,mosquito,competent,lab_experiment,yes,,Puerto Rico,Poole-Smith et al. 2015,primary vector; higher DENV-4 infection and transmission than Ae. mediovittatus
```

Extraction notes
- Narrow laboratory comparison from Puerto Rico; strong direct competence evidence for both species.
- `Aedes mediovittatus` is mixed because competence varies by serotype and DENV-4 performance was lower than `Ae. aegypti`.

## Rakotonirina et al. - 2022 - MALDI-TOF MS An effective tool for a global surveillance of dengue vector species.pdf

```csv
dengue,Aedes aegypti,mosquito,competent,narrative_review,,,global,Rakotonirina et al. 2022,primary vector of DENV
dengue,Aedes albopictus,mosquito,competent,narrative_review,,,global,Rakotonirina et al. 2022,primary vector of DENV
dengue,Aedes polynesiensis,mosquito,competent,narrative_review,,,Pacific region,Rakotonirina et al. 2022,main vector in the Pacific region
dengue,Aedes scutellaris,mosquito,unclear,narrative_review,,,Pacific region,Rakotonirina et al. 2022,vector or highly suspected vector of DENV
dengue,Aedes pseudoscutellaris,mosquito,unclear,narrative_review,,,Pacific region,Rakotonirina et al. 2022,vector or highly suspected vector of DENV
dengue,Aedes malayensis,mosquito,competent,narrative_review,,,,Rakotonirina et al. 2022,vector of DENV and CHIKV
```

Extraction notes
- Broad identification paper with review-style background; partial support only.
- The `Aedes scutellaris` and `Aedes pseudoscutellaris` rows are the most uncertain because the source frames them as vectors or highly suspected vectors.

## Schaffner and Mathis - 2014 - Dengue and dengue vectors in the WHO European region past, present, and scenarios for the future.pdf

```csv
dengue,Aedes aegypti,mosquito,competent,narrative_review,yes,,WHO European region,Schaffner and Mathis 2014,principal urban vector worldwide; historical transmission in Europe
dengue,Aedes albopictus,mosquito,competent,narrative_review,yes,,WHO European region,Schaffner and Mathis 2014,main established vector in Europe; secondary vector worldwide
dengue,Aedes cretinus,mosquito,unclear,narrative_review,,,Greece and Turkey,Schaffner and Mathis 2014,potential contributor if competent; none subjected to experimental infection
dengue,Aedes japonicus,mosquito,competent,narrative_review,,,central Europe,Schaffner and Mathis 2014,high vector competence in laboratory; potential dengue vector
dengue,Aedes vittatus,mosquito,unclear,narrative_review,,,India; western Mediterranean,Schaffner and Mathis 2014,occasionally mentioned as a dengue vector in India; role likely restricted in Europe
```

Extraction notes
- Broad WHO European region review; partial support only.
- `Aedes cretinus` and `Aedes vittatus` are uncertain because the paper mentions them as possible or occasional vectors without direct experimental support in the review itself.

## ALL TOGETHER

```csv
dengue,Aedes furcifer,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes taylori,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes luteocephalus,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes vittatus,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes aegypti,mosquito,mixed,narrative_review,,,Senegal,Althouse et al. 2012,main vector in Senegal; studies cited as differing widely in vector competence
dengue,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Barrera 2025,principal vector worldwide; high vector competence and vector capacity
dengue,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Barrera 2025,less efficient than Ae. aegypti; still involved in dengue transmission
dengue,Aedes polynesiensis,mosquito,competent,narrative_review,yes,,South Pacific Islands,Barrera 2025,can transmit dengue; main vector on the islands
dengue,Aedes mediovittatus,mosquito,unclear,narrative_review,,,Puerto Rico and US Virgin Islands,Barrera 2025,potential vector; not yet incriminated in actual outbreaks
dengue,Aedes aegypti,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent competence; most studies show higher temperatures increase infection dissemination and transmission
dengue,Aedes albopictus,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent competence; most studies show higher temperatures increase infection dissemination and transmission
dengue,Aedes aegypti,mosquito,competent,field_plus_lab,yes,,Cambodia,Doeurk et al. 2024,main vector; 2023 Phnom Penh DENV-1 experiment showed infection dissemination and saliva positivity
dengue,Aedes albopictus,mosquito,competent,field_plus_lab,yes,,Cambodia,Doeurk et al. 2024,main vector; 2023 Phnom Penh DENV-1 experiment showed infection dissemination and saliva positivity
dengue,Aedes malayensis,mosquito,unclear,narrative_review,,,Cambodia,Doeurk et al. 2024,potential vector in Cambodia; cited as competent in Thailand Laos and Singapore
dengue,Aedes scutellaris,mosquito,unclear,narrative_review,,,Cambodia,Doeurk et al. 2024,potential vector in Cambodia; cited as competent in Thailand Laos and Singapore
dengue,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Islam et al. 2021,primary vector
dengue,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Islam et al. 2021,primary vector
dengue,Aedes polynesiensis,mosquito,unclear,narrative_review,,,global,Islam et al. 2021,secondary vector in some regions
dengue,Aedes niveus,mosquito,unclear,narrative_review,,,global,Islam et al. 2021,secondary vector in some regions
dengue,Aedes aegypti,mosquito,competent,systematic_review,yes,yes,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Aedes albopictus,mosquito,competent,systematic_review,yes,yes,Americas,Jones et al. 2020,experimentally competent and field infection observed
dengue,Culex quinquefasciatus,mosquito,not_competent,systematic_review,,yes,Americas,Jones et al. 2020,field infection observed but not identified as experimentally competent
dengue,Aedes aegypti,mosquito,unclear,systematic_review,,yes,global,Kristan et al. 2025,field DENV prevalence meta-analysis; primary vector role but no direct competence test
dengue,Aedes albopictus,mosquito,mixed,systematic_review,,yes,global,Kristan et al. 2025,field DENV prevalence meta-analysis; overall lower prevalence than Ae. aegypti but post-2000 differences were not significant
dengue,Aedes mediovittatus,mosquito,mixed,lab_experiment,yes,,Puerto Rico,Poole-Smith et al. 2015,competent for DENV-1 to -3; DENV-4 transmission lower than Ae. aegypti
dengue,Aedes aegypti,mosquito,competent,lab_experiment,yes,,Puerto Rico,Poole-Smith et al. 2015,primary vector; higher DENV-4 infection and transmission than Ae. mediovittatus
dengue,Aedes aegypti,mosquito,competent,narrative_review,,,global,Rakotonirina et al. 2022,primary vector of DENV
dengue,Aedes albopictus,mosquito,competent,narrative_review,,,global,Rakotonirina et al. 2022,primary vector of DENV
dengue,Aedes polynesiensis,mosquito,competent,narrative_review,,,Pacific region,Rakotonirina et al. 2022,main vector in the Pacific region
dengue,Aedes scutellaris,mosquito,unclear,narrative_review,,,Pacific region,Rakotonirina et al. 2022,vector or highly suspected vector of DENV
dengue,Aedes pseudoscutellaris,mosquito,unclear,narrative_review,,,Pacific region,Rakotonirina et al. 2022,vector or highly suspected vector of DENV
dengue,Aedes malayensis,mosquito,competent,narrative_review,,,,Rakotonirina et al. 2022,vector of DENV and CHIKV
dengue,Aedes aegypti,mosquito,competent,narrative_review,yes,,WHO European region,Schaffner and Mathis 2014,principal urban vector worldwide; historical transmission in Europe
dengue,Aedes albopictus,mosquito,competent,narrative_review,yes,,WHO European region,Schaffner and Mathis 2014,main established vector in Europe; secondary vector worldwide
dengue,Aedes cretinus,mosquito,unclear,narrative_review,,,Greece and Turkey,Schaffner and Mathis 2014,potential contributor if competent; none subjected to experimental infection
dengue,Aedes japonicus,mosquito,competent,narrative_review,,,central Europe,Schaffner and Mathis 2014,high vector competence in laboratory; potential dengue vector
dengue,Aedes vittatus,mosquito,unclear,narrative_review,,,India; western Mediterranean,Schaffner and Mathis 2014,occasionally mentioned as a dengue vector in India; role likely restricted in Europe
```
