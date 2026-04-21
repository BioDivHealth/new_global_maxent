# Zika Vector Competence Extractions By PDF

## Bisia et al. - 2023 - Secondary vectors of Zika Virus, a systematic review of laboratory vector competence studies.pdf

```csv
zika,Aedes aegypti,mosquito,competent,systematic_review,yes,,global,Bisia et al. 2023,primary vector
zika,Aedes albopictus,mosquito,competent,systematic_review,yes,,global,Bisia et al. 2023,established secondary vector
zika,Aedes japonicus,mosquito,mixed,systematic_review,yes,,temperate regions,Bisia et al. 2023,transmission at higher temperature
zika,Aedes detritus,mosquito,mixed,systematic_review,yes,,temperate regions,Bisia et al. 2023,transmission at higher temperature
zika,Aedes vexans,mosquito,mixed,systematic_review,yes,,temperate regions,Bisia et al. 2023,transmission at higher temperature
zika,Aedes notoscriptus,mosquito,mixed,systematic_review,yes,,Australia,Bisia et al. 2023,local Australian vector; temperature-sensitive
zika,Aedes camptorhynchus,mosquito,mixed,systematic_review,yes,,Australia,Bisia et al. 2023,local Australian vector; temperature-sensitive
zika,Culex quinquefasciatus,mosquito,not_competent,systematic_review,no,,global,Bisia et al. 2023,not found competent
```

Extraction notes

- Broad systematic review of laboratory vector competence studies; partial support only.
- I kept the paper's explicit secondary-vector conclusions and the clear negative assessment for `Culex quinquefasciatus`.
- The temperature-dependent species are marked `mixed` because the review ties them to higher-temperature transmission, not uniform competence across conditions.

## Delrieu et al. - 2023 - Temperature and transmission of chikungunya, dengue, and Zika viruses A systematic review of experi.pdf

```csv
zika,Aedes aegypti,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent ZIKV transmission across studies
zika,Aedes albopictus,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent ZIKV transmission across studies
```

Extraction notes

- Broad systematic review, but the Zika content is narrowly centered on `Aedes aegypti` and `Aedes albopictus`.
- I used `mixed` because the review explicitly reports higher-temperature increases in many studies but also lower/no-effect findings in others.

## Dodson et al. - 2018 - Vector competence of selected North American Anopheles and Culex mosquitoes for Zika virus.pdf

```csv
zika,Aedes aegypti,mosquito,competent,lab_experiment,yes,,United States,Dodson et al. 2018,positive control; infected, disseminated, and transmitted
zika,Anopheles freeborni,mosquito,not_competent,lab_experiment,no,,United States,Dodson et al. 2018,unable to be infected
zika,Anopheles quadrimaculatus,mosquito,not_competent,lab_experiment,no,,United States,Dodson et al. 2018,unable to be infected
zika,Culex tarsalis,mosquito,not_competent,lab_experiment,no,,United States,Dodson et al. 2018,unable to be infected
```

Extraction notes

- Narrow four-species laboratory study from North America.
- The negative rows are explicit: the three comparison species were unable to become infected, while `Aedes aegypti` served as the transmitting positive control.

## Epelboin et al. - 2017 - Zika virus An updated review of competent or naturally infected mosquitoes.pdf

```csv
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,yes,global,Epelboin et al. 2017,main vector; field infection also reported
zika,Aedes albopictus,mosquito,competent,narrative_review,yes,yes,global,Epelboin et al. 2017,competent and field infected
zika,Aedes africanus,mosquito,unclear,narrative_review,,yes,Africa,Epelboin et al. 2017,field infection only
zika,Aedes hensilli,mosquito,unclear,narrative_review,,yes,Yap,Epelboin et al. 2017,field association; no transmission confirmed
zika,Aedes luteocephalus,mosquito,competent,narrative_review,yes,yes,Senegal,Epelboin et al. 2017,saliva positive in review summary
zika,Aedes vittatus,mosquito,competent,narrative_review,yes,yes,Senegal,Epelboin et al. 2017,saliva positive in review summary
zika,Aedes notoscriptus,mosquito,mixed,narrative_review,mixed,,Australia,Epelboin et al. 2017,one positive and one negative study in review
zika,Culex quinquefasciatus,mosquito,mixed,narrative_review,mixed,yes,global,Epelboin et al. 2017,controversial; some positive, many negative studies
zika,Culex pipiens,mosquito,not_competent,narrative_review,no,,global,Epelboin et al. 2017,review treats as refractory overall
```

Extraction notes

- Broad review covering both laboratory competence and natural infection.
- I kept only the clearest source-grounded competence rows and the explicit controversy/negative signal for `Culex quinquefasciatus` and `Culex pipiens`.
- Several other mosquitoes are mentioned in the paper as field detections only; I left those out here to avoid inflating a competence list from natural infection alone.

## Evans et al. - 2017 - Data-driven identification of potential Zika virus vectors.pdf

```csv
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Evans et al. 2017,known vector in Table 1
zika,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Evans et al. 2017,known vector in Table 1
zika,Aedes furcifer,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes vittatus,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes taylori,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes luteocephalus,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes hensilli,mosquito,competent,narrative_review,yes,,Yap,Evans et al. 2017,known vector in Table 1
zika,Culex quinquefasciatus,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
zika,Culex pipiens,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
zika,Psorophora ferox,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
zika,Runchomyia frontosa,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
```

Extraction notes

- Model-based prediction paper, not a primary competence experiment.
- I treated the Table 1 `Yes` entries as source-supported known vectors and kept several model-predicted taxa as `unclear` because the paper itself does not establish competence for them.

## Gardner et al. - 2017 - Vector status of Aedes species determines geographical risk of autochthonous Zika virus establishmen.pdf

```csv
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Gardner et al. 2017,primary vector
zika,Aedes albopictus,mosquito,unclear,narrative_review,,,global,Gardner et al. 2017,relative competence treated as scenario uncertainty
```

Extraction notes

- This is a risk-model paper, so it is not a direct competence assay.
- The paper clearly treats `Aedes aegypti` as the competent baseline vector and `Aedes albopictus` as the uncertain secondary-vector scenario.

## Gutiérrez-Bugallo et al. - 2019 - Vector-borne transmission and evolution of Zika virus.pdf

```csv
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,yes,global,Gutiérrez-Bugallo et al. 2019,urban cycle
zika,Aedes albopictus,mosquito,mixed,narrative_review,mixed,yes,global,Gutiérrez-Bugallo et al. 2019,potential major vector; mixed evidence across studies
zika,Aedes vexans,mosquito,mixed,narrative_review,mixed,yes,Mexico,Gutiérrez-Bugallo et al. 2019,field detection and lab transmission
zika,Aedes notoscriptus,mosquito,mixed,narrative_review,yes,,Australia,Gutiérrez-Bugallo et al. 2019,lab transmission in review summary
zika,Aedes camptorhynchus,mosquito,mixed,narrative_review,yes,,Australia,Gutiérrez-Bugallo et al. 2019,lab transmission in review summary
zika,Aedes hensilli,mosquito,unclear,narrative_review,,yes,Yap Island,Gutiérrez-Bugallo et al. 2019,outbreak-associated; no confirmed transmission
zika,Aedes polynesiensis,mosquito,unclear,narrative_review,,yes,French Polynesia,Gutiérrez-Bugallo et al. 2019,outbreak-associated; no confirmed transmission
zika,Culex quinquefasciatus,mosquito,mixed,narrative_review,mixed,yes,global,Gutiérrez-Bugallo et al. 2019,debated role
zika,Culex coronator,mosquito,unclear,narrative_review,,yes,Mexico,Gutiérrez-Bugallo et al. 2019,field detection only
zika,Culex tarsalis,mosquito,unclear,narrative_review,,yes,Mexico,Gutiérrez-Bugallo et al. 2019,field detection only
zika,Anopheles gambiae s.l.,mosquito,unclear,narrative_review,,yes,Africa,Gutiérrez-Bugallo et al. 2019,infected only once
zika,Aedes apicoargenteus,mosquito,unclear,narrative_review,,yes,global,Gutiérrez-Bugallo et al. 2019,infected only once
```

Extraction notes

- Broad review; partial support only.
- I kept the review's main urban and outbreak-associated taxa, plus the explicit field-detection rows that matter for competence interpretation.
- `Culex quinquefasciatus` remains mixed because the review describes a debated literature rather than a settled negative.

## Huang et al. - 2016 - Culex Species Mosquitoes and Zika Virus.pdf

```csv
zika,Culex pipiens,mosquito,not_competent,lab_experiment,no,,United States,Huang et al. 2016,refractory in California and New Jersey populations
zika,Culex quinquefasciatus,mosquito,not_competent,lab_experiment,no,,United States,Huang et al. 2016,refractory in Florida population
```

Extraction notes

- Narrow laboratory study focused on `Culex` populations from the United States.
- The paper is strongly negative and does not support a Zika transmission role for the tested populations.

## Ledermann et al. - 2014 - Aedes hensilli as a Potential Vector of Chikungunya and Zika Viruses.pdf

```csv
zika,Aedes hensilli,mosquito,mixed,lab_experiment,no,no,Yap Island,Ledermann et al. 2014,high infection and dissemination but no field virus isolate and no transmission demonstrated
```

Extraction notes

- Narrow Yap outbreak paper with complementary laboratory infections.
- I kept `Aedes hensilli` as mixed because the study supports infection and dissemination, but it does not demonstrate transmission and the field material itself was negative.

## Lourenço-de-Oliveira and Failloux - 2017 - Lessons learned on Zika virus vectors.pdf

```csv
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Lourenço-de-Oliveira and Failloux 2017,main vector
zika,Culex pipiens,mosquito,not_competent,narrative_review,no,,global,Lourenço-de-Oliveira and Failloux 2017,review says populations proved incompetent
zika,Culex quinquefasciatus,mosquito,not_competent,narrative_review,no,,global,Lourenço-de-Oliveira and Failloux 2017,review says populations proved incompetent
zika,Aedes hensilli,mosquito,mixed,narrative_review,mixed,,Yap,Lourenço-de-Oliveira and Failloux 2017,disseminated infections but no confirmed transmission
zika,Aedes polynesiensis,mosquito,mixed,narrative_review,mixed,,French Polynesia,Lourenço-de-Oliveira and Failloux 2017,weak competence; no infectious saliva
```

Extraction notes

- Short viewpoint-style review, not a primary experiment.
- The review is very explicit that domestic `Culex` populations studied to date do not transmit ZIKV, while the two Pacific `Aedes` taxa remain incomplete/weak rather than fully established vectors.

## Richard et al. - 2016 - Vector Competence of French Polynesian Aedes aegypti and Aedes polynesiensis for Zika Virus.pdf

```csv
zika,Aedes aegypti,mosquito,competent,lab_experiment,yes,,French Polynesia,Richard et al. 2016,late transmission but infectious saliva detected
zika,Aedes polynesiensis,mosquito,mixed,lab_experiment,no,,French Polynesia,Richard et al. 2016,weak competence; no infectious saliva at sampled time points
```

Extraction notes

- Narrow two-species laboratory study.
- `Aedes aegypti` is clearly positive, while `Aedes polynesiensis` remains weak/mixed because infection and dissemination were observed but infectious saliva was not.

## Schulz and Becker - 2018 - Mosquitoes as Arbovirus Vectors From Species Identification to Vector Competence.pdf

```csv
zika,Aedes albopictus,mosquito,mixed,narrative_review,mixed,,Europe,Schulz and Becker 2018,temperature-dependent competence; 18 C not competent, 27 C competent
zika,Culex pipiens,mosquito,not_competent,narrative_review,no,,Italy,Schulz and Becker 2018,not competent at tested temperatures
zika,Culex pipiens biotype pipiens,mosquito,not_competent,narrative_review,no,,Germany,Schulz and Becker 2018,not competent at tested temperatures
zika,Culex pipiens biotype molestus,mosquito,not_competent,narrative_review,no,,Germany,Schulz and Becker 2018,not competent at tested temperatures
zika,Culex torrentium,mosquito,not_competent,narrative_review,no,,Germany,Schulz and Becker 2018,not competent at tested temperatures
```

<oai-mem-citation>
<citation_entries>
MEMORY.md:446-455|note=[kept competence and network evidence separate]
MEMORY.md:582-590|note=[preserved negative evidence and avoided field-only overclaiming]
MEMORY.md:619-619|note=[avoided overstating role labels from review language]
</citation_entries>
<rollout_ids>
</rollout_ids>
</oai-mem-citation>

Extraction notes

- Broad review chapter; only the explicit Zika competence statements are kept.
- `Aedes albopictus` is conditional on temperature, while the `Culex` taxa are explicit negative rows.

## s41559-019-0836-z.pdf

```csv
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,yes,global,Gutierrez-Bugallo et al. 2019,urban cycle
zika,Aedes albopictus,mosquito,mixed,narrative_review,mixed,yes,global,Gutierrez-Bugallo et al. 2019,potential major vector; mixed evidence across studies
zika,Aedes vexans,mosquito,mixed,narrative_review,mixed,yes,Mexico,Gutierrez-Bugallo et al. 2019,field detection and lab transmission
zika,Aedes notoscriptus,mosquito,mixed,narrative_review,yes,,Australia,Gutierrez-Bugallo et al. 2019,lab transmission in review summary
zika,Aedes camptorhynchus,mosquito,mixed,narrative_review,yes,,Australia,Gutierrez-Bugallo et al. 2019,lab transmission in review summary
zika,Aedes hensilli,mosquito,unclear,narrative_review,,yes,Yap Island,Gutierrez-Bugallo et al. 2019,outbreak-associated; no confirmed transmission
zika,Aedes polynesiensis,mosquito,unclear,narrative_review,,yes,French Polynesia,Gutierrez-Bugallo et al. 2019,outbreak-associated; no confirmed transmission
zika,Culex quinquefasciatus,mosquito,mixed,narrative_review,mixed,yes,global,Gutierrez-Bugallo et al. 2019,debated role
zika,Culex coronator,mosquito,unclear,narrative_review,,yes,Mexico,Gutierrez-Bugallo et al. 2019,field detection only
zika,Culex tarsalis,mosquito,unclear,narrative_review,,yes,Mexico,Gutierrez-Bugallo et al. 2019,field detection only
zika,Anopheles gambiae s.l.,mosquito,unclear,narrative_review,,yes,Africa,Gutierrez-Bugallo et al. 2019,infected only once
zika,Aedes apicoargenteus,mosquito,unclear,narrative_review,,yes,global,Gutierrez-Bugallo et al. 2019,infected only once
```

Extraction notes

- Duplicate underlying study for `Gutiérrez-Bugallo et al. 2019`.
- I kept the same row set here for source coverage, but these rows are omitted from `## ALL TOGETHER` to avoid duplicating the underlying study.

## ALL TOGETHER

```csv
zika,Aedes aegypti,mosquito,competent,systematic_review,yes,,global,Bisia et al. 2023,primary vector
zika,Aedes albopictus,mosquito,competent,systematic_review,yes,,global,Bisia et al. 2023,established secondary vector
zika,Aedes japonicus,mosquito,mixed,systematic_review,yes,,temperate regions,Bisia et al. 2023,transmission at higher temperature
zika,Aedes detritus,mosquito,mixed,systematic_review,yes,,temperate regions,Bisia et al. 2023,transmission at higher temperature
zika,Aedes vexans,mosquito,mixed,systematic_review,yes,,temperate regions,Bisia et al. 2023,transmission at higher temperature
zika,Aedes notoscriptus,mosquito,mixed,systematic_review,yes,,Australia,Bisia et al. 2023,local Australian vector; temperature-sensitive
zika,Aedes camptorhynchus,mosquito,mixed,systematic_review,yes,,Australia,Bisia et al. 2023,local Australian vector; temperature-sensitive
zika,Culex quinquefasciatus,mosquito,not_competent,systematic_review,no,,global,Bisia et al. 2023,not found competent
zika,Aedes aegypti,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent ZIKV transmission across studies
zika,Aedes albopictus,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,temperature-dependent ZIKV transmission across studies
zika,Aedes aegypti,mosquito,competent,lab_experiment,yes,,United States,Dodson et al. 2018,positive control; infected, disseminated, and transmitted
zika,Anopheles freeborni,mosquito,not_competent,lab_experiment,no,,United States,Dodson et al. 2018,unable to be infected
zika,Anopheles quadrimaculatus,mosquito,not_competent,lab_experiment,no,,United States,Dodson et al. 2018,unable to be infected
zika,Culex tarsalis,mosquito,not_competent,lab_experiment,no,,United States,Dodson et al. 2018,unable to be infected
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,yes,global,Epelboin et al. 2017,main vector; field infection also reported
zika,Aedes albopictus,mosquito,competent,narrative_review,yes,yes,global,Epelboin et al. 2017,competent and field infected
zika,Aedes africanus,mosquito,unclear,narrative_review,,yes,Africa,Epelboin et al. 2017,field infection only
zika,Aedes hensilli,mosquito,unclear,narrative_review,,yes,Yap,Epelboin et al. 2017,field association; no transmission confirmed
zika,Aedes luteocephalus,mosquito,competent,narrative_review,yes,yes,Senegal,Epelboin et al. 2017,saliva positive in review summary
zika,Aedes vittatus,mosquito,competent,narrative_review,yes,yes,Senegal,Epelboin et al. 2017,saliva positive in review summary
zika,Aedes notoscriptus,mosquito,mixed,narrative_review,mixed,,Australia,Epelboin et al. 2017,one positive and one negative study in review
zika,Culex quinquefasciatus,mosquito,mixed,narrative_review,mixed,yes,global,Epelboin et al. 2017,controversial; some positive, many negative studies
zika,Culex pipiens,mosquito,not_competent,narrative_review,no,,global,Epelboin et al. 2017,review treats as refractory overall
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Evans et al. 2017,known vector in Table 1
zika,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Evans et al. 2017,known vector in Table 1
zika,Aedes furcifer,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes vittatus,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes taylori,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes luteocephalus,mosquito,competent,narrative_review,yes,,Africa,Evans et al. 2017,known vector in Table 1
zika,Aedes hensilli,mosquito,competent,narrative_review,yes,,Yap,Evans et al. 2017,known vector in Table 1
zika,Culex quinquefasciatus,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
zika,Culex pipiens,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
zika,Psorophora ferox,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
zika,Runchomyia frontosa,mosquito,unclear,narrative_review,,no,global,Evans et al. 2017,GBM predicted; not a known vector
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Gardner et al. 2017,primary vector
zika,Aedes albopictus,mosquito,unclear,narrative_review,,,global,Gardner et al. 2017,relative competence treated as scenario uncertainty
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,yes,global,Gutierrez-Bugallo et al. 2019,urban cycle
zika,Aedes albopictus,mosquito,mixed,narrative_review,mixed,yes,global,Gutierrez-Bugallo et al. 2019,potential major vector; mixed evidence across studies
zika,Aedes vexans,mosquito,mixed,narrative_review,mixed,yes,Mexico,Gutierrez-Bugallo et al. 2019,field detection and lab transmission
zika,Aedes notoscriptus,mosquito,mixed,narrative_review,yes,,Australia,Gutierrez-Bugallo et al. 2019,lab transmission in review summary
zika,Aedes camptorhynchus,mosquito,mixed,narrative_review,yes,,Australia,Gutierrez-Bugallo et al. 2019,lab transmission in review summary
zika,Aedes hensilli,mosquito,unclear,narrative_review,,yes,Yap Island,Gutierrez-Bugallo et al. 2019,outbreak-associated; no confirmed transmission
zika,Aedes polynesiensis,mosquito,unclear,narrative_review,,yes,French Polynesia,Gutierrez-Bugallo et al. 2019,outbreak-associated; no confirmed transmission
zika,Culex quinquefasciatus,mosquito,mixed,narrative_review,mixed,yes,global,Gutierrez-Bugallo et al. 2019,debated role
zika,Culex coronator,mosquito,unclear,narrative_review,,yes,Mexico,Gutierrez-Bugallo et al. 2019,field detection only
zika,Culex tarsalis,mosquito,unclear,narrative_review,,yes,Mexico,Gutierrez-Bugallo et al. 2019,field detection only
zika,Anopheles gambiae s.l.,mosquito,unclear,narrative_review,,yes,Africa,Gutierrez-Bugallo et al. 2019,infected only once
zika,Aedes apicoargenteus,mosquito,unclear,narrative_review,,yes,global,Gutierrez-Bugallo et al. 2019,infected only once
zika,Culex pipiens,mosquito,not_competent,lab_experiment,no,,United States,Huang et al. 2016,refractory in California and New Jersey populations
zika,Culex quinquefasciatus,mosquito,not_competent,lab_experiment,no,,United States,Huang et al. 2016,refractory in Florida population
zika,Aedes hensilli,mosquito,mixed,lab_experiment,no,no,Yap Island,Ledermann et al. 2014,high infection and dissemination but no field virus isolate and no transmission demonstrated
zika,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Lourenço-de-Oliveira and Failloux 2017,main vector
zika,Culex pipiens,mosquito,not_competent,narrative_review,no,,global,Lourenço-de-Oliveira and Failloux 2017,review says populations proved incompetent
zika,Culex quinquefasciatus,mosquito,not_competent,narrative_review,no,,global,Lourenço-de-Oliveira and Failloux 2017,review says populations proved incompetent
zika,Aedes hensilli,mosquito,mixed,narrative_review,mixed,,Yap,Lourenço-de-Oliveira and Failloux 2017,disseminated infections but no confirmed transmission
zika,Aedes polynesiensis,mosquito,mixed,narrative_review,mixed,,French Polynesia,Lourenço-de-Oliveira and Failloux 2017,weak competence; no infectious saliva
zika,Aedes aegypti,mosquito,competent,lab_experiment,yes,,French Polynesia,Richard et al. 2016,late transmission but infectious saliva detected
zika,Aedes polynesiensis,mosquito,mixed,lab_experiment,no,,French Polynesia,Richard et al. 2016,weak competence; no infectious saliva at sampled time points
zika,Aedes albopictus,mosquito,mixed,narrative_review,mixed,,Europe,Schulz and Becker 2018,temperature-dependent competence; 18 C not competent, 27 C competent
zika,Culex pipiens,mosquito,not_competent,narrative_review,no,,Italy,Schulz and Becker 2018,not competent at tested temperatures
zika,Culex pipiens biotype pipiens,mosquito,not_competent,narrative_review,no,,Germany,Schulz and Becker 2018,not competent at tested temperatures
zika,Culex pipiens biotype molestus,mosquito,not_competent,narrative_review,no,,Germany,Schulz and Becker 2018,not competent at tested temperatures
zika,Culex torrentium,mosquito,not_competent,narrative_review,no,,Germany,Schulz and Becker 2018,not competent at tested temperatures
```
