## Van den Hurk et al. 2010

```csv
Chikungunya virus,Aedes aegypti,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"efficient laboratory vector; transmission 64%"
Chikungunya virus,Ae. albopictus,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"efficient laboratory vector; transmission 32%"
Chikungunya virus,Ae. notoscriptus,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"midgut escape barrier; transmission 20%"
Chikungunya virus,Ae. procax,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"highly susceptible; transmission 64%"
Chikungunya virus,Ae. vigilax,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"highest transmission in study; 76%"
Chikungunya virus,Coquillettidia linealis,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"transmission 75%; efficient laboratory vector"
Chikungunya virus,Culex annulirostris,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"poor laboratory vector but transmission observed; 12%"
Chikungunya virus,Cx. quinquefasciatus,mosquito,not_competent,lab_experiment,no,,"Queensland, Australia",van den Hurk et al. 2010,"refractory to infection"
Chikungunya virus,Cx. sitiens,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"poor laboratory vector; low infection/dissemination/transmission; 4%"
Chikungunya virus,Verrallina funerea,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"transmission 12%; dissemination 28%"
```

Extraction notes: Narrow primary laboratory study of 10 Australian mosquito taxa. This is close to a complete competence read for the species tested. I kept the Culex taxa explicit, including the refractory Cx. quinquefasciatus and the low-but-positive Cx. annulirostris and Cx. sitiens.

## Ledermann et al. 2014

```csv
Chikungunya virus,Aedes hensilli,mosquito,unclear,field_plus_lab,,no,"Yap Island, Federated States of Micronesia",Ledermann et al. 2014,"field surveys were negative; lab infection 62% and dissemination 80%, but no transmission assay"
```

Extraction notes: Narrow outbreak-focused field-plus-lab study. The field material was negative for virus, so I did not label natural infection as present. Lab susceptibility and dissemination were shown, but transmission was not measured, so the row stays unclear.

## Delrieu et al. 2023

```csv
Chikungunya virus,Aedes aegypti,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,"temperature-dependent competence; most studies found higher temperatures increased infection/dissemination/transmission, but some differed"
Chikungunya virus,Aedes albopictus,mosquito,mixed,systematic_review,mixed,,global,Delrieu et al. 2023,"temperature-dependent competence; most studies found higher temperatures increased infection/dissemination/transmission, but two studies found lower competence at high temperature"
```

Extraction notes: Broad systematic review of experimental temperature studies. I treated the source-level result as mixed because the included studies were not uniform, even though the overall trend was toward higher temperature increasing transmission.

## Coffey et al. 2014

```csv
Chikungunya virus,Aedes furcifer,mosquito,competent,narrative_review,yes,,Africa,Coffey et al. 2014,"sylvatic vector; transmission shown in the reviewed studies"
Chikungunya virus,Aedes fulgens,mosquito,competent,narrative_review,yes,,"South Africa",Coffey et al. 2014,"review table reports transmission"
Chikungunya virus,Aedes hensilli,mosquito,unclear,narrative_review,,,"Micronesia",Coffey et al. 2014,"probable vector in Yap; lab susceptibility and dissemination cited, but no transmission assay in the cited study"
Chikungunya virus,Aedes vittatus,mosquito,competent,narrative_review,yes,,"Senegal",Coffey et al. 2014,"review table reports experimental transmission"
Chikungunya virus,Eretmapodites chrysogaster,mosquito,competent,narrative_review,yes,,"not stated",Coffey et al. 2014,"review table reports transmission"
Chikungunya virus,Opifex fuscus,mosquito,competent,narrative_review,yes,,"New Zealand",Coffey et al. 2014,"highly competent at transmitting CHIKV from India"
Chikungunya virus,Culex quinquefasciatus,mosquito,not_competent,narrative_review,no,,global,Coffey et al. 2014,"poor/refractory in the reviewed studies"
Chikungunya virus,Culex pipiens,mosquito,not_competent,narrative_review,no,,"France",Coffey et al. 2014,"review table reports 0 infection and 0 transmission"
Chikungunya virus,Ornithodoros savignyi,tick,not_competent,narrative_review,no,,"South Africa",Coffey et al. 2014,"review table reports 0 infection and no transmission"
```

Extraction notes: Broad narrative review with a large compiled table. I extracted a conservative subset of the explicit species-level claims: the clearest positive vectors, the hensilli uncertainty, and the explicit negative taxa. This is useful support, but not a complete row-by-row reissue of the whole table.

## Higgs and Vanlandingham 2015

```csv
Chikungunya virus,Aedes aegypti,mosquito,competent,narrative_review,yes,,global,Higgs and Vanlandingham 2015,"primary vector; review states A. aegypti and A. albopictus are the urban CHIKV vectors"
Chikungunya virus,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Higgs and Vanlandingham 2015,"major vector; E1-A226V adaptation highlighted"
Chikungunya virus,Aedes hensilli,mosquito,unclear,narrative_review,,,"Yap Island, Federated States of Micronesia",Higgs and Vanlandingham 2015,"probable outbreak vector; lab susceptibility discussed, but transmission not directly shown in the source text"
Chikungunya virus,Opifex fuscus,mosquito,competent,narrative_review,yes,,"New Zealand",Higgs and Vanlandingham 2015,"review says this species is highly competent at transmitting CHIKV from India"
```

Extraction notes: Broad narrative review. I kept the central urban-vector statements and one clearly highlighted alternate vector. The source also reviews much of the same primary literature already captured elsewhere in this file, so I avoided re-listing the full table row set here.

## Campbell et al. 2015

```csv
Chikungunya virus,Aedes albopictus,mosquito,competent,narrative_review,yes,,global,Campbell et al. 2015,"distribution paper with a brief competence statement: chikungunya is readily transmitted by Ae. albopictus, at least in some cases"
```

Extraction notes: Broad distribution-modeling paper, not a competence study. I only kept the one explicit chikungunya competence statement about Ae. albopictus; the rest of the paper is about climate-driven distribution, not transmission experiments.

## ALL TOGETHER

```csv
Chikungunya virus,Aedes aegypti,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"efficient laboratory vector; transmission 64%"
Chikungunya virus,Ae. albopictus,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"efficient laboratory vector; transmission 32%"
Chikungunya virus,Ae. notoscriptus,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"midgut escape barrier; transmission 20%"
Chikungunya virus,Ae. procax,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"highly susceptible; transmission 64%"
Chikungunya virus,Ae. vigilax,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"highest transmission in study; 76%"
Chikungunya virus,Coquillettidia linealis,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"transmission 75%; efficient laboratory vector"
Chikungunya virus,Culex annulirostris,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"poor laboratory vector but transmission observed; 12%"
Chikungunya virus,Cx. quinquefasciatus,mosquito,not_competent,lab_experiment,no,,"Queensland, Australia",van den Hurk et al. 2010,"refractory to infection"
Chikungunya virus,Cx. sitiens,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"poor laboratory vector; low infection/dissemination/transmission; 4%"
Chikungunya virus,Verrallina funerea,mosquito,competent,lab_experiment,yes,,"Queensland, Australia",van den Hurk et al. 2010,"transmission 12%; dissemination 28%"
Chikungunya virus,Aedes hensilli,mosquito,unclear,field_plus_lab,,no,"Yap Island, Federated States of Micronesia",Ledermann et al. 2014,"field surveys were negative; lab infection 62% and dissemination 80%, but no transmission assay"
Chikungunya virus,Aedes furcifer,mosquito,competent,narrative_review,yes,,Africa,Coffey et al. 2014,"sylvatic vector; transmission shown in the reviewed studies"
Chikungunya virus,Aedes fulgens,mosquito,competent,narrative_review,yes,,"South Africa",Coffey et al. 2014,"review table reports transmission"
Chikungunya virus,Aedes vittatus,mosquito,competent,narrative_review,yes,,"Senegal",Coffey et al. 2014,"review table reports experimental transmission"
Chikungunya virus,Eretmapodites chrysogaster,mosquito,competent,narrative_review,yes,,"not stated",Coffey et al. 2014,"review table reports transmission"
Chikungunya virus,Opifex fuscus,mosquito,competent,narrative_review,yes,,"New Zealand",Coffey et al. 2014,"highly competent at transmitting CHIKV from India"
Chikungunya virus,Culex pipiens,mosquito,not_competent,narrative_review,no,,"France",Coffey et al. 2014,"review table reports 0 infection and 0 transmission"
Chikungunya virus,Ornithodoros savignyi,tick,not_competent,narrative_review,no,,"South Africa",Coffey et al. 2014,"review table reports 0 infection and no transmission"
```
