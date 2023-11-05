# Rumen microbiome predictions
Repository for scripts utilized in the manuscript

## Structure of analyses

- *Feature engineering*:
- All taxonomy levels (phylum, class, order, family, and genus) in raw counts, relative abundance, and centered-log ratio (CLR) were used for the model
- 22 alpha-diversity indexes from the default of the microbiome package were also used in the model

- *Ensemble method*:
- Method performed to identify microbes significant across several differential abundance analysis (DAA) tests, including linear regression (LR) and machine learning (ML)
- Results are a summary of variables significant in:
- ALDEx2
- MaAslin2
- ANCOM-BC
- LinDA
- ML as raw microbial counts
- ML as relative abundance
- ML as CLR
- LR as raw microbial counts
- LR as relative abundance
- LR as CLR

## More info can be found by directly addressing the first and corresponding authors

- Hugo Monteiro: hmonteiro@ucdavis.edu
- Fabio Lima: falima@ucdavis.edu
