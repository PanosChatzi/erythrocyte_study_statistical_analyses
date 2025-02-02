## Erythrocyte glycolytic and redox metabolism affects muscle oxygenation and exercise performance: a randomized double-blind crossover study in humans

This repository corresponds to our paper, titled "**Erythrocyte glycolytic and redox metabolism affects muscle oxygenation and exercise performance: a randomized double-blind crossover study in humans**", where we investigated the role of erythrocyte redox metabolism in exercise fatigue using *in vivo*, *ex vivo* and computational analyses.

### Contact

-   Panagiotis N. Chatzinikolaou
-   Twitter handle: <https://twitter.com/PanosChatz1>
-   Personal website: [panos-chatz.netlify.app](https://panos-chatz.netlify.app/)
-   Email: [chatzinpn\@phed-sr.auth.gr](mailto:chatzinpn@phed-sr.auth.gr)
-   ORCID ID: <https://orcid.org/0000-0002-8136-1638>

### Repository structure
- **analysis_docs/**
  - Contains Quarto Markdown (`.qmd`) files documenting:
    - `01_Data_Preparation.qmd`: Data cleaning, processing and reshaping from wide to long format for analysis.
    - `02_Figures.qmd`: Code for visualizing the data in figures and panels.
    - `03_Statistics.qmd`: Statistical analyses (repeated measures anova, post-hocs, parametric and non-parametric t-tests, and effect sizes).

- **data/**
  - Includes raw and processed data files:
    - `complete_database.csv`: Dataset in wide format in CSV.
    - `tidyData.RData`: Data in tidy/long format.
  - `README.md`: Documentation for the data files.

- **quantitative_analysis/**
  - Scripts and results for quantitative analysis:
    - `ODC plots.R`: R script for calculating p50 values and plotting Oxygen Dissociation Curves.
    - `p50_analysis.R`: Statistical analysis of p50 values.

- **r_docs/**
  - Custom helper functions for statistics:
    - `MyCohensEffSizes.R`: Functions for calculating Cohen's effect sizes.
    - `MyStatsFunctions.R`: Summary statistics functions.
  - `README.md`: Documentation for the R scripts.


### License

This project uses a [CC-BY 4.0](http://creativecommons.org/licenses/by/4.0/) license. All code is additionally licensed under an [MIT license](https://github.com/PanosChatzi/Erythrocyte-Metabolism/blob/main/LICENSE).
