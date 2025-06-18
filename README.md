
### Description
PLACO+ implements a variant-level formal statistical test of pleiotropy of two traits using summary-level GWAS data, and can account for potential correlation across traits, such as that arising due to measurement of traits on the same set of individuals or due to shared controls in case-control studies. The R function `placo.plus()` implements this pleiotropic association test. PLACO+ may also be used on summary-level data from family-based studies such as trios. For details of this statistical method, please refer/cite:

1. Park, J., Ray, D. (2025+) "A robust pleiotropy method with applications to lipid traits and to inflammatory bowel disease subtypes with sample overlap". *In revision*

PLACO, originally proposed in 2020 and meant for two independent/uncorrelated traits, is a special case of PLACO+. The R function `placo()` is recommended for uncorrelated traits since `placo.plus` can take longer time. For details of this method, please refer/cite:

2. Ray, D., Chatterjee, N. (2020) "A powerful method for pleiotropic analysis under composite null hypothesis identifies novel shared loci between Type 2 Diabetes and Prostate Cancer". *PLoS Genetics* 16(12): e1009218, https://doi.org/10.1371/journal.pgen.1009218

3. Ray, D. et al. (2021) "Pleiotropy method reveals genetic overlap between orofacial clefts at multiple novel loci from GWAS of multi-ethnic trios". *PLoS Genetics* 17(7): e1009584, https://doi.org/10.1371/journal.pgen.1009584

**Key Words:** Composite null hypothesis; GWAS summary statistics; Intersection-union test; Meta-analysis; Multiple traits; Overlapping samples; Pleiotropy

### Requirements
R (>= 3.0.1)


### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE")
```
It is recommended to download/copy the stand-alone R program in this repository, save it in your local directory of choice and `source()` it from your local directory. When a new version of the software is available, older versions may be removed from this repository, and the above `devtools::source_url()` technique may not work.


### Changes
Version 0.2.0 - June 18, 2025
> A more general pleiotropy test PLACO+ is incorporated in this release.

Version 0.1.1 - August 30, 2020
> First public release of the software.


### Notes
1. PLACO+ and its software is designed to test pleiotropic association of two traits (categorical and/or continuous) from a single study or from two studies. It only requires single-trait GWAS summary statistics. 
    * Caution: PLACO+ is a single-variant association test; so not expected to work well for rare variants (i.e., genetic variants with very low allele-frequencies).

2. PLACO+ uses the summary statistics for all variants genome-wide to estimate correlation of the traits. If two studies have overlapping samples/individuals (which may or may not be known), the estimated correlation matrix reflects this overlap. 
    * Caution: PLACO+ is a general version of PLACO originally proposed in 2020. While PLACO+ is meant for any two traits, PLACO is specifically meant for two uncorrelated traits. If two traits are uncorrelated, PLACO+ and PLACO are identical but we recommend using PLACO for uncorrelated traits to save on computation time.
    * Caution: PLACO+ requires genome-wide summary data to infer pleiotropic association of each variant, and cannot be used when summary data on only a handful of candidate genetic variants are available.
    * Caution: PLACO+ cannot be applied in a hypothesis driven fashion (e.g., test pleiotropy of only a set of variants known to be significantly associated with one trait).

3. Since PLACO+ uses only summary statistics, it is assumed that all necessary covariate/confounder adjustments were performed when the single-trait summary statistics were obtained.
    * Caution: Harmonize the same effect allele across the two studies/traits so that Z-scores from the two datasets can be jointly analyzed appropriately using PLACO+.
    * Recommendation: Remove variants with ![equation](https://latex.codecogs.com/svg.image?Z^2%3E80) (squared Z-scores above 80) for any or both traits, similar to recommendation for LD-score regression techniques, since extremely large effect sizes can influence PLACO+ to show spurious signals.

4. PLACO+ is particularly useful for traits with small or modest sample sizes where the identification of pleiotropy by leveraging information from a correlated trait can lead to the identification of novel genetic associations. While PLACO+ is robust to moderate skewness in sample sizes for the two traits, it may show spurious pleiotropy in the presence of heavily skewed sample sizes.
    * Caution: For both traits with massive sample sizes, PLACO+ can have reduced power due to potential removal of pleiotropic variants with large magnitude (due to ![equation](https://latex.codecogs.com/svg.image?Z^2%3E80) filtering). For pleiotropic variants with large effect sizes, a naive approach of declaring pleiotropy based on genome-wide significance of both trait-specific GWAS p-values may be a better choice. However, pleiotropic variants with modest effect sizes will still be missed by such a naive approach. A practical yet powerful alternative in such scenarios might be a 2-step approach: in the first step, PLACO+ is applied as recommended on variants with ![equation](https://latex.codecogs.com/svg.image?Z^2\leq&space;80) for both traits; in the second step, variants removed from PLACO+ analysis can be further screened using the naive approach of declaring pleiotropy based on genome-wide significance of the trait-specific p-values. 

5. PLACO+ does not require unrelatedness of samples. When samples are related, PLACO+ can use the summary statistics from [EMMAX](https://genome.sph.umich.edu/wiki/EMMAX) (or other univariate mixed model framework) to appropriately test for genetic associations.

6. PLACO+ does not assume homogeneity of genetic effects of the two traits. 

7. PLACO+ can only detect statistical association of a variant with two traits, and cannot distinguish between the various types of pleiotropy such as biological or horizontal or vertical/mediated.

8. If you receive an error message like `the integral is probably divergent`, try reducing the absolute tolerance parameter `AbsTol`.

9. For more details on using this R program, please refer to the Supplementary Methods of Park and Ray (2025+). 
