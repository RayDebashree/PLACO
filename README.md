
### Description
PLACO implements a variant-level formal statistical test of pleiotropy of two traits using summary-level GWAS data, and can account for potential correlation across traits, such as that arising due to shared controls in case-control studies. The R function `placo` implements this pleiotropic association test. For details of this statistical method, please refer/cite:

Ray, D., Chatterjee, N. "A powerful method for pleiotropic analysis under composite null hypothesis identifies novel shared loci between Type 2 Diabetes and Prostate Cancer". *PLoS Genetics* 16(12): e1009218, https://doi.org/10.1371/journal.pgen.1009218

**Key Words:** Composite null hypothesis; GWAS summary statistics; Intersection-union test; Meta-analysis; Multiple traits; Overlapping samples; Pleiotropy

### Requirements
R (>= 3.0.1)


### How to Install within R
```{r}
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")
```
It is recommended to download/copy the stand-alone R program in this repository, save it in your local directory of choice and `source()` it from your local directory. When a new version of the software is available, older versions may be removed from this repository, and the above `devtools::source_url()` technique may not work.


### Changes
Version 0.1.1 - August 30, 2020
> First public release of the software.


### Notes
1. PLACO and its software is designed to test pleiotropic association of two traits (categorical and/or continuous) from a single study or from two studies. It only requires single-trait GWAS summary statistics. 
    * Caution: PLACO is a single-variant association test; so not expected to work well for rare variants (i.e., genetic variants with very low allele-frequencies).

2. PLACO uses the summary statistics for all variants genome-wide to estimate correlation matrix of the traits. If two studies have overlapping samples/individuals (which may or may not be known), the estimated correlation matrix reflects this overlap. After decorrelating the Z-scores using this correlation matrix, PLACO may be applied.
    * Caution: PLACO does not work well if the two traits are strongly correlated. We advocate using PLACO for two uncorrelated or moderately correlated traits.
    * Caution: PLACO requires genome-wide summary data to infer pleiotropic association of each variant, and cannot be used when summary data on only a handful of candidate genetic variants are available.
    * Caution: PLACO cannot be applied in a hypothesis driven fashion (e.g., test pleiotropy of only a set of variants known to be significantly associated with one trait).

3. Since PLACO uses only summary statistics, it is assumed that all necessary covariate/confounder adjustments were performed when the single-trait summary statistics were obtained.
    * Caution: Harmonize the same effect allele across the two studies/traits so that Z-scores from the two datasets can be jointly analyzed appropriately using PLACO.
    * Recommendation: Remove variants with ![equation](http://www.sciweavers.org/tex2img.php?eq=%20Z%5E%7B2%7D%3E80%20&bc=White&fc=Black&im=png&fs=12&ff=arev&edit=0) (squared Z-scores above 80), similar to recommendation for LD-score regression techniques, since extremely large effect sizes can influence PLACO to show spurious signals.

4. PLACO does not require unrelatedness of samples. When samples are related (e.g., in family-based GWAS), PLACO can use the summary statistics from [EMMAX](https://genome.sph.umich.edu/wiki/EMMAX) (or other univariate mixed model framework) to appropriately test for genetic associations.

5. PLACO does not assume homogeneity of genetic effects of the two traits. 

6. PLACO can only detect statistical association of a variant with two traits, and cannot distinguish between the various types of pleiotropy such as biological or horizontal or vertical/mediated.

7. If you receive an error message like `the integral is probably divergent`, try reducing the absolute tolerance parameter `AbsTol`.
