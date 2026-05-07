# Multiple Imputation Pipeline – ESS Hungary Round 11

This project implements a full **Multiple Imputation by Chained Equations (MICE)** pipeline on survey data from the European Social Survey (ESS) Round 11, Hungary. The goal is to handle missing data properly before running regression analyses, and to compare results between MI-pooled estimates and complete-case analysis.

---

## 📦 Required Packages

```r
install.packages(c("haven", "mice", "VIM", "naniar"))
```

| Package | Purpose |
|---------|---------|
| `haven` | Load SPSS `.sav` files |
| `mice` | Multiple imputation (MICE algorithm) |
| `VIM` | Missing data visualisation |
| `naniar` | MCAR test and missing data summaries |

---

## 📁 Data

**Source:** [European Social Survey (ESS) Round 11 – Hungary](https://ess.sikt.no/en/country/321b06ad-1b98-4b7d-93ad-ca8a24e8788a/hu/)

**Sample:** Hungarian respondents from ESS Round 11

**Variables used:**

| Variable | Description |
|----------|-------------|
| `agea` | Age of respondent (calculated) |
| `gndr` | Gender |
| `edlvdahu` | Highest level of education – Hungary (ISCED) |
| `hinctnta` | Household total net income, all sources (decile) |
| `wkhtot` | Total hours normally worked per week (incl. overtime) |
| `wrkctra` | Employment contract type |
| `wkhct` | Hours per week in main job, contracted |
| `emplrel` | Employment relation |
| `stfedu` | State of education in country nowadays |
| `health` | Subjective general health |

---

## 🔍 Missing Data Analysis

### Missing rates per variable

```
agea     gndr  edlvdahu  hinctnta    wkhtot   wrkctra     wkhct   emplrel    stfedu    health
 1.7      0.0       0.1      23.1      12.0      12.6       9.2       7.5       7.4       0.1
```

`hinctnta` (household income) has the highest missingness at **23.1%**, followed by employment-related variables.

### Missing data pattern

![Missing data overview](output/figures/01_aggr_missing.png)

The left panel shows the proportion of missing values per variable. The right panel shows the combinations of missingness across variables — income and employment variables tend to be missing together.

### Little's MCAR Test

```
statistic = 2473, df = 330, p-value < 0.001, missing patterns = 46
```

The test strongly rejects the null hypothesis that data is Missing Completely At Random (MCAR). This confirms that **multiple imputation is more appropriate than Complete-case analysis**.

---

## ⚙️ Imputation Setup

The MICE algorithm was run with **m = 5 imputations** and **20 iterations**, using variable-appropriate methods:

| Variable | Type | Method |
|----------|------|--------|
| `agea`, `hinctnta`, `wkhtot`, `wkhct`, `stfedu` | Continuous | Predictive Mean Matching (`pmm`) |
| `edu4`, `health_r` | Ordered categorical | Proportional odds model (`polr`) |
| `gender` | Binary | Logistic regression (`logreg`) |
| `perm_contract`, `emplrel_f` | Nominal categorical | Polytomous regression (`polyreg`) |

---

## 📊 Diagnostics

### Convergence plots

The convergence plots show that the mean and standard deviation of imputed values mix well across the 5 chains over 20 iterations, indicating the algorithm has converged.

![Convergence – agea, edu4, hinctnta](output/figures/02_convergence_1.png)

![Convergence – wkhtot, wkhct, perm_contract](output/figures/03_convergence_2.png)

![Convergence – emplrel_f, stfedu, health_r](output/figures/04_convergence_3.png)

### Density plots

Blue lines show the observed data distribution; red lines show the 5 imputed datasets. Close overlap indicates plausible imputations.

![Density plots](output/figures/05_density.png)

`hinctnta` and `stfedu` show good overlap. `wkhct` contains some extreme values (contracted hours reported as 555+) which are visible as outliers but preserved in the imputation.

### Strip plots

![Strip plots](output/figures/06_stripplot.png)

Imputed values (red, imputation numbers 2–6) follow a similar distribution to observed values (blue, imputation 1), confirming imputation quality.

### Stability across imputations

```
             [,1]   [,2]   [,3]   [,4]   [,5]   SD
agea       50.455 50.531 50.596 50.576 50.683  0.084
hinctnta    5.526  5.521  5.512  5.541  5.554  0.017
wkhtot     43.687 43.705 43.550 43.472 43.378  0.140
stfedu      4.139  4.144  4.145  4.109  4.161  0.019
```

Standard deviations across the 5 imputed datasets are very small, confirming stable and consistent imputations.

---

## 📈 Results: MI vs Complete-Case Regression

**Outcome:** Household income decile (`hinctnta`)  
**Predictors:** Education level (`edu4`), subjective health (`health_r`), age (`agea`)

### MI Pooled estimates (Rubin's Rules)

```
term           estimate   std.error   statistic      p.value
(Intercept)    6.042      0.279       21.66          < 0.001
edu4.L         2.354      0.157       14.97          < 0.001
edu4.Q        -0.253      0.154       -1.64           0.114
edu4.C         0.366      0.126        2.91           0.007
health_r.L     1.848      0.308        6.00          < 0.001
health_r.C    -0.542      0.203       -2.67           0.008
agea          -0.018      0.004       -4.23          < 0.001
```

### Complete-Case estimates (listwise deletion)

```
term           estimate   std.error   t value    p.value
(Intercept)    5.844      0.308       18.997     < 0.001
edu4.L         2.335      0.178       13.120     < 0.001
health_r.L     2.056      0.384        5.352     < 0.001
agea          -0.014      0.005       -2.864      0.004

R-squared: 0.261 (1,273 obs after deletion)
```

### Key findings

- **Education** is the strongest predictor of income in both models (linear trend: β ≈ 2.35 in MI vs 2.34 in CC).
- **Health** has a significant positive linear effect on income, but MI yields a smaller estimate (β = 1.85) compared to complete-case (β = 2.06) — suggesting complete-case analysis may overestimate this effect.
- **Age** shows a small negative effect on income. MI gives a stronger estimate (β = −0.018) than complete-case (β = −0.014), likely because listwise deletion removes older, lower-income respondents disproportionately.
- Multiple imputation recovers information from **~23% missing income data**, producing more reliable and generalisable estimates.

---

## 📤 Output Files

| File | Description |
|------|-------------|
| `ESS_HU_mids.rds` | Full `mids` object (all 5 imputations, for further analysis in R) |
| `ESS_HU_imputed_long.csv` | Long-format dataset with all 5 imputed datasets |
| `ESS_HU_imputed_complete.rds` | Single completed dataset (imputation 1) |
| `ESS_HU_codebook.csv` | Variable codebook |

---

## ▶️ How to Run

1. Place `ESS country data - HU.sav` in your working directory
2. Open `MI_pipeline.R` in RStudio
3. Run the script top to bottom (install packages if prompted)
4. Output files will be saved in your working directory

---

## 📚 References

- van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software*, 45(3), 1–67.
- European Social Survey (2024). ESS Round 11 – Hungary. Norwegian Social Science Data Services.
- Little, R. J. A. (1988). A test of Missing Completely at Random for multivariate data with missing values. *Journal of the American Statistical Association*, 83(404), 1198–1202.
