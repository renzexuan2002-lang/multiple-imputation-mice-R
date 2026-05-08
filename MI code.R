# ============================================================
# Multiple Imputation Pipeline (MICE)
# ESS Hungary Round 11
# ============================================================


# 0. Packages -----------------------------------------------
pkgs <- c("haven", "mice", "VIM", "naniar")

for (p in pkgs) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}


# 1. Load data ----------------------------------------------
data_raw <- read_sav("ESS country data - HU.sav")
ess <- subset(data_raw, essround == 11)

vars <- c("agea","gndr","edlvdahu","hinctnta",
          "wkhtot","wrkctra","wkhct",
          "emplrel","stfedu","health")

ess <- ess[, vars]


# 2. Codebook -----------------------------------------------
lapply(ess, function(x) sort(unique(x)))  #check the values labels

codebook <- data.frame(
  variable = c("agea","gndr","edlvdahu","hinctnta",
               "wkhtot","wrkctra","wkhct","emplrel","stfedu","health"),
  
  label = c(
    "Age of respondent (calculated)",
    "Gender",
    "Highest level of education – Hungary (ISCED)",
    "Household total net income, all sources (decile)",
    "Total hours normally worked per week (incl. overtime)",
    "Employment contract type",
    "Hours per week in main job, contracted",
    "Employment relation",
    "State of education in country nowadays",
    "Subjective general health"
  ),
  
  na_codes = c("999","9","5555,7777,8888,9999","77,88,99",
               "666,777,888,999",
               "6,7,8,9",
               "555,666,777,888,999",
               "6,7,8,9",
               "77,88,99",
               "7,8,9")
)

print(codebook, right = FALSE)
write.csv(codebook, "ESS_HU_codebook.csv", row.names = FALSE)


# 3. Recode missing values ----------------------------------
na_codes <- list(
  hinctnta = c(77,88,99),
  wkhtot   = c(666,777,888,999),
  stfedu   = c(77,88,99),
  wrkctra  = c(6,7,8,9),
  wkhct    = c(555,666,777,888,999),
  emplrel  = c(6,7,8,9),
  health   = c(7,8,9),
  gndr     = c(9),
  agea     = c(999),
  edlvdahu = c(5555,7777,8888,9999)
)

for (v in names(na_codes)) {
  ess[[v]][ess[[v]] %in% na_codes[[v]]] <- NA
}


# 4. Missing data overview ----------------------------------
print(round(colMeans(is.na(ess)) * 100, 1))

VIM::aggr(ess, col = c("steelblue","tomato"),
          numbers = TRUE, sortVars = TRUE, cex.axis = 0.7)

mice::md.pattern(ess, rotate.names = TRUE)

ess_num <- as.data.frame(lapply(ess, function(x)
  if (inherits(x, "haven_labelled")) as.numeric(x) else x))

print(naniar::mcar_test(ess_num))


# 5. Convert labelled variables -----------------------------
ess <- data.frame(lapply(ess, function(x)
  if (inherits(x, "haven_labelled")) as.numeric(x) else x))


# 6. Feature engineering ------------------------------------
ess$gender <- factor(ess$gndr, c(1,2), c("Male","Female"))

ess$edu4 <- cut(ess$edlvdahu,
                breaks = c(-Inf,3,6,9,Inf),
                labels = c("Lower","Lower_secondary","Upper_secondary","Tertiary"),
                ordered_result = TRUE)

ess$health_r <- factor(6 - ess$health, levels = 1:5,
                       labels = c("Very_bad","Bad","Fair","Good","Very_good"),
                       ordered = TRUE)

ess$perm_contract <- factor(ess$wrkctra, c(1,2,3),
                            c("Unlimited","Limited","No_contract"))

ess$emplrel_f <- factor(ess$emplrel, c(1,2,3),
                        c("Employee","Self_employed","Work for Family_business"))


# 7. Imputation dataset -------------------------------------
ess_imp <- ess[, c("agea","gender","edu4","hinctnta","wkhtot",
                   "wkhct","perm_contract","emplrel_f","stfedu","health_r")]


# 8. MICE setup ---------------------------------------------
ini  <- mice(ess_imp, maxit = 0, print = FALSE)
ini
meth <- ini$method
meth
meth[c("agea","hinctnta","wkhtot","wkhct","stfedu")] <- "pmm"
meth["edu4"]          <- "polr"
meth["health_r"]      <- "polr"
meth["gender"]        <- "logreg"
meth["perm_contract"] <- "polyreg"
meth["emplrel_f"]     <- "polyreg"


# 9. Run imputation -----------------------------------------
set.seed(123)
imp <- mice(ess_imp, m = 5, maxit = 20, method = meth, printFlag = FALSE)


# 10. Diagnostics -------------------------------------------
plot(imp)
densityplot(imp, ~ hinctnta + wkhtot + wkhct + stfedu)
stripplot(imp, hinctnta + wkhtot + wkhct + stfedu ~ .imp, pch = 20, cex = 0.7)


# 11. Stability check ---------------------------------------
imp_means <- sapply(1:5, function(i) {
  d <- complete(imp, i)
  d <- as.data.frame(lapply(d, function(x) if (is.factor(x)) as.integer(x) else x))
  colMeans(d)
})

print(round(imp_means, 3))
print(round(apply(imp_means, 1, sd), 4))


# 12. Export ------------------------------------------------
saveRDS(imp, "ESS_HU_mids.rds")

imp_long <- complete(imp, "long", include = TRUE)
write.csv(imp_long, "ESS_HU_imputed_long.csv", row.names = FALSE)

saveRDS(complete(imp, 1), "ESS_HU_imputed_complete.rds")


# 13. Regression: MI pooled vs complete-case ----------------
fit_mi <- with(imp, lm(hinctnta ~ edu4 + health_r + agea))
print(summary(pool(fit_mi)))

cc <- na.omit(ess_imp)
fit_cc <- lm(hinctnta ~ edu4 + health_r + agea, data = cc)
print(summary(fit_cc))
