# Load databases ----
covid19_result       <- loadCovid19Result(ukb)
sequela_table        <- loadSequelaTable()
genetic_data         <- loadGeneticData(ukb)
waves_data           <- loadWaveData(ukb)

# Create LC cohort ----
health_questionnaire <- loadHealthQuestionnaire(ukb)
x <- covid19_result |>
  select("eid","specdate","result") |>
  recordAttrition() |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  inner_join(
    health_questionnaire |>
      filter(!is.na(questionnaire_started)),
    by = "eid"
  ) |>
  recordAttrition("Restrict to participants that answered the health and well-being web-questionnaire") |>
  filter(questionnaire_started != as.Date("1999-01-01")) |>
  recordAttrition("Restrict to participants that answered yes/no to all the questions from the health and well-being web-questionnaire.") |>
  filter(specdate < questionnaire_started) |>
  recordAttrition("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction before answering the questionnaire.")

x1 <- x |> group_by(eid) |>
  filter(specdate == max(specdate)) |>
  ungroup() |>
  restoreAttrition(x) |>
  recordAttrition("Restrict to the latest record.")

x2 <- x1 |>
  filter(((specdate+washout_period) < questionnaire_started) & (specdate > (questionnaire_started-365))) |>
  recordAttrition(paste0("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction between 1 year and ",washout_period," days before answering the questionnaire."))

x3 <- x2 |>
  rowwise() |>
  mutate("symptom" = across(starts_with("symptom")) |> sum()) |>
  ungroup() |>
  restoreAttrition(x2) %>%
  mutate(length = do.call(pmax, c(select(., starts_with("length")), na.rm = TRUE))) |>
  select(-starts_with("length_")) |>
  filter((questionnaire_started-length) > specdate) |>
  recordAttrition("Restrict to people with no persistent symptoms") |>
  inner_join(
    loadGeneticData(ukb), by = "eid"
  ) |>
  filter(!is.na(pc1)) |>
  recordAttrition("Restrict to people with available genetic data") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with equal sex and genetic sex recorded")


longCovid_cases_cohort <- x3 |>
  filter(symptom >= n_symptoms) |>
  recordAttrition("Restrict to people that reported at least one symptom of the questionnaire")

longCovid_controls_cohort <- x3 |>
  filter(symptom < n_symptoms & symptom >= 0) |>
  recordAttrition("Restrict to people that did not report any symptom of the questionnaire")

longCovid_cohort <- longCovid_cases_cohort |>
  select("eid", "specdate", "questionnaire_started", "year_birth") |>
  mutate(state = 1) |>
  union_all(
    longCovid_controls_cohort |>
      select("eid", "specdate", "questionnaire_started", "year_birth") |>
      mutate(state = 0)
  ) 

longCovid_cohort <- longCovid_cohort |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  select(-c("year_of_birth")) |>
  filter(!is.na(pc1)) |>
  recordAttrition("Restrict to people with available genetic data") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with equal sex and genetic sex recorded") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate")


longCovid_attrition <- attr(longCovid_cohort, "cohort_attrition")

longCovid_cohort <- as.phe(longCovid_cohort |> 
                             select("IID" = "eid", "state", "age" = "age_when_infected", "sex", "age2", "agesex", "batch", starts_with("pc")) |> 
                             mutate("FID" = IID) |> 
                             relocate("FID"), "FID", "IID")

# Create LC cohort (With subtypes)----
health_questionnaire_wide <- loadHealthAndWellBeingQuestionnaireWideSubtypes(ukb)

cohort <- covid19_result %>%
  dplyr::select(eid, specdate, result) %>%
  recordAttrition() %>%
  filter(result == 1) %>%
  dplyr::select(-"result") %>%
  distinct() %>%
  recordAttrition("Restrict to participants with a positive COVID-19 test result") %>%
  inner_join(
    health_questionnaire_wide %>%
      filter(!is.na(questionnaire_started)),
    by = "eid"
  ) %>%
  recordAttrition("Restrict to participants that answered the health and well-being web-questionnaire") |>
  filter(specdate < questionnaire_started) %>%
  recordAttrition("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction before answering the questionnaire.")

cohort_latest <- cohort %>%
  group_by(eid) %>%
  filter(specdate == max(specdate)) %>%
  ungroup() %>%
  restoreAttrition(cohort) %>%
  recordAttrition("Restrict to the latest record.")

cohort_filtered <- cohort_latest %>%
  filter(
    ((specdate + washout_period) < questionnaire_started) &
      (specdate > (questionnaire_started - 365))
  ) %>%
  recordAttrition(paste0("Restrict to participants with a SARS-CoV-2 infection confirmed by a positive polymerase chain test reaction between 1 year and ",washout_period," days before answering the questionnaire."))

# Identify symptom columns
symptom_cols <- gsub("symptom_", "", colnames(cohort_filtered %>% dplyr::select(starts_with("symptom"))))

# Process symptom variables
for (symptom in symptom_cols) {
  var <- paste0("symptom_", symptom)
  var_len <- paste0("length_", symptom)
  
  cohort_filtered <- cohort_filtered %>%
    mutate(!!sym(var) := if_else(!!sym(var) < 0 | is.na(!!sym(var)), 0, !!sym(var))) |> # if people didn't answer the questionnaire for this symptom -> control
    filter(!((questionnaire_started - !!sym(var_len)) < specdate & !is.na(!!sym(var_len)))) |>
    # mutate(!!sym(var) := if_else(
    #   (questionnaire_started - !!sym(var_len)) < specdate & # pre-existing symptom
    #     !is.na(!!sym(var_len)), # and the length is not NA and the symptom is not NA
    #   NA,
    #   !!sym(var))) |>
    mutate(!!sym(var_len) := if_else(!!sym(var) == 0, NA, !!sym(var_len)))
}

cohort_aggregated <- cohort_filtered %>%
  rowwise() %>%
  mutate(symptom = sum(c_across(starts_with("symptom_")), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(length = do.call(pmin, c(dplyr::select(., starts_with("length_")), na.rm = TRUE))) %>%
  dplyr::select(-starts_with("length_")) 

longCovid_cohort_wide_subtypes <- cohort_aggregated %>%
  restoreAttrition(cohort_filtered) %>%
  # filter(((questionnaire_started - length) > specdate) | symptom == 0) |>
  recordAttrition("Restrict to participants with no symptoms or non-persistent symptoms") |>
  mutate(symptom = if_else(symptom > 1, 1, symptom)) |>
  select("eid", "specdate", 
         "questionnaire_started", 
         "ent" = "symptom_ent", 
         "neuro"  = "symptom_neurological_symptoms",
         "cardio" = "symptom_respiratory_chest_symptoms",
         "fatigue" = "symptom_general_systemic_symptoms",
         "overall" = "symptom") 

longCovid_cohort_wide_subtypes |>
  summarise("ent_cases"   = sum(ent),   "ent_controls" = sum(ent == 0),
            "neuro_cases" = sum(neuro), "neuro_controls" = sum(neuro == 0),
            "cardio_cases"  = sum(cardio),  "cardio_controls" = sum(cardio == 0),
            "fatigue_cases" = sum(fatigue), "fatigue_controls" = sum(fatigue == 0),
            "overall_cases" = sum(overall), "overall_controls" = sum(overall == 0)) |>
  tidyr::pivot_longer(everything())


longCovid_cohort_wide_subtypes <- longCovid_cohort_wide_subtypes |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with the same sex and genetic sex recorded") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate")

longCovid_cohort_wide_subtypes_attrition <- attr(longCovid_cohort_wide_subtypes, "cohort_attrition")

gwas_file <- longCovid_cohort_wide_subtypes |>
  mutate(FID = eid,
         IID = eid) |>
  mutate("genetic_ethnic_grouping" = if_else(is.na(genetic_ethnic_grouping), "Other", "Causcasian")) |>
  select("FID", "IID", 
         "ent", "neuro", "cardio", "fatigue", "overall", 
         "G_batch" = "batch", "G_ethnic" = "genetic_ethnic_grouping",
         "G_PCA_1" = "pc1", "G_PCA_2" = "pc2", "G_PCA_3" = "pc3", "G_PCA_4" = "pc4",        
         "G_PCA_5" = "pc5", "G_PCA_6" = "pc6", "G_PCA_7" = "pc7", "G_PCA_8" = "pc8",  
         "G_PCA_9" = "pc9", "G_PCA_10" = "pc10", "G_PCA_11" = "pc11", "G_PCA_12" = "pc12", 
         "G_PCA_13" = "pc13", "G_PCA_14" = "pc14", "G_PCA_15" = "pc15", "G_PCA_16" = "pc16", 
         "G_PCA_17" = "pc17", "G_PCA_18" = "pc18",  "G_PCA_19" = "pc19", "G_PCA_20" = "pc20",
         "BS_age"  = "age_when_infected",
         "BS_age2" = "age2",
         "BS_agesex" = "agesex",
         "BS_sex" = "sex")

write_tsv(x = gwas_file |> select("FID", "IID", starts_with("BS"), starts_with("G_")),
          file = paste0(dir_data,"/Results/covariates.txt"))
write_tsv(x = gwas_file |> select("FID", "IID", "ent", "cardio", "neuro", "fatigue", "overall"),
          file = paste0(dir_data,"/Results/phenotype.txt"))
write_tsv(x = gwas_file |> select("FID", "IID"),
          file = paste0(dir_data,"/Results/sample_included.txt"), col_names = FALSE)

rm(list = c("health_questionaire", 
            "health_questionnaire_wide", "longCovid_cases_wide", "longCovid_controls_wide", 
            "cohort", "cohort_latest", "cohort_filtered", "cohort_aggregated", "symptom_cols", "var", "var_len"))


# Create PACS cohort ----
t <- hes |>
  recordAttrition() |>
  mutate(diag_icd10   = if_else(!diag_icd10 %in% sequela_table$icd10_code, NA, diag_icd10),
         episode_date = if_else(is.na(diag_icd10), NA, episode_date)) |>
  distinct() |>
  recordAttrition("Restrict to PACS records") |>
  inner_join(
    covid19_result |>
      select("eid","specdate","result"),
    by = "eid",
    relationship = "many-to-many"
  ) |>
  recordAttrition("Restrict to participants with COVID-19 linkage data") |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  mutate(episode_date = if_else(episode_date >= (specdate-365), episode_date, NA)) |>
  mutate(episode_date = if_else(episode_date <= (specdate+365), episode_date, NA)) |>
  mutate(diag_icd10 = if_else(is.na(episode_date), NA, diag_icd10)) |>
  distinct() |>
  recordAttrition("Restrict the study period between one year before and after testing positive for COVID-19.") |>
  mutate(exclusion = if_else(
    episode_date >= (specdate-365) & episode_date <= specdate & (!is.na(diag_icd10)),
    1,0
  ))

t1 <- t |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t) |>
  filter(exclusion == 0) |>
  recordAttrition("Restrict to participants that did not have any PACS event one year prior the COVID-19 infection") |>
  mutate(exclusion = if_else((!is.na(diag_icd10)) &
                               (episode_date > specdate) &
                               (episode_date <= (specdate+washout_period)), 1, 0))
t2 <- t1 |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t1) |>
  filter(exclusion == 0) |>
  recordAttrition(paste0("Restrict to participants that did not have any PACS event after ", washout_period, " days of the infection")) |>
  mutate(diagnoses = if_else(!is.na(diag_icd10), 1, 0))

t3 <- t2 |>
  group_by(eid) |>
  mutate(cases = sum(diagnoses)) |>
  ungroup() |>
  restoreAttrition(t2)

pacs_cases_cohort <- t3 |>
  filter(cases >= n_symptoms) |>
  filter(!is.na(diag_icd10)) |>
  recordAttrition("Restrict to cases")

pacs_controls_cohort <- t3 |>
  filter(cases >= 0 & cases < n_symptoms) |>
  distinct() |>
  recordAttrition("Restrict to controls")

pacs_cohort  <- pacs_cases_cohort |>
  union_all(pacs_controls_cohort) |>
  select("eid", "specdate", "state" = "cases") |>
  group_by(eid) |>
  mutate(specdate = min(specdate),
         state = if_else(state < n_symptoms, 0, 1)) |>
  distinct() |>
  ungroup() |>
  restoreAttrition(pacs_cases_cohort) |>
  recordAttrition("Restrict to the earliest record")

rm(list = c("t", "t1", "t2", "t3", "pacs_cases_cohort", "pacs_controls_cohort"))

# Prepare .phe files and save cohorts -----


longCovid_cohort_wide <- longCovid_cohort_wide |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with the same sex and genetic sex recorded") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate")

longCovid_wide_attrition <- attr(longCovid_cohort_wide, "cohort_attrition")

longCovid_cohort_wide <- as.phe(longCovid_cohort_wide |> 
                             select("IID" = "eid", "state", "age" = "age_when_infected", "sex", "age2", "agesex", "batch", starts_with("pc")) |> 
                             mutate("FID" = IID) |> 
                             relocate("FID"), "FID", "IID")

pacs_cohort <- pacs_cohort |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  filter(!is.na(pc1)) |> 
  recordAttrition("Restrict to people with the same sex and genetic sex recorded and all principal components") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate")

pacs_cohort <- as.phe(pacs_cohort |> 
                        select("IID" = "eid", "state", "age" = "age_when_infected", "sex", "age2", "agesex", "batch", starts_with("pc")) |> 
                        mutate("FID" = IID) |> 
                        relocate("FID"), "FID", "IID")

pacs_attrition <- attr(pacs_cohort, "cohort_attrition")

save(longCovid_cohort,
     longCovid_cohort_wide,
     pacs_cohort,
     longCovid_attrition,
     longCovid_wide_attrition,
     pacs_attrition,
     file = paste0(dir_data,"Results/cohortsData.Rdata"))

write.phe(paste0(dir_data,"/Results/LongCovid_cohort.phe"), longCovid_cohort)
write.phe(paste0(dir_data,"/Results/LongCovid_wide_cohort.phe"), longCovid_cohort_wide)
write.phe(paste0(dir_data,"/Results/Pacs_cohort.phe"), pacs_cohort)

rm(list = c("longCovid_cohort", "longCovid_cohort_wide", "pacs_cohort", "longCovid_attrition", "longCovid_wide_attrition", "pacs_attrition"))

# Create LC validation cohorts ----
t1 <- waves_data |>
  recordAttrition() |>
  filter(across(ends_with("_antibody_test_result"), ~!is.na(.))) |>
  recordAttrition("Participants with no missing antibody test data within the first 6 waves") |>
  filter(across(ends_with("self-reported_date_that_antibody_test_sample_was_collected"), ~ !. %in% c(as.Date("1900-01-01"), as.Date("1999-01-01")))) |>
  recordAttrition("Participants with a valid self-reported date that antibody test sample was collected") |>
  mutate("infection" = rowSums(across(ends_with("_antibody_test_result")), na.rm = TRUE)) |> 
  filter(infection > 0) |>
  recordAttrition("Participants with at least one reported COVID-19 infection") 

# Keep only those dates with a positive test result
for(i in 0:5){
  date_col <- paste0("w",i,"_self-reported_date_that_antibody_test_sample_was_collected")
  inf_col  <- paste0("w",i,"_antibody_test_result")
  t1[,date_col][!t1[,inf_col]] <- NA
}

# Keep the earliest infection date
t1 <- t1 |>
  rowwise() |>
  mutate("infection_date" = min(c_across(ends_with("self-reported_date_that_antibody_test_sample_was_collected")), na.rm = TRUE)) |>
  ungroup() |>
  select(-c(contains("antibody_test_result"), contains("self-reported_date_that_antibody_test_sample_was_collected"),
            contains("method_of_returning")))

# Merge all symptoms for each date
for(i in 0:8){
  t1 <- t1 |>
    rowwise() |>
    mutate(!!paste0("symptom_w",i) := max(c_across(starts_with(paste0("w",i,"_covid-19_symptom_"))), na.rm = TRUE)) |>
    ungroup() |>
    select(-starts_with(paste0("w",i,"_covid-19_symptom_"))) |>
    mutate(!!paste0("date_w",i) := .data[[paste0("w",i,"_covid-19_date_questionnaire_results_were_received_by_post")]]-30)
}

t1 <- t1 |> 
  select(-starts_with("w")) |>
  # Replace Inf and -Inf by NA
  mutate(across(everything(), ~gsub("-Inf|Inf", NA,.)))

# Determine LC cases
for(i in 0:8){
  t1 <- t1 |>
    mutate("diff_days" = as.numeric(difftime(.data[[paste0("date_w",i)]], infection_date))) |>
    mutate(!!paste0("case_w",i) := if_else(diff_days > 30 & .data[[paste0("symptom_w",i)]] == 1, 1, 0)) |>
    select(-c(paste0("symptom_w",i), paste0("date_w",i)))
}

t1 <- t1 |>
  rowwise() |>
  mutate("case" = max(c_across(starts_with("case_")), na.rm = TRUE)) |>
  select(-c(starts_with("case_"),"diff_days", "infection"))

t1 <- t1 |>
  mutate(eid = as.numeric(eid)) |>
  addGwasCovariates(ukb) |> 
  mutate("age_when_infected" = year(infection_date) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  filter(is.na(heterozygosity)) |>
  rename("state" = "case")

gwas_file <- t1 |>
  mutate(FID = eid,
         IID = eid) |>
  mutate("genetic_ethnic_grouping" = if_else(is.na(genetic_ethnic_grouping), "Other", "Causcasian")) |>
  select("FID", "IID", 
         "overall" = "state",
         "G_batch" = "batch", "G_ethnic" = "genetic_ethnic_grouping",
         "G_PCA_1" = "pc1", "G_PCA_2" = "pc2", "G_PCA_3" = "pc3", "G_PCA_4" = "pc4",        
         "G_PCA_5" = "pc5", "G_PCA_6" = "pc6", "G_PCA_7" = "pc7", "G_PCA_8" = "pc8",  
         "G_PCA_9" = "pc9", "G_PCA_10" = "pc10", "G_PCA_11" = "pc11", "G_PCA_12" = "pc12", 
         "G_PCA_13" = "pc13", "G_PCA_14" = "pc14", "G_PCA_15" = "pc15", "G_PCA_16" = "pc16", 
         "G_PCA_17" = "pc17", "G_PCA_18" = "pc18",  "G_PCA_19" = "pc19", "G_PCA_20" = "pc20",
         "BS_age"  = "age_when_infected",
         "BS_age2" = "age2",
         "BS_agesex" = "agesex",
         "BS_sex" = "sex")

write_delim(x = gwas_file |> select("FID", "IID", starts_with("BS"), starts_with("G_")),
            file = paste0(dir_data,"/Results/covariates_validation.txt"))
write_delim(x = gwas_file |> select("FID", "IID","overall"),
            file = paste0(dir_data,"/Results/phenotype_validation.txt"))
write_delim(x = gwas_file |> select("FID", "IID"),
            file = paste0(dir_data,"/Results/sample_included_validation.txt"), col_names = FALSE)

t1_cohort <- as.phe(t1 |> 
                      select("IID" = "eid", "state", "age" = "age_when_infected", "sex", "age2", "agesex", "batch", starts_with("pc")) |> 
                      mutate("FID" = IID) |> 
                      relocate("FID"), "FID", "IID")

write.phe(paste0(dir_data,"/Results/LongCovid_validation_cohort.phe"), t1_cohort)

# Create PACS narrower cohorts ----
sequela_arterial_table  <- loadSequelaArterialTable()
sequela_venous_table    <- loadSequelaVenousTable()

t <- hes |>
  recordAttrition() |>
  mutate(diag_icd10   = if_else(!diag_icd10 %in% sequela_arterial_table$icd10_code, NA, diag_icd10),
         episode_date = if_else(is.na(diag_icd10), NA, episode_date)) |>
  distinct() |>
  recordAttrition("Restrict to arterial PACS records") |>
  inner_join(
    covid19_result |>
      select("eid","specdate","result"),
    by = "eid",
    relationship = "many-to-many"
  ) |>
  recordAttrition("Restrict to participants with COVID-19 linkage data") |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  mutate(episode_date = if_else(episode_date >= (specdate-365), episode_date, NA)) |>
  mutate(episode_date = if_else(episode_date <= (specdate+365), episode_date, NA)) |>
  mutate(diag_icd10 = if_else(is.na(episode_date), NA, diag_icd10)) |>
  distinct() |>
  recordAttrition("Restrict the study period between one year before and after testing positive for COVID-19.") |>
  mutate(exclusion = if_else(
    episode_date >= (specdate-365) & episode_date <= specdate & (!is.na(diag_icd10)),
    1,0
  ))

t1 <- t |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t) |>
  filter(exclusion == 0) |>
  recordAttrition("Restrict to participants that did not have any PACS event one year prior the COVID-19 infection") |>
  mutate(exclusion = if_else((!is.na(diag_icd10)) &
                               (episode_date > specdate) &
                               (episode_date <= (specdate+washout_period)), 1, 0))
t2 <- t1 |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t1) |>
  filter(exclusion == 0) |>
  recordAttrition(paste0("Restrict to participants that did not have any PACS event after ", washout_period, " days of the infection")) |>
  mutate(diagnoses = if_else(!is.na(diag_icd10), 1, 0))

t3 <- t2 |>
  group_by(eid) |>
  mutate(cases = sum(diagnoses)) |>
  ungroup() |>
  restoreAttrition(t2)

pacs_arterial_cases_cohort <- t3 |>
  filter(cases >= n_symptoms) |>
  filter(!is.na(diag_icd10)) |>
  recordAttrition("Restrict to cases")

pacs_arterial_controls_cohort <- t3 |>
  filter(cases >= 0 & cases < n_symptoms) |>
  distinct() |>
  recordAttrition("Restrict to controls")

pacs_arterial_cohort  <- pacs_arterial_cases_cohort |>
  union_all(pacs_arterial_controls_cohort) |>
  select("eid", "specdate", "state" = "cases") |>
  group_by(eid) |>
  mutate(specdate = min(specdate),
         state = if_else(state < n_symptoms, 0, 1)) |>
  distinct() |>
  ungroup() |>
  restoreAttrition(pacs_arterial_cases_cohort) |>
  recordAttrition("Restrict to the earliest record")

pacs_arterial_cohort <- pacs_arterial_cohort |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with the same sex and genetic sex recorded") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate")

pacs_arterial_attrition <- attr(pacs_arterial_cohort, "cohort_attrition")

pacs_arterial_cohort <- as.phe(pacs_arterial_cohort |> 
                        select("IID" = "eid", "state", "age" = "age_when_infected", "sex", "age2", "agesex", "batch", starts_with("pc")) |> 
                        mutate("FID" = IID) |> 
                        relocate("FID"), "FID", "IID")

rm(list = c("t", "t1", "t2", "t3", "pacs_arterial_cases_cohort", "pacs_arterial_controls_cohort"))

t <- hes |>
  recordAttrition() |>
  mutate(diag_icd10   = if_else(!diag_icd10 %in% sequela_venous_table$icd10_code, NA, diag_icd10),
         episode_date = if_else(is.na(diag_icd10), NA, episode_date)) |>
  distinct() |>
  recordAttrition("Restrict to venous PACS records") |>
  inner_join(
    covid19_result |>
      select("eid","specdate","result"),
    by = "eid",
    relationship = "many-to-many"
  ) |>
  recordAttrition("Restrict to participants with COVID-19 linkage data") |>
  filter(result == 1) |>
  select(-"result") |>
  distinct() |>
  recordAttrition("Restrict to participants with a positive COVID-19 test result") |>
  mutate(episode_date = if_else(episode_date >= (specdate-365), episode_date, NA)) |>
  mutate(episode_date = if_else(episode_date <= (specdate+365), episode_date, NA)) |>
  mutate(diag_icd10 = if_else(is.na(episode_date), NA, diag_icd10)) |>
  distinct() |>
  recordAttrition("Restrict the study period between one year before and after testing positive for COVID-19.") |>
  mutate(exclusion = if_else(
    episode_date >= (specdate-365) & episode_date <= specdate & (!is.na(diag_icd10)),
    1,0
  ))

t1 <- t |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t) |>
  filter(exclusion == 0) |>
  recordAttrition("Restrict to participants that did not have any PACS event one year prior the COVID-19 infection") |>
  mutate(exclusion = if_else((!is.na(diag_icd10)) &
                               (episode_date > specdate) &
                               (episode_date <= (specdate+washout_period)), 1, 0))
t2 <- t1 |>
  group_by(eid) |>
  mutate(exclusion = max(exclusion)) |>
  ungroup() |>
  restoreAttrition(t1) |>
  filter(exclusion == 0) |>
  recordAttrition(paste0("Restrict to participants that did not have any PACS event after ", washout_period, " days of the infection")) |>
  mutate(diagnoses = if_else(!is.na(diag_icd10), 1, 0))

t3 <- t2 |>
  group_by(eid) |>
  mutate(cases = sum(diagnoses)) |>
  ungroup() |>
  restoreAttrition(t2)

pacs_venous_cases_cohort <- t3 |>
  filter(cases >= n_symptoms) |>
  filter(!is.na(diag_icd10)) |>
  recordAttrition("Restrict to cases")

pacs_venous_controls_cohort <- t3 |>
  filter(cases >= 0 & cases < n_symptoms) |>
  distinct() |>
  recordAttrition("Restrict to controls")

pacs_venous_cohort  <- pacs_venous_cases_cohort |>
  union_all(pacs_venous_controls_cohort) |>
  select("eid", "specdate", "state" = "cases") |>
  group_by(eid) |>
  mutate(specdate = min(specdate),
         state = if_else(state < n_symptoms, 0, 1)) |>
  distinct() |>
  ungroup() |>
  restoreAttrition(pacs_venous_cases_cohort) |>
  recordAttrition("Restrict to the earliest record")

pacs_venous_cohort <- pacs_venous_cohort |>
  addGwasCovariates(ukb) |>
  mutate("age_when_infected" = year(specdate) - year_of_birth) |>
  mutate("age2" = age_when_infected^2) |>
  mutate("agesex" = age_when_infected*sex) |>
  select(-c("year_of_birth")) |>
  filter(genetic_sex == sex) |>
  recordAttrition("Restrict to people with the same sex and genetic sex recorded") |>
  filter(is.na(sex_chromosome_aneuploidy)) |>
  recordAttrition("Restrict to people with no sex chromosome aneuploidy") |>
  filter(is.na(heterozygosity)) |>
  recordAttrition("Restrict to people that are no outliers for heterozygosity or missing rate")

pacs_venous_attrition <- attr(pacs_venous_cohort, "cohort_attrition")

pacs_venous_cohort <- as.phe(pacs_venous_cohort |> 
                                 select("IID" = "eid", "state", "age" = "age_when_infected", "sex", "age2", "agesex", "batch", starts_with("pc")) |> 
                                 mutate("FID" = IID) |> 
                                 relocate("FID"), "FID", "IID")

rm(list = c("t", "t1", "t2", "t3", "pacs_venous_cases_cohort", "pacs_venous_controls_cohort"))

save(pacs_arterial_cohort,
     pacs_venous_cohort,
     pacs_arterial_attrition,
     pacs_venous_attrition,
     file = paste0(dir_data,"Results/cohortsData_extraPacs.Rdata"))

write.phe(paste0(dir_data,"/Results/Pacs_arterial_cohort.phe"), pacs_arterial_cohort)
write.phe(paste0(dir_data,"/Results/Pacs_venous_cohort.phe"), pacs_venous_cohort)

rm(list = c("pacs_arterial_cohort", "pacs_arterial_attrition", "pacs_venous_cohort", "pacs_venous_attrition"))



  
