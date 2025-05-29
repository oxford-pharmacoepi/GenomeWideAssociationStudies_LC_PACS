# Validation table -----
for(name in c("LongCovid", "PACS", "Neuro", "Cardio", "ENT", "Fatigue")){
  path <- paste0(dir_results, "FUMA/",name,"/")
  
  gwas <- read_table(paste0(dir_data,"Results/GWAS/",name,".txt"))
  leadSnps    <- read_delim(paste0(path,'leadSNPs.txt'))
  annovar     <- read_delim(paste0(path,'annov.txt'))
  gwascatalog <- read_delim(paste0(path,'gwascatalog.txt'))
  genomicRL   <- read_delim(paste0(path,'GenomicRiskLoci.txt'))
  
  tab <- leadSnps |>
    select('uniqID',
           'CHR' = 'chr',
           'BP'  = 'pos',
           'SNP' = 'rsID') |>
    left_join(gwas |>
                mutate(OR = exp(BETA)) |>
                mutate(lower = exp(BETA-1.96*SE),
                       higher = exp(BETA+1.96*SE)) |>
                select('CHR' = 'CHROM',
                       'BP'  = 'GENPOS',
                       'SNP' = 'ID',
                       'EA'  = 'ALLELE1',
                       'OA'  = 'ALLELE0',
                       'EAF' = 'A1FREQ',
                       'OR'  = 'OR',
                       "lower", "higher",
                       'SE'  = 'SE',
                       "N"   = "N",
                       'PVAL' = 'PVAL'),
              by = c('CHR','BP','SNP')) |>
    left_join(annovar |>
                group_by(uniqID, 'Function' = annot) |>
                summarise('Gene' = paste0(symbol, collapse=" - "), 
                          'Distance (b)' = paste0(dist, collapse=" - "))|>
                ungroup(),
              by = "uniqID") |>
    left_join(
      gwascatalog |> rename("SNP" = "IndSigSNP") |>
        select('SNP', 'Trait') |>
        distinct() |>
        group_by(SNP) |>
        summarise('Previous reported trait' = paste0(Trait, collapse = ", ")),
      by = 'SNP') |>
    select(-uniqID) |>
    mutate('EAF' = round(EAF, digits=2),
           'OR'  = round(OR,  digits=2),
           "lower" = round(lower, digits = 2),
           "higher" = round(higher, digits = 2),
           'SE'  = round(SE,  digits=2),
           'PVAL' = formatC(PVAL, format = "e", digits = 1)) |> 
    flextable() |>
    bg(i = 1, bg = "#EFEFEF", part = "head") %>%
    bold(bold = TRUE, part="header") %>%
    align(align = 'center', part = "all") %>%
    align(j = 14, align = 'left', part = 'body') %>%
    width(j = c('CHR','EAF','EA','OA','OR','SE'), width = 1.2, unit = 'cm') %>%
    width(j = c('SNP','BP','PVAL'), width = 1.7, unit = 'cm') %>%
    width(j = c('Function','Gene','Distance (b)'), width = 2, unit = 'cm') %>%
    width(j = c('Previous reported trait'), width = 5, unit = 'cm') %>%
    hline(part = "body") |>
    fontsize(size = 11, part = "all") |>
    font(fontname = "Calibri", part = "all")
  
  save_as_image(tab, path = paste0(dir_results, "GWAS_Tables/",name,".png"), res = 400)
  save_as_docx(tab, path = paste0(dir_results, "GWAS_Tables/",name,".docx"))
  
}

# Validation figure ----
validation <- read_table(paste0(dir_data,"Results/GWAS/LongCovid_validation.txt"))
gwas <- read_table(paste0(dir_data,"Results/GWAS/LongCovid.txt"))
name <- "LongCovid"

path <- paste0(dir_results, "FUMA/",name,"/")
leadSnps    <- read_delim(paste0(path,'leadSNPs.txt'))
annovar     <- read_delim(paste0(path,'annov.txt'))
genomicRL   <- read_delim(paste0(path,'GenomicRiskLoci.txt'))

genomicRL <- genomicRL |> 
  select(SNP = rsID) %>%
  left_join(leadSnps |>
              select('uniqID',
                     'CHR' = 'chr',
                     'BP'  = 'pos',
                     'SNP' = 'rsID'),
            by = 'SNP') |>
  left_join(annovar |>
              group_by(uniqID) |>
              summarise('symbol' = paste0(symbol, collapse=" - ")) |>
              ungroup(),
            by = "uniqID")  |>
  select(rsID = SNP, symbol)

tab <- genomicRL |> 
  rename(ID = rsID) |>
  left_join(gwas, by = "ID") |>
  rename_all(~paste0(.,"_Main")) |>
  rename("ID" = "ID_Main") |>
  left_join(
    validation |>
      select("ID", 
             "ALLELE0_Val" = "ALLELE0", 
             "ALLELE1_Val" = "ALLELE1", 
             "A1FREQ_Val" = "A1FREQ", 
             "N_Val" = "N", 
             "BETA_Val" = "BETA", 
             "SE_Val" = "SE", 
             "LOG10P_Val" = "LOG10P"),
    by = "ID"
  ) |>
  mutate(A1FREQ_Val = if_else(ALLELE0_Val != ALLELE0_Main, 1-A1FREQ_Val, A1FREQ_Val),
         BETA_Val = if_else(ALLELE0_Val != ALLELE0_Main, -BETA_Val, BETA_Val)) |>
  mutate(PVAL_Val = 10^(-LOG10P_Val)) |>
  select(CHROM_Main, ID, symbol_Main, A1FREQ_Main, A1FREQ_Val, BETA_Main, BETA_Val, SE_Main, SE_Val, PVAL_Main, PVAL_Val) 

tab_long <- tab |>
  mutate(col = case_when(
    sign(BETA_Main) == sign(BETA_Val) & PVAL_Val < 0.05 ~ "Fully validated",
    sign(BETA_Main) == sign(BETA_Val) & PVAL_Val < 0.1 & PVAL_Val >= 0.05 ~ "Partially validated",
    sign(BETA_Main) == sign(BETA_Val) & PVAL_Val > 0.1 ~ "Partially validated",
    sign(BETA_Main) != sign(BETA_Val) ~ "Not validated"
  )) |>
  mutate(col = factor(col, levels = c("Fully validated", "Partially validated", "Not validated"))) |>
  pivot_longer(cols = c(starts_with("BETA"), starts_with("SE"), starts_with("A1FREQ"), starts_with("PVAL")),
               names_to = c(".value", "group"),
               names_pattern = "(.*)_(.*)") |>
  mutate(group = gsub("Val","Validation", group)) 

tab_long <- tab_long |>
  mutate(ID = factor(ID, levels = rev(unique(tab_long$ID)))) |>
  mutate(group = factor(group, levels = c("Main","Validation"))) |>
  arrange(ID, group) |>
  mutate(row_num = match(ID, unique(ID))) |>
  mutate(row_num_offset = ifelse(group == "Main", row_num + 0.1, row_num - 0.1)) |>
  mutate(bg = as.character(row_num %% 2))

p <- tab_long |>
  ggplot(aes(y = row_num_offset, x = exp(BETA), color = col, shape = group)) +
  geom_point(size = 3) +
  geom_linerange(aes(xmin = exp(BETA - 1.96 * SE),
                     xmax = exp(BETA + 1.96*SE))) +
  scale_y_continuous(breaks = tab_long$row_num,
                     labels = tab_long$ID,
                     sec.axis = sec_axis(~.,
                                         breaks = tab_long$row_num,
                                         labels = tab_long$symbol_Main),
                     expand = c(0,0)) +
  scale_x_log10(limits = c(0.4, 3.5),
                labels = c(0.5, 1, 1.75, 2),
                breaks = c(0.5, 1, 1.75, 3),
                expand = c(0,0),
                oob = scales::rescale_none) +
  xlab("Odds Ratio (OR)") +
  ylab("") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.box = "vertical",
        legend.spacing.y = unit(0.001, "lines"),
        strip.text.x = element_text(size = 11, face = "bold"),
        axis.title.x = element_text(),
        axis.text.x  = element_text(),
        axis.text.y  = element_text(size = 10),
        axis.title.y = element_blank(),
        text = element_text(family = 'Calibri')) +
  geom_vline(aes(xintercept = 1), 
             linetype = "dashed") +
  scale_color_manual(values = c("Fully validated" = "#7BB2B7",
                                "Partially validated" = "#CF731F",
                                "Not validated" = "#DAB92F")) +
  guides(
    color = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 2, nrow = 1),
    fill = "none"
  ) +
  geom_rect(aes(xmin = 0.4, xmax = 3.5,
                ymin = row_num-0.5, ymax = row_num + 0.5, fill = bg),
            alpha = 0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = c("1" = "#c1c1c1","0" = "white"))
    
ggsave(paste0(dir_results,'GWAS_Tables/Validation.png'), plot = p, width = 20, height = 30, dpi = 600, units = 'cm')

# Table one -----
ukb_t1 <- loadUKBTableOne(ukb)

cohort <- as_tibble(read_delim(paste0(dir_data,"/Results/LongCOVID_cohort.phe"))) |>
  select("eid" = "FID", "age", "Overall" = "state")  |>
  mutate(`Overall` = if_else(`Overall` == 1, "Cases", "Controls")) |>
  left_join(ukb_t1, by = "eid") |>
  mutate(Sex = if_else(sex == 0, "Female", "Male")) |>
  select(-c("sex", "year_of_birth")) |>
  mutate(`Overall` = as.factor(`Overall`)) |>
  select("eid", "Overall", "Age" = "age", "Sex", "Body mass index" = "body_mass_index",
         "Index of multiple deprivation" = "index_of_multiple_deprivation", "Genetic ethnic background" = "genetic_ethnic_grouping")

tableOne_lc <- CreateTableOne(data = cohort,
                           vars = colnames(cohort[2:length(colnames(cohort))])) |>
  print(showAllLevels = TRUE) 

tb_lc <- tibble(data.frame(tableOne_lc)) |>
  mutate(Variables = rownames(tableOne_lc)) |>
  relocate(Variables, .before = "level") |>
  rename("Level" = "level") |>
  mutate(Variables = if_else(Variables == "n", "N", Variables))

cohort <- as_tibble(read_delim(paste0(dir_data,"/Results/Pacs_cohort.phe"))) |>
  select("eid" = "FID", "age", "Overall" = "state")  |>
  mutate(`Overall` = if_else(`Overall` == 1, "Cases", "Controls")) |>
  left_join(ukb_t1, by = "eid") |>
  mutate(Sex = if_else(sex == 0, "Female", "Male")) |>
  select(-c("sex", "year_of_birth")) |>
  mutate(`Overall` = as.factor(`Overall`)) |>
  select("eid", "Overall", "Age" = "age", "Sex", "Body mass index" = "body_mass_index",
         "Index of multiple deprivation" = "index_of_multiple_deprivation", "Genetic ethnic background" = "genetic_ethnic_grouping")

tableOne_pacs <- CreateTableOne(data = cohort,
                              vars = colnames(cohort[2:length(colnames(cohort))])) |>
  print(showAllLevels = TRUE) 

tb_pacs <- tibble(data.frame(tableOne_pacs)) |>
  mutate(Variables = rownames(tableOne_pacs)) |>
  relocate(Variables, .before = "level") |>
  rename("Level" = "level") |>
  mutate(Variables = if_else(Variables == "n", "N", Variables))


tb_lc |>
  rename("Long COVID" = "Overall") |>
  left_join(tb_pacs |>
              rename("PACS" = "Overall")) |>
  flextable() |>
  width(j = 1, width = 6, unit = "cm") |>
  width(j = 3, width = 3, unit = "cm") |>
  width(j = 4, width = 3, unit = "cm") |>
  save_as_docx(path = paste0(dir_results, "Tables/TableOne.docx"))


cohort <- as_tibble(read_delim(paste0(dir_data,"/Results/LongCOVID_validation_cohort.phe"))) |>
  select("eid" = "FID", "age", "Overall" = "state")  |>
  mutate(`Overall` = if_else(`Overall` == 1, "Cases", "Controls")) |>
  left_join(ukb_t1, by = "eid") |>
  mutate(Sex = if_else(sex == 0, "Female", "Male")) |>
  select(-c("sex", "year_of_birth")) |>
  mutate(`Overall` = as.factor(`Overall`)) |>
  select("eid", "Overall", "Age" = "age", "Sex", "Body mass index" = "body_mass_index",
         "Index of multiple deprivation" = "index_of_multiple_deprivation", "Genetic ethnic background" = "genetic_ethnic_grouping")
tableOne_validation <- CreateTableOne(data = cohort,
                                vars = colnames(cohort[2:length(colnames(cohort))])) |>
  print(showAllLevels = TRUE) 
tb_validation <- tibble(data.frame(tableOne_pacs)) |>
  mutate(Variables = rownames(tableOne_pacs)) |>
  relocate(Variables, .before = "level") |>
  rename("Level" = "level") |>
  mutate(Variables = if_else(Variables == "n", "N", Variables))
