
quality_control <- function(gwas){
  gwas |> 
    dplyr::filter(!is.na(LOG10P)) |> 
    dplyr::filter(nchar(ALLELE0) == 1, nchar(ALLELE1) == 1) |>
    dplyr::filter(INFO >= 0.8) |>
    dplyr::mutate(dif = 1-A1FREQ) |>
    dplyr::filter(!(A1FREQ < 0.01 | dif < 0.01)) |>
    dplyr::select(-"dif")   
}

fuma_format <- function(gwas, filename, dir_data){
  gwas <- gwas |> as_tibble() 
  colnames(gwas) <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ",
                      "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "PVAL")
  
  gwas <- gwas |> quality_control()
  write_delim(gwas, file = paste0(dir_data,"Results/GWAS/", filename, ".txt"),
              col_names = TRUE,
              progress = TRUE,
              quote = "none", 
              delim = "\tab")
}


# gwas <- read_table(paste0(dir_data,"Results/GWAS/LongCovid_assoc.regenie.merged.txt"))
# fuma_format(gwas, "LongCovid", dir_data)
# 
# gwas <- read_table(paste0(dir_data,"Results/GWAS/Pacs_assoc.regenie.merged.txt"))
# fuma_format(gwas, "Pacs", dir_data)

gwas <- read_table(paste0(dir_data,"Results/GWAS/LC_Subtypes/assoc.cardio.regenie.merged.txt"))
fuma_format(gwas, "Cardio", dir_data)

gwas <- read_table(paste0(dir_data,"Results/GWAS/LC_Subtypes/assoc.ent.regenie.merged.txt"))
fuma_format(gwas, "ENT", dir_data)

gwas <- read_table(paste0(dir_data,"Results/GWAS/LC_Subtypes/assoc.fatigue.regenie.merged.txt"))
fuma_format(gwas, "Fatigue", dir_data)

gwas <- read_table(paste0(dir_data,"Results/GWAS/LC_Subtypes/assoc.neuro.regenie.merged.txt"))
fuma_format(gwas, "Neuro", dir_data)

gwas <- read_table(paste0(dir_data,"Results/GWAS/LC_Subtypes/assoc.overall.regenie.merged.txt"))
fuma_format(gwas, "Overall", dir_data)

gwas <- read_table(paste0(dir_data,"Results/GWAS/"))
fuma_format(gwas, "Overall", dir_data)

plot_helper <- function(gwas, name) {
  # Separate chromosomes x and xy
  gwas23 <- gwas |> filter(CHROM == 23) 
  x <- which(gwas23$GENPOS[1:length(gwas23$GENPOS)-1] > gwas23$GENPOS[2:length(gwas23$GENPOS)])
  
  gwas <- gwas |>
    filter(CHROM != 23) |>
    rbind(
      gwas23 |> mutate(CHROM = if_else(row_number() > x, 24, .data$CHROM))
    ) 
  
  gwas |>
    write_delim(here(paste0("Results/",name,".txt")), delim = "\t")
  
  count_snps <- tibble(
    step = "initial",
    counts = gwas |> tally() |> pull()
  )
  
  gwas <- gwas |>
    dplyr::mutate(P = 10^(-LOG10P)) %>%
    filter(!is.na(LOG10P)) 
  
  count_snps <- union_all(
    count_snps,
    tibble(
      step = "LOG10P not NA",
      counts = gwas |> tally() |> pull()
    ))
  
  if(("INFO" %in% colnames(gwas))) {
    gwas <- gwas |>
      filter(INFO > 0.8)   
    
    count_snps <- union_all(
      count_snps,
      tibble(
        step = "INFO greater than 0.8",
        counts = gwas |> tally() |> pull()
      ))
    
    gwas <- gwas |>
      mutate(dif = 1-A1FREQ) |>
      filter(!(A1FREQ < 0.01 | dif < 0.01)) |>
      select(-"dif")   
    
    count_snps <- union_all(
      count_snps,
      tibble(
        step = "MAF greater than 0.01",
        counts = gwas |> tally() |> pull()
      ))
  } 
  
  gwas |>
    filter(LOG10P > 7) |>
    print()
  
  count_snps |>
    write_csv(here(paste0("Results/",name,"_post_processing_filtering.txt")))
  
  plot1 <- getManhattanPlot(gwas, 8) #ylim?
  plot2 <- getQQPlot(gwas, x_lim = 6.5, y_lim = 8)
  
  plot1 + plot2
  ggsave(here(paste0("Results/",name,".png")), width = 16)
  
}

plot_helper(gwas, "LongCovidKim") # SNP rs200784681 (CHROM 12, GENPOS 129417693)
plot_helper(gwas_wide, "LongCovid_wide") # SNP rs6453106 (CHROM 5, GENPOS 74298946)
plot_helper(gwas_old, "LongCovid_old") # SNP rs200784681 (CHROM 12, GENPOS 129417693)
plot_helper(gwas_pacs, "Pacs") # no SNP 
plot_helper(gwas_pacs_a, "Pacs_arterial") # no SNP 
plot_helper(gwas_pacs_v, "Pacs_venous") # no SNP 

# no difference 0.7 SNPs ??

