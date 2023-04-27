## Pairwise T test to compare 2 sample groups

    #Example: compare stemFinder vs. CytoTRACE single-cell Spearman values
    #perform Bonferroni p value correction
stat.test <- df_results[df_results$method %in% c('stemFinder','CytoTRACE'),] %>% pairwise_t_test(spear_all ~ method, paired = T, p.adjust.method = 'bonferroni')
