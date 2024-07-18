# Formatting results of robustness testing

df_robust = as.data.frame(list_all) #convert list to df

# Add formatted dataset names
load("data_names.rda")
df_robust$dataset_forplot = NA
for(d in 1:length(datasets_forplot)){
  df_robust[df_robust$dataset == datasets_forplot[[d]],]$dataset_forplot = names(datasets_forplot[d])
}

# Calculate negative deviation from performance at default parameter value
    ##example: downsampling 
df_robust$neg_deviation_auc = NA
df_robust$neg_deviation_spearss = NA
df_robust$neg_deviation_spearp = NA
df_robust$neg_deviation_pctrecov = NA

for (d in df_robust$dataset){
  df_robust[df_robust$dataset == d,]$neg_deviation_auc = df_robust[df_robust$dataset == d,]$auc_sF - df_robust[df_robust$ds_ratio == 1 & df_robust$dataset == d, ]$auc_sF
  df_robust[df_robust$dataset == d,]$neg_deviation_spearss = df_robust[df_robust$dataset == d,]$spear_all_sF - df_robust[df_robust$ds_ratio == 1 & df_robust$dataset == d, ]$spear_all_sF
  df_robust[df_robust$dataset == d,]$neg_deviation_spearp = df_robust[df_robust$dataset == d,]$spear_pheno_sF - df_robust[df_robust$ds_ratio == 1 & df_robust$dataset == d, ]$spear_pheno_sF
  df_robust[df_robust$dataset == d,]$neg_deviation_pctrecov = df_robust[df_robust$dataset == d,]$pct.recov_sF - df_robust[df_robust$ds_ratio == 1 & df_robust$dataset == d, ]$pct.recov_sF
}


