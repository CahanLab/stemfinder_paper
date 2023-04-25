### One-way ANOVA and Tukey's HSD
# To compare multiple sample groups in robustness testing
# Example: robustness to K value in KNN 

#ANOVA one way: tells you if any groups are signif diff
df_robust_k$kratio = as.factor(df_robust_k$kratio) #tukey function requires groups as factors
one.way = aov(data = df_robust_k[df_robust_k$UMI == 1,], spear_all_sF ~ kratio)
summary(one.way)

#Tukey's HSD: tells you which groups are signif diff
tukey = TukeyHSD(one.way)
tukey


