
library(car)
library(lsmeans)
library(lsr)
library(BayesFactor)

data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220509.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")

## Preprocess 

#filter the data with include variable
# - extract rows with 1 in Include_FS_studies coloumn
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]

# - extract rows with 1 in Include_FS_studies_euler_outliers_sibpairs_out
data_csv_filtered <- data_csv_filtered[c(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out),]

#make new variables with shorter and contained names
# - tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$site = as.factor(datab$MRI_site_v11)
datab$diag = as.factor(datab$ksads_any_diag_excl_elim_lft_v11) #Axis.1_diag_v11) 
datab$age = as.numeric(datab$MRI_age)
datab$group = as.factor(datab$HighRiskStatus_v11)

datab$diag_alt = datab$HRS_SZ_K_axis1 # HAVE TO CHANGE THIS
datab$diag_alt[datab$diag_alt == 3] = "K"
datab$diag_alt[datab$diag_alt == 2] = "SZ+"
datab$diag_alt[datab$diag_alt == 1] = "SZ-"
datab$diag_alt = as.factor(datab$diag_alt)
contrast_signs = c(1,1, 1) # for K/BP-/BP+
contrast_signs = c(-1,-1, 1) # for K/SZ-/SZ+
contrast_signs = c(1,1, -1) # for K/BP/SZ
contrast_signs = c(1,1, -1) # for K/BP/SZ 

# lms for any axis I vs. no axis I 
# for each risk group, and for females and males separately 

datab_SZ <- datab[datab$HighRiskStatus_v11 == 'SZ',]
datab_BP <- datab[datab$HighRiskStatus_v11 == 'BP',]
datab_PBC <- datab[datab$HighRiskStatus_v11 == 'K',]

datab_SZ_F <- datab_SZ[datab_SZ$sex == 0,] 
datab_SZ_M <- datab_SZ[datab_SZ$sex == 1,]

datab_BP_F <- datab_BP[datab_BP$sex == 0,] 
datab_BP_M <- datab_BP[datab_BP$sex == 1,]

datab_PBC_F <- datab_PBC[datab_PBC$sex == 0,] 
datab_PBC_M <- datab_PBC[datab_PBC$sex == 1,]

datab_list <- list(datab_SZ_F = datab_SZ_F, 
                   datab_SZ_M = datab_SZ_M, 
                   datab_BP_F = datab_BP_F, 
                   datab_BP_M = datab_BP_M, 
                   datab_PBC_F = datab_PBC_F, 
                   datab_PBC_M = datab_PBC_M)

# test area

datab_SZ_K_M <- datab[datab$sex==1,]
datab_SZ_K_M <- datab_SZ_K_M[datab_SZ_K_M$HighRiskStatus_v11 != 'BP',]

print('BrainTotalVol')
print(summary(lm(BrainTotalVol ~ diag_alt + age + site + TotalEulerNumber, data=datab_SZ_K_M)))
model_temp = lm(BrainTotalVol ~ diag_alt + age + site + TotalEulerNumber, data=datab_SZ_K_M)
Anova(model_temp,Type='III')
print(lsmeans(model_temp,pairwise ~ diag_alt, adjust="none"))

# end test area

# Group differences in: "BrainTotalVol", "CortexVol", "total_area", "eICV_samseg"

setwd('~/Desktop/VIA11_project/')
sink("pairwise_lsmeans.txt")
for(subgroup in 1:length(datab_list)){
  print(paste('Running models for ',names(datab_list[subgroup])),sep='')
  
  df_temp <- as.data.frame(cbind(datab_list[[subgroup]]$BrainTotalVol,
                   datab_list[[subgroup]]$diag,
                   datab_list[[subgroup]]$age,
                   datab_list[[subgroup]]$site,
                   datab_list[[subgroup]]$TotalEulerNumber,
                   datab_list[[subgroup]]$diag_alt))
  colnames(df_temp) <- c('brain','diag','age','site','euler','diag_alt')
  df_temp$diag <- as.factor(df_temp$diag)
  df_temp$site <- as.factor(df_temp$site)
  df_temp$diag_alt <- as.factor(df_temp$diag_alt)
  
  print('BrainTotalVol')
  print(summary(lm(brain ~ diag + age + site + euler, data=df_temp)))
  model_temp = lm(brain ~ diag + age + site + euler, data=df_temp)
  print(lsmeans(model_temp,pairwise ~ diag))
  
  print('CortexVol')
  df_temp$brain <- datab_list[[subgroup]]$CortexVol
  print(summary(lm(brain ~ diag + age + site + euler, data=df_temp)))
  model_temp = lm(brain ~ diag + age + site + euler, data=df_temp)
  print(lsmeans(model_temp,pairwise ~ diag))
  
  print('total_area')
  df_temp$brain <- datab_list[[subgroup]]$total_area
  print(summary(lm(brain ~ diag + age + site + euler, data=df_temp)))
  model_temp = lm(brain ~ diag + age + site + euler, data=df_temp)
  print(lsmeans(model_temp,pairwise ~ diag))
  
  print('eICV_samseg')
  df_temp$brain <- datab_list[[subgroup]]$eICV_samseg
  print(summary(lm(brain ~ diag + age + site + euler, data=df_temp)))
  model_temp = lm(brain ~ diag + age + site + euler, data=df_temp)
  print(lsmeans(model_temp,pairwise ~ diag))
}
sink()

# pairwise comparisons of demographic and symptom variables 

# Variables: CBCL_totsc_cg_v11, CBCL_int_cg_v11, CBCL_ext_cg_v11, CGASx_v11, diag 

# groups to compare
group1 <- "SZ"
group2 <- "BP"

group_df <- subset(datab, group %in% c(group1, group2))
t.test(CGASx_v11 ~ group, data = group_df)

group_df$group <- droplevels(group_df$group)
contingency_table <- table(group_df$diag,group_df$group)
chisq.test(contingency_table)

# imaging measure analyses

# in datab variables of interest are:
# split dfs by "sex"
# control for "age", "site", "TotalEulerNumber", "diag"
# interested in effect of "group"
# effect on: "BrainTotalVol", "CortexVol", "total_area", "eICV_samseg"
# also possibly: "mean_thickness"

# run model -- replication of previous results

model_bvol_glob = lm(eICV_samseg ~ group*sex + age + site + TotalEulerNumber, data=datab)
#Anova(model_bvol_glob,type="III")

mf1 <- 'eICV_samseg ~ group + age + sex + site + TotalEulerNumber'
mf2 <- 'eICV_samseg ~ age + sex + site + TotalEulerNumber'
bf_iterations = 100000
bf1 = lmBF(formula = eval(parse(text=mf1)), data=datab)
bf2 = lmBF(formula = eval(parse(text=mf2)), data=datab)
bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
bf_g = as.numeric(bf_res[1])
bf_g_error = as.numeric(bf_res[2])
bf_g

etaSquared(model_bvol_glob, type = 3, anova = T)


# split by sex

datab_female <- datab[datab$sex == 0,] # 0 = female
datab_male <- datab[datab$sex == 1,] # 1 = male

# run models -- replication of previous results

model_bvol_glob = lm(BrainTotalVol ~ group + age + site + TotalEulerNumber, data=datab_female)
Anova(model_bvol_glob,type="III")
lsmeans(model_bvol_glob,pairwise ~ group)

# run models -- new

setwd('~/Desktop')
sink('ancova_results_summary.txt')

# BrainTotalVol, female
print('BrainTotalVol, FEMALE')
model_bvol_glob_1 = lm(BrainTotalVol ~ group + age + site + TotalEulerNumber + diag, data=datab_female)
Anova(model_bvol_glob_1,type="III")
lsmeans(model_bvol_glob_1,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_1,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_1), edf=df.residual(model_bvol_glob_1))

# BrainTotalVol, male
print('BrainTotalVol, MALE')
model_bvol_glob_2 = lm(BrainTotalVol ~ group + age + site + TotalEulerNumber + diag, data=datab_male)
Anova(model_bvol_glob_2,type="III")
lsmeans(model_bvol_glob_2,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_2,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_2), edf=df.residual(model_bvol_glob_2))

# CortexVol, female
print('CortexVol, FEMALE')
model_bvol_glob_3 = lm(CortexVol ~ group + age + site + TotalEulerNumber + diag, data=datab_female)
Anova(model_bvol_glob_3,type="III")
lsmeans(model_bvol_glob_3,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_3,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_3), edf=df.residual(model_bvol_glob_3))

# CortexVol, male
print('CortexVol, MALE')
model_bvol_glob_4 = lm(CortexVol ~ group + age + site + TotalEulerNumber + diag, data=datab_male)
Anova(model_bvol_glob_4,type="III")
lsmeans(model_bvol_glob_4,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_4,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_4), edf=df.residual(model_bvol_glob_4))

# total_area, female
print('total_area, FEMALE')
model_bvol_glob_5 = lm(total_area ~ group + age + site + TotalEulerNumber + diag, data=datab_female)
Anova(model_bvol_glob_5,type="III")
lsmeans(model_bvol_glob_5,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_5,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_5), edf=df.residual(model_bvol_glob_5))

# total_area, male
print('total_area, MALE')
model_bvol_glob_6 = lm(total_area ~ group + age + site + TotalEulerNumber + diag, data=datab_male)
Anova(model_bvol_glob_6,type="III")
lsmeans(model_bvol_glob_6,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_6,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_6), edf=df.residual(model_bvol_glob_6))

# eICV_samseg, female
print('eICV_samseg, FEMALE')
model_bvol_glob_7 = lm(eICV_samseg ~ group + age + site + TotalEulerNumber + diag, data=datab_female)
Anova(model_bvol_glob_7,type="III")
lsmeans(model_bvol_glob_7,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_7,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_7), edf=df.residual(model_bvol_glob_7))

# eICV_samseg, male
print('eICV_samseg, MALE')
model_bvol_glob_8 = lm(eICV_samseg ~ group + age + site + TotalEulerNumber + diag, data=datab_male)
Anova(model_bvol_glob_8,type="III")
lsmeans(model_bvol_glob_8,pairwise ~ group)
emm_eff = emmeans(model_bvol_glob_8,specs="group")
eff_size(emm_eff, sigma=sigma(model_bvol_glob_8), edf=df.residual(model_bvol_glob_8))

sink()

# put results in df

tab <- data.frame(matrix(0,8,11))
models_list <- list(model_bvol_glob_1,model_bvol_glob_2,model_bvol_glob_3,model_bvol_glob_4,model_bvol_glob_5,model_bvol_glob_6,model_bvol_glob_7,model_bvol_glob_8)

for(i in 1:length(models_list)){
  temp <- as.data.frame(lapply(models_list[i], Anova, type = "III"))
  tab[i,1] <- temp$F.value[2]
  tab[i,2] <- temp$Pr..F.[2]
  tab[i,3] <- temp$F.value[3]
  tab[i,4] <- temp$Pr..F.[3]
  tab[i,5] <- temp$F.value[4]
  tab[i,6] <- temp$Pr..F.[4]
  tab[i,7] <- temp$F.value[5]
  tab[i,8] <- temp$Pr..F.[5]
  tab[i,9] <- temp$F.value[6]
  tab[i,10] <- temp$Pr..F.[6]
  if(i %% 2 == 0) {
    tab[i,11] <- 1
  } else {
    tab[i,11] <- 0
  }
}



