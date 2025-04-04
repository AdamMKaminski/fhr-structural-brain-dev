

# read in master df (for FHR analysis)
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
df <- read.csv('master_df.csv')
df_1 <- df[df$time=='T1',]
df_2 <- df[df$time=='T2',]

# read in demos (for dropout analysis)
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
demos <- as.data.frame(read_xlsx('VIA11-15_longitudinal_demograph_clinic_wide.xlsx'))

# adjust variable types
# age:
demos$via11_mri_age <- as.numeric(demos$via11_mri_age)
demos$via15_mri_age <- as.numeric(demos$via15_mri_age)
# site:
demos$via11_mri_site[demos$via11_mri_site=="NA"] <- NA
demos$via11_mri_site <- as.factor(demos$via11_mri_site)
demos$via15_mri_site[demos$via15_mri_site=="NA"] <- NA
demos$via15_mri_site <- as.factor(demos$via15_mri_site)
# sex:
demos$sex_string[demos$sex_string=="NA"] <- NA
demos$sex_string <- as.factor(demos$sex_string)
# FHR group:
demos$fhr_group_string[demos$fhr_group_string=="NA"] <- NA
demos$fhr_group_string <- as.factor(demos$fhr_group_string)
# ADHD KSADS:
demos$KSADS_adhd_lft_v11[demos$KSADS_adhd_lft_v11=="NA"] <- NA
demos$KSADS_adhd_lft_v11 <- as.factor(demos$KSADS_adhd_lft_v11)
# Any KSADS:
demos$KSADS_any_diag_lft_v11[demos$KSADS_any_diag_lft_v11=="NA"] <- NA
demos$KSADS_any_diag_lft_v11 <- as.factor(demos$KSADS_any_diag_lft_v11)

# Create variable for whether sub is in or out of dataset (for VIA11 and VIA15 separately)
demos$in_analysis_11 <- as.numeric(demos$famlbnr %in% df_1$sub_id)
demos$in_analysis_11[demos$in_analysis_11==1] <- "In"
demos$in_analysis_11[demos$in_analysis_11==0] <- "Out"
demos$in_analysis_11 <- as.factor(demos$in_analysis_11)
demos$in_analysis_15 <- as.numeric(demos$famlbnr %in% df_2$sub_id)
demos$in_analysis_15[demos$in_analysis_15==1] <- "In"
demos$in_analysis_15[demos$in_analysis_15==0] <- "Out"
demos$in_analysis_15 <- as.factor(demos$in_analysis_15)

# set cols to compare between in and out subjects
imp_cols_1 <- c('via11_mri_age','via11_mri_site','sex_string',
                'fhr_group_string','KSADS_adhd_lft_v11','KSADS_any_diag_lft_v11')
imp_cols_2 <- c('via15_mri_age','via15_mri_site')

sink("dropout_stats.txt")

for(col in imp_cols_1){
  print(paste('Testing ',col,sep=''))
  print(paste("In N=",length(demos[,col][(!is.na(demos[,col])&demos$in_analysis_11=="In")]),sep=""))
  print(paste("Out N=",length(demos[,col][(!is.na(demos[,col])&demos$in_analysis_11=="Out")]),sep=""))
  if(is.factor(demos[,col])){
    temp_table <- table(demos$in_analysis_11,demos[,col])
    print(temp_table)
    print(chisq.test(temp_table))
  } else if(is.numeric(demos[,col])){
    print(tapply(demos[,col], demos$in_analysis_11, function(x) c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE))))
    print(t.test(demos[,col]~demos$in_analysis_11))
  }}

for(col in imp_cols_2){
  print(paste('Testing ',col,sep=''))
  print(paste("In N=",length(demos[,col][(!is.na(demos[,col])&demos$in_analysis_15=="In")]),sep=""))
  print(paste("Out N=",length(demos[,col][(!is.na(demos[,col])&demos$in_analysis_15=="Out")]),sep=""))
  if(is.factor(demos[,col])){
    temp_table <- table(demos$in_analysis_15,demos[,col])
    print(temp_table)
    print(chisq.test(temp_table))
  } else if(is.numeric(demos[,col])){
    print(tapply(demos[,col], demos$in_analysis_15, function(x) c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE))))
    print(t.test(demos[,col]~demos$in_analysis_15))
  }}

# Do the same for CBCL and CGAS
allkey_11 <- read.csv('VIA11_allkey_050624.csv',sep=';')
allkey_15 <- read.csv('VIA15_allkey_291124.csv',sep=',')

# Create variable for whether sub is in or out of dataset (for VIA11 and VIA15 separately)
allkey_11$in_analysis_11 <- as.numeric(allkey_11$famlbnr %in% df_1$sub_id)
allkey_11$in_analysis_11[allkey_11$in_analysis_11==1] <- "In"
allkey_11$in_analysis_11[allkey_11$in_analysis_11==0] <- "Out"
allkey_11$in_analysis_11 <- as.factor(allkey_11$in_analysis_11)
allkey_15$in_analysis_15 <- as.numeric(allkey_15$famlbnr %in% df_2$sub_id)
allkey_15$in_analysis_15[allkey_15$in_analysis_15==1] <- "In"
allkey_15$in_analysis_15[allkey_15$in_analysis_15==0] <- "Out"
allkey_15$in_analysis_15 <- as.factor(allkey_15$in_analysis_15)

# set cols to compare between in and out subjects
imp_cols_1 <- c('CBCL_ext_cg_v11','CBCL_int_cg_v11','CBCL_totsc_cg_v11','CGASx_v11')
imp_cols_2 <- c('CBCL_totsc_cg_v15','CGASx_v15')

for(col in imp_cols_1){
  print(paste('Testing ',col,sep=''))
  print(paste("In N=",length(allkey_11[,col][(!is.na(allkey_11[,col])&allkey_11$in_analysis_11=="In")]),sep=""))
  print(paste("Out N=",length(allkey_11[,col][(!is.na(allkey_11[,col])&allkey_11$in_analysis_11=="Out")]),sep=""))
  print(tapply(allkey_11[,col], allkey_11$in_analysis_11, function(x) c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE))))
  print(t.test(allkey_11[,col]~allkey_11$in_analysis_11))
}

for(col in imp_cols_2){
  print(paste('Testing ',col,sep=''))
  print(paste("In N=",length(allkey_15[,col][(!is.na(allkey_15[,col])&allkey_15$in_analysis_15=="In")]),sep=""))
  print(paste("Out N=",length(allkey_15[,col][(!is.na(allkey_15[,col])&allkey_15$in_analysis_15=="Out")]),sep=""))
  print(tapply(allkey_15[,col], allkey_15$in_analysis_15, function(x) c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE))))
  print(t.test(allkey_15[,col]~allkey_15$in_analysis_15))
}
sink()



# validation
mean(as.numeric(demos$via11_mri_age[demos$in_analysis_11==1]),na.rm=TRUE)
mean(df$age[df$time=='T1'])
mean(as.numeric(demos$via15_mri_age[demos$in_analysis_15==1]),na.rm=TRUE)
mean(df$age[df$time=='T2'])

