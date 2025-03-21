
# analysis

library(readODS)
library(readxl)
library(tidyr)
library(lme4)
library(lmerTest)
library(ggplot2)

# ==== FUNCTIONS ====

# function to display file name and get user input on whether to read it
read_file_or_not <- function(file_name) {
  selected_input <- readline(prompt = paste("Read ",file_name,"? (Enter Y or N) ",sep=""))
  if(selected_input == 'Y'){
    return(TRUE)
  }
  else if(selected_input == 'N'){
    return(FALSE)
  }
  else {
    print('Please enter Y or N')
    break
  }
}

# function to report on quantity of longitudinal data
# df must have "sub_id" and "time" cols, where time is either "ses01" or "ses02"
report <- function(sublist_df) {
  uniq_subs <- length(unique(sublist_df$sub_id))
  long_subs <- sum(duplicated(sublist_df$sub_id))
  if(all(c("T1","T2") %in% sublist_df$time)){
    t1_subs <- sum(sublist_df$time[!duplicated(sublist_df$sub_id) & !duplicated(sublist_df$sub_id, fromLast = TRUE)]=="T1")
    t2_subs <- sum(sublist_df$time[!duplicated(sublist_df$sub_id) & !duplicated(sublist_df$sub_id, fromLast = TRUE)]=="T2")
  } else if(all(c("ses01","ses02") %in% sublist_df$time)){
    t1_subs <- sum(sublist_df$time[!duplicated(sublist_df$sub_id) & !duplicated(sublist_df$sub_id, fromLast = TRUE)]=="ses01")
    t2_subs <- sum(sublist_df$time[!duplicated(sublist_df$sub_id) & !duplicated(sublist_df$sub_id, fromLast = TRUE)]=="ses02")
  } else {
    return("something off in time col")
  }
  return(cat('Total number of unique subjects: ',uniq_subs,
      '\nTotal number of subjects with long data (T1 and T2): ',long_subs,
      '\nTotal number of subjects with T1 data only: ',t1_subs,
      '\nTotal number of subjects with T2 data only: ',t2_subs,
      '\nTotal T1 data points: ',long_subs + t1_subs,
      '\nTotal T2 data points: ', long_subs + t2_subs,
      sep=''))
}

# function to display column numbers and names and get user input for selection
select_columns_console <- function(df) {
  # Display column numbers and names to the user
  cat("Columns in the dataframe:\n")
  for (i in 1:ncol(df)) {
    cat(i, ": ", colnames(df)[i], "\n", sep = "")
  }
  
  # Ask the user to input the column numbers they want to select
  selected_input <- readline(prompt = "Enter the column numbers you want to select, separated by commas: ")
  
  # Convert the input into a numeric vector
  selected_indices <- as.numeric(strsplit(selected_input, ",")[[1]])
  
  # Get the selected column names based on the user input
  selected_columns <- colnames(df)[selected_indices]
  
  # Return the selected column names
  return(selected_columns)
}

# function to match subs and time points in a new df to a reference df
# dfs are in long format; ref_df is "sublist", and new_df is whatever is new
# colnames must be "sub_id" and "time" with 'T1' or 'T2'
match_sub_time <- function(ref_df,new_df) {
  # if ref_df doesn't have the id_T column already
  if(!('id_T' %in% colnames(ref_df))) {
    ref_df$time[ref_df$time=='ses01'] <- 'T1'
    ref_df$time[ref_df$time=='ses02'] <- 'T2'
    ref_df$id_T <- paste(ref_df$sub_id,ref_df$time,sep="_") }
  
  # keep only correct time-points
  new_df$id_T <- paste(new_df$sub_id,new_df$time,sep="_")
  new_df <- new_df[new_df$id_T %in% ref_df$id_T,]
  if(all(new_df$sub_id==ref_df$sub_id)){
    if(all(new_df$id_T==ref_df$id_T)){
      print('solid matching dfs') 
      print('success')
      return(list(ref_df = ref_df, new_df = new_df))}
    } else {stop('Error: dfs do not match !')}}

# ==== siblings ====

# READ; NO NEED TO RUN THIS SECTION

# this is how I created a df to find sub_ids who were siblings as well as who
# they were matched to (i.e., who their sibling was). From there, I worked in
# excel to determine which sibling to include -- based on, 1) is there 
# longitudinal data for one and not the other? 2) does one have better image
# quality ratings? and finally 3) if they're equal on 1 and 2, just use the 
# sibling that came in first. 

# The results of this work have been saved in sub_list.ods (see next section)
# THEREFORE, NO NEED TO RUN THIS SECTION

# read in complete list of subs
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1')
sublist <- read_ods('sub_list.ods')
sublist <- sublist[sublist$image_incl==1,]
sublist <- sublist[sublist$qc_incl==1,]
sublist <- sublist[sublist$base_out_incl==1,]

# exclude sibling pairs ***
# Go into allkey (can use from VIA11 project you were helping Kathrine on) use sib variables there
# - Use the sib with more data
# - how many actual sib pairs are included?
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220509.csv"
via11_data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")
all(sublist$sub_id %in% via11_data_csv$famlbnr)

# sib variables are:
# Sib_ID_v11
# Sib_pair_no_v11 - Use this one, pairs each sibling
# Sibpair_v11
sib_dat <- via11_data_csv[,c('famlbnr','Sib_pair_no_v11')]
sib_dat <- sib_dat[sib_dat$famlbnr %in% sublist$sub_id,]
all(unique(sib_dat$famlbnr) == unique(sublist$sub_id))

x <- sib_dat$Sib_pair_no_v11
included_pairs <- unique(x[duplicated(x) & !is.na(x)])
sib_dat <- sib_dat[sib_dat$Sib_pair_no_v11 %in% included_pairs,]
sib_dat <- sib_dat[order(sib_dat$famlbnr),]

# ==== combining data; CREATE DF master_df.csv FOR ANALYSIS ==== 

# read in complete list of subs
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1')
sublist <- read_ods('sub_list.ods')

# only keep keepers,
# exclude based on:
# original image include (i.e., was it run through FS at all? OR is it usable or CORTICAL?)
# my qc ratings -- based on external primarily (i.e., pre, postcentral misaligned)
# stat outliers -- based on internal primarily (base and long data)
# euler outliers
sublist <- sublist[sublist$image_incl==1,]
sublist <- sublist[sublist$qc_incl=="1",]
sublist <- sublist[sublist$base_out_incl=="1",]
sublist <- sublist[sublist$long_out_incl!="0",]
sublist <- sublist[sublist$euler_incl=="1",]    
sublist <- sublist[sublist$twin_pair_include=='1',]

# report final subject list:
sublist <- as.data.frame(sublist[,c(1,2)])
rownames(sublist) <- NULL
report(sublist)

# compare these subs with data you read in, does it all match?

# BASIC DEMOGRAPHICS
# get demo variables to control for
# only include subs in sublist
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
demos <- as.data.frame(read_xlsx('VIA11-15_longitudinal_demograph_clinic_wide.xlsx'))
demos <- demos[demos$famlbnr %in% sublist$sub_id,]
# demos is in wide, must convert to long
# first get only important vars
# and clean them:
# ADD CONTROL VARIABLES HERE:   ... what else to control for? Psychotic-like experiences? Weight? Height?? Handedness? 
imp_cols <- c('via11_mri_age','via15_mri_age',
              'via11_mri_site','via15_mri_site',
              'sex_string','sex_code',
              'fhr_group_string','fhr_group_code',
              'KSADS_adhd_lft_v11','KSADS_any_diag_lft_v11')
# cleaning variables (i.e., string "NA" should be actual NA; numbers should be numeric):
for(imp_col in imp_cols){
  demos[[imp_col]] <- ifelse(demos[[imp_col]]=="NA",NA,demos[[imp_col]])
  if(imp_col == 'via11_mri_age' | imp_col == 'via15_mri_age'){
    demos[[imp_col]] <- as.numeric(demos[[imp_col]]) }}
# get new df
demos <- demos[,c("famlbnr",imp_cols)]
# pivot from wide to long
# IMPORTANT:
# !!!!!
# CHANGE IF ADDING ADDITIONAL CONTROL VARIABLES
# !!!!!
colnames(demos)
# very important naming convention for each column
colnames(demos) <- c('sub_id','age_T1','age_T2',
                     'site_T1','site_T2','sex_str',
                     'sex','FHR_str','FHR',
                     'VIA11_adhd_lft','VIA11_any_dx_lft')
demos_long <- as.data.frame(pivot_longer(demos, 
                    cols = c(age_T1, age_T2, site_T1, site_T2), 
                    names_to = c(".value", "time"),
                    names_sep = "_"))
# adjust and double check dfs are matching
df_list <- match_sub_time(sublist,demos_long)
sublist <- df_list$ref_df
demos_long <- df_list$new_df

# BRAIN DATA
# set directory for data
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
brain_long <- sublist[,c("sub_id","time")]
# loop through brain data files and select columns to add to your df, e.g.,
# lh_WhiteSurfArea_area
# rh_WhiteSurfArea_area
# lh_MeanThickness_thickness
# rh_MeanThickness_thickness
# BrainSegVolNotVent
# eTIV
# long-ish loop, double checks to make sure IDs and TIME matches across dfs
for(f in list.files(pattern='_long_')) {
  
  if(read_file_or_not(f)) {
    temp_brain_dat <- read.csv(f)
    
    temp_brain_dat$session_id[temp_brain_dat$session_id=='ses01'] <- 'T1'
    temp_brain_dat$session_id[temp_brain_dat$session_id=='ses02'] <- 'T2'
    temp_brain_dat$id_T <- paste(temp_brain_dat$famlbnr,temp_brain_dat$session_id,sep="_")
    temp_brain_dat <- temp_brain_dat[temp_brain_dat$id_T %in% sublist$id_T,]
    temp_brain_dat <- temp_brain_dat[order(temp_brain_dat$famlbnr),]
    
    if(all(temp_brain_dat$famlbnr==sublist$sub_id)) {
      print('brain data sub ids match sublist ids')
      if(all(temp_brain_dat$session_id==sublist$time)) {
          print('brain data TIME matches sublist time, proceeding...')
        
          selected_cols <- select_columns_console(temp_brain_dat)
          print(selected_cols)
          
          for(m in selected_cols) {
            brain_long <- cbind(brain_long,temp_brain_dat[,m])
            colnames(brain_long)[ncol(brain_long)] <- colnames(temp_brain_dat[m])
          } 
    
          print(paste('Added cols: ',selected_cols,sep=''))
      } else { 
        print('brain data TIME does not match sublsit TIME, there is an issue') 
        } 
    } else {
      print('brain data IDs do not match sublist ids, there is an issue')
      print(paste('issue for ',f,sep=''))
    }
  }
}
# average the lh and rh measures to get the global
cols <- colnames(brain_long)
lh_cols <- grep("^lh_", cols, value = TRUE)
rh_cols <- grep("^rh_", cols, value = TRUE)
# loop: determine lh, rh pairs, then for each pair average them
for (lh in lh_cols) {
  base_name <- sub("^lh_", "", lh)
  rh_match <- grep(paste0("^rh_", base_name), rh_cols, value = TRUE)
  if (length(rh_match) > 0) {
    print(paste('averaging ',lh,' and ',rh_match,sep=''))
    brain_long[[paste(base_name,'_global',sep='')]] <- (brain_long[,lh] + brain_long[,rh_match])/2 }}
# adjust and double check dfs are matching
# the new df MUST have a "sub_id" col and a "time" col with "T1" and "T2"
df_list <- match_sub_time(sublist,brain_long)
sublist <- df_list$ref_df
brain_long <- df_list$new_df

# EULER NUMBER
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
fs_qa_long <- read.csv('VIA11-15_fs741_long_qa_measures_restructured.csv')
fs_qa_long <- fs_qa_long[fs_qa_long$famlbnr %in% sublist$sub_id,]
fs_qa_long <- fs_qa_long[order(fs_qa_long$famlbnr),]
fs_qa_long$sub_id <- fs_qa_long$famlbnr
fs_qa_long$time <- fs_qa_long$session_id
fs_qa_long$time[fs_qa_long$time=='ses01'] <- 'T1'
fs_qa_long$time[fs_qa_long$time=='ses02'] <- 'T2'
df_list <- match_sub_time(sublist,fs_qa_long)
sublist <- df_list$ref_df
fs_qa_long <- df_list$new_df

# VIA11 DENSITY
setwd('/home/adamk/Desktop/VIA11_project/')
ak <- read_xlsx('VIA11_allkey-091221_WB-MRI-20220329_MRI-QC_combi_20220427.xlsx')
ak <- as.data.frame(ak[,c('famlbnr','REG_degurba_EU_v7')])
# T1 and T2 are NOT true here. It's a fixed variable (for my sake here...)
colnames(ak) <- c('sub_id','urb_T1')
ak$urb_T2 <- ak$urb_T1
ak <- as.data.frame(pivot_longer(ak, cols = c(urb_T1, urb_T2),names_to = c(".value", "time"),names_sep = "_"))
# adjust and make dfs matching
df_list <- match_sub_time(sublist,ak)
sublist <- df_list$ref_df
ak <- df_list$new_df
colnames(ak)[3] <- 'VIA7_density'

# CREATE MASTER DF FROM THE DFS YOU HAVE CREATED UP TO THIS POINT
# IF YOU NEED TO DO FROM SCRATCH (i.e., NEW SUB EXCLUSION)
master_df <- sublist[,c("sub_id","time")]
df_long_list <- list(demos_long,brain_long,fs_qa_long,ak)
# select cols from the dfs you've created
# sub_id and time are ALREADY in the master_df
for(dat_fr in df_long_list){
  match_sub_time(master_df,dat_fr)
  selected_cols <- select_columns_console(dat_fr)
  print(selected_cols)
  
  master_df <- cbind(master_df,dat_fr[,c(selected_cols)])
  if(length(selected_cols) == 1) {
    colnames(master_df)[ncol(master_df)] <- selected_cols }
}
# IF YOU WANT
# write the master df to your directory:
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
write.csv(master_df,file='master_df.csv',row.names=FALSE)

# IF YOU ALREADY HAVE A STARTING POINT
# read in data you already combined
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
master_df <- read.csv('master_df_TEST.csv')
if(all(master_df$id_T==ak$id_T)){
  master_df$VIA11_urbanicity <- ak$VIA11_density }
if(all(master_df$id_T==demos_long$id_T)){
  master_df$VIA11_adhd_lft <- demos_long$VIA11_adhd_lft 
  master_df$VIA11_any_dx_lft <- demos_long$VIA11_any_dx_lft }
# IF YOU WANT:
# write the master df to your directory
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
write.csv(master_df,file='master_df_TEST.csv',row.names=FALSE)

# ==== START - read in ====

# read in data
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
master_df <- read.csv('master_df.csv')

master_df$FHR <- as.factor(master_df$FHR)
master_df$FHR <- relevel(master_df$FHR, ref = 3)
master_df$FHR_str <- as.factor(master_df$FHR_str)
master_df$FHR_str <- relevel(master_df$FHR_str, ref = "PBC")

# ==== replications ====

# can you basically replicate previous VIA 11 results?
# it's a different sample subset so wouldn't expect identical results
# to replicate:
# - group-by-sex interactions in brain volume, cortical volume, and surface area 
# - FHR-SZ males - smaller volumes and surface area relative to male controls 
# - FHR-BP females - larger brain and cortical volumes than female controls

master_df_T1 <- master_df[master_df$time=="T1",]

# FHR-SZ males should show smaller BRAIN VOLUME than controls
summary(lm(BrainSegVolNotVent ~ FHR_str*sex + age + site + mean_eulnum_cross, data=master_df_T1)) # yes, signif interaction
temp <- master_df_T1[(master_df_T1$FHR_str!='BP' & master_df_T1$sex_str=='male'),]
t.test(temp$BrainSegVolNotVent~temp$FHR_str) # yes, signif t-test
temp <- master_df_T1[(master_df_T1$FHR_str!='BP' & master_df_T1$sex_str=='female'),]
t.test(temp$BrainSegVolNotVent~temp$FHR_str) # yes, NOT a signif t-test

# FHR-SZ males should show smaller CORTICAL VOLUME than controls
summary(lm(eTIV ~ FHR_str*sex + age + site + mean_eulnum_cross, data=master_df_T1)) # yes, signif interaction
temp <- master_df_T1[(master_df_T1$FHR_str!='BP' & master_df_T1$sex_str=='male'),]
t.test(temp$eTIV~temp$FHR_str) # yes, signif t-test
temp <- master_df_T1[(master_df_T1$FHR_str!='BP' & master_df_T1$sex_str=='female'),]
t.test(temp$eTIV~temp$FHR_str) # yes, NOT a signif t-test

# FHR-SZ males should show smaller SURFACE AREA than controls
summary(lm(WhiteSurfArea_area_global ~ FHR_str*sex + age + site + mean_eulnum_cross, data=master_df_T1)) # yes, signif interaction
temp <- master_df_T1[(master_df_T1$FHR_str!='BP' & master_df_T1$sex_str=='male'),]
t.test(temp$WhiteSurfArea_area_global~temp$FHR_str) # yes, signif t-test
temp <- master_df_T1[(master_df_T1$FHR_str!='BP' & master_df_T1$sex_str=='female'),]
t.test(temp$WhiteSurfArea_area_global~temp$FHR_str) # yes, NOT a signif t-test

# FHR-BP females should show larger BRAIN VOLUME than controls
summary(lm(BrainSegVolNotVent ~ FHR_str*sex + age + site + mean_eulnum_cross, data=master_df_T1)) # yes, signif interaction
temp <- master_df_T1[(master_df_T1$FHR_str!='SZ' & master_df_T1$sex_str=='male'),]
t.test(temp$BrainSegVolNotVent~temp$FHR_str) # yes, NOT a signif t-test
temp <- master_df_T1[(master_df_T1$FHR_str!='SZ' & master_df_T1$sex_str=='female'),]
t.test(temp$BrainSegVolNotVent~temp$FHR_str) # yes but barely (p=0.047), signif t-test; in the paper I think it was weak too

# FHR-BP females should show larger CORTICAL VOLUME than controls
summary(lm(eTIV ~ FHR_str*sex + age + site + mean_eulnum_cross, data=master_df_T1)) # yes, signif interaction
temp <- master_df_T1[(master_df_T1$FHR_str!='SZ' & master_df_T1$sex_str=='male'),]
t.test(temp$eTIV~temp$FHR_str) # yes, NOT a signif t-test
temp <- master_df_T1[(master_df_T1$FHR_str!='SZ' & master_df_T1$sex_str=='female'),]
t.test(temp$eTIV~temp$FHR_str) # not really... a marginal effect (p=0.07); in the paper I think it was weak too

# should be no evidence for interaction effects for cortical thickness
summary(lm(MeanThickness_thickness_global ~ FHR_str*sex + age + site + mean_eulnum_cross, data=master_df_T1)) # yes, signif interaction
# yep, same finding, VERY much not signif...

# ==== visualizations ====

colnames(master_df)
ggplot(master_df, aes(x = time, y = eTIV, fill = FHR_str)) +
  geom_violin(trim = FALSE, position = position_dodge(0.8)) +   
  scale_fill_brewer(palette = "Set3") +  
  labs(x = "Time", y = "eTIV", fill = "Group", color = "Group") +   
  theme_minimal()

# ==== analysis 1 ====

# read in data
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
master_df <- read.csv('master_df_TEST.csv')


# !!!! Create false variable for testing statistical model !!!!

master_df$FHR_SHUFFLE <- sample(master_df$FHR)
table(master_df$FHR_SHUFFLE, master_df$FHR)

# analysis
# clean a little
master_df$sex_str[master_df$sex_str=="male"] <- 'M'
master_df$sex_str[master_df$sex_str=="female"] <- 'F'
# make sure variables are correct classes
master_df$sex <- as.factor(master_df$sex)
master_df$FHR <- as.factor(master_df$FHR)
master_df$FHR_SHUFFLE <- as.factor(master_df$FHR_SHUFFLE) # the false variable !
# change the reference to PBC (if desired)
master_df$FHR <- relevel(master_df$FHR, ref = "3")
master_df$FHR_SHUFFLE <- relevel(master_df$FHR_SHUFFLE, ref = "3") # the false variable !

# linear mixed effects model with random intercept
model <- lmer(BrainSegVolNotVent ~ time*FHR_SHUFFLE + age + sex_str + site +
                mean_eulnum_cross + (1 | sub_id), data = master_df)
summary(model)


# ==== end








