

# outlier

library(readODS)
library(readxl)
library(tidyr)
library(ggplot2)
library(dplyr)

# ==== reading in data ====

# read in your list of subs, only keep keepers
# this is to compare with data you read in, does it all match?
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
sublist <- read_ods('sub_list.ods')
sublist <- sublist[sublist$image_incl==1,]
sublist <- sublist[sublist$qc_incl==1,]
sublist <- as.data.frame(sublist[,c(1,2)])
rownames(sublist) <- NULL

# is there longitudinal data?
sublist$long <- ifelse(sublist$sub_id %in% sublist$sub_id[duplicated(sublist$sub_id)],'both',sublist$time) 
# remove duplicates to look at base data
sublist_base <- sublist[!duplicated(sublist$sub_id), ]
rownames(sublist_base) <- NULL

# get demo variables to control for
# match to sublist_base
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
demos <- as.data.frame(read_xlsx('VIA11-15_longitudinal_demograph_clinic_wide.xlsx'))
demos <- demos[demos$famlbnr %in% sublist_base$sub_id,]
demos$sub_id <- demos$famlbnr
all(demos$sub_id == sublist_base$sub_id)

# filter demos out to match availability of imaging data
# important variables to clean:
imp_cols <- c('via11_mri_age','via15_mri_age','via11_mri_site','via15_mri_site','sex_code')
# first cleaning loop:
for(imp_col in imp_cols){
  demos[[imp_col]] <- ifelse(demos[[imp_col]]=="NA",NA,demos[[imp_col]])
  if(imp_col == 'via11_mri_age' | imp_col == 'via15_mri_age'){
    demos[[imp_col]] <- as.numeric(demos[[imp_col]])
  }}
# second cleaning loop:
if(all(sublist_base$sub_id == demos$famlbnr)){
  for(sub in 1:dim(demos)[1]){
    if(sublist_base$long[sub]=='both'){
      # do nothing
    } else if(sublist_base$long[sub]=='ses01'){
        demos$via15_mri_age[sub] <- NA
        demos$via15_mri_site[sub] <- NA
    } else if(sublist_base$long[sub]=='ses02'){
        demos$via11_mri_age[sub] <- NA
        demos$via11_mri_site[sub] <- NA
    }} 
  print('done cleaning demos')}

# create a few variables for base outliers, method 2 below
# create mixed age variable
demos$mixed_age[sublist_base$long=='ses01'] <- demos$via11_mri_age[sublist_base$long=='ses01']
demos$mixed_age[sublist_base$long=='ses02'] <- demos$via15_mri_age[sublist_base$long=='ses02']
demos$mixed_age[sublist_base$long=='both'] <- (demos$via11_mri_age[sublist_base$long=='both'] + demos$via15_mri_age[sublist_base$long=='both'])/2

# create a mixed site variable
demos$mixed_site[sublist_base$long=='ses01'] <- paste(demos$via11_mri_site[sublist_base$long=='ses01'],'_NA',sep='')
demos$mixed_site[sublist_base$long=='ses02'] <- paste('NA_',demos$via15_mri_site[sublist_base$long=='ses02'],sep='')
demos$mixed_site[sublist_base$long=='both'] <- paste(demos$via11_mri_site[sublist_base$long=='both'],'_',demos$via15_mri_site[sublist_base$long=='both'],sep='')

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

# clean generated data tables (e.g., remove duplicate subs if they exist)
clean_up <- function(df) {
  # remove duplicated subs
  to_check <- c(110, 168, 203, 429, 442, 511)
  for(num in to_check){
    if(sum(df$famlbnr==num)==2){
      df <- df[!(df$famlbnr == num & !duplicated(df$famlbnr)), ]
    }
  }
  df <- df[!df$famlbnr %in% c(176,232,296),] # remove subjects which failed cortical QC
  df <- df[order(df$famlbnr),]
  rownames(df) <- NULL
  return(df)
}

# check if dfs match
df_check <- function(df1,df2) {
  if(dim(df1)[1] == dim(df2)[1]){
    for(i in 1:dim(df1)[1]){
      if(df1$sub_id[i] != df2$famlbnr[i]){
        return(paste('problem with row ',i,'sublist: ',df1$sub_id[i],', data: ',df2$famlbnr[i],sep=""))
      }
    }
    return('finished, no problems')
  }
  return('not matching sizes')
}

# get outliers
get_outs <- function(varib) {
  avg <- mean(varib, na.rm = TRUE)
  std <- sd(varib, na.rm = TRUE)
  
  return(as.numeric((varib > avg + 2*std | varib < avg - 2*std)))
}

# ==== FOR BASE, method 1: separate models ====

# set directory for data
setwd('/mnt/projects/VIA_BIDS/BIDS_VIA11-15_derivatives/stat_tables/ANALYSIS_CORTICAL/')
# loop through files, initialize df to keep outliers
outlier_flags <- data.frame(sub_id = sublist_base$sub_id)
for(f in list.files(pattern='_base_')){
  # decide whether or not to read in file
  # if yes,
  if(read_file_or_not(f)){
    # read in file
    measures <- read.csv(f)
    
    # clean it
    measures_clean <- clean_up(measures)
    
    # double check it looks good compared to your sublist, defined above
    print('comparing data to sublist_base...')
    print(df_check(sublist_base,measures_clean))
    print('comparing data to demos...')
    print(df_check(demos,measures_clean))
    readline(prompt = "press enter to continue")
    
    # and selected the columns to look at for outliers
    selected_cols <- select_columns_console(measures_clean)
    print(selected_cols)
    
    for(m in selected_cols){
      
      # regress out age, sex, and site
      
      # must be repeated for subjects with ses01, ses02, and longitudinal data separately
      subs_ses01 <- sublist_base$sub_id[sublist_base$long=='ses01']
      subs_ses02 <- sublist_base$sub_id[sublist_base$long=='ses02']
      subs_both <- sublist_base$sub_id[sublist_base$long=='both']
      
      vari_ses01 <- measures_clean[[m]][sublist_base$long=='ses01']
      vari_ses02 <- measures_clean[[m]][sublist_base$long=='ses02']
      vari_both <- measures_clean[[m]][sublist_base$long=='both']
    
      df_ses01 <- cbind(vari_ses01, demos[sublist_base$long=='ses01',c('via11_mri_age','via11_mri_site',
                                     'sex_code')])
      
      df_ses02 <- cbind(vari_ses02, demos[sublist_base$long=='ses02',c('via15_mri_age','via15_mri_site',
                                     'sex_code')])
      
      df_both <- cbind(vari_both, demos[sublist_base$long=='both',c('via11_mri_age','via15_mri_age',
                              'via11_mri_site','via15_mri_site',
                              'sex_code')])
      
      # run models
      model_ses01 <- lm(vari_ses01 ~ via11_mri_age + via11_mri_site + sex_code, data = df_ses01)
      
      model_ses02 <- lm(vari_ses02 ~ via15_mri_age + via15_mri_site + sex_code, data = df_ses02)
      
      model_both <- lm(vari_both ~ via11_mri_age + via15_mri_age + via11_mri_site + 
                          via15_mri_site + sex_code, data = df_both)
      
      vari_resid_ses01 <- residuals(model_ses01)
      vari_resid_ses02 <- residuals(model_ses02)
      vari_resid_both <- residuals(model_both)
      
      subs_ses01 <- as.data.frame(cbind(subs_ses01, get_outs(vari_resid_ses01)))
      subs_ses02 <- as.data.frame(cbind(subs_ses02, get_outs(vari_resid_ses02)))
      subs_both <- as.data.frame(cbind(subs_both, get_outs(vari_resid_both)))
      
      # add ses01 outliers
      match_idx <- match(outlier_flags$sub_id, subs_ses01$subs_ses01)
      non_na <- !is.na(match_idx) & !is.na(subs_ses01$V2[match_idx])
      outlier_flags[[m]][non_na] <- subs_ses01$V2[match_idx][non_na]
 
      # add ses02 outliers
      match_idx <- match(outlier_flags$sub_id, subs_ses02$subs_ses02)
      non_na <- !is.na(match_idx) & !is.na(subs_ses02$V2[match_idx])
      outlier_flags[[m]][non_na] <- subs_ses02$V2[match_idx][non_na]
      
      # add both outliers
      match_idx <- match(outlier_flags$sub_id, subs_both$subs_both)
      non_na <- !is.na(match_idx) & !is.na(subs_both$V2[match_idx])
      outlier_flags[[m]][non_na] <- subs_both$V2[match_idx][non_na]
    }
  }
}

# ==== FOR BASE, method 2: mixed age and site vars - CORRECT ==== 

# set directory for data
setwd('/mnt/projects/VIA_longitudin/adam/data/ANALYSIS_CORTICAL/')
# initialize df to keep outliers
outlier_flags <- data.frame(sub_id = sublist_base$sub_id)
# loop through files
for(f in list.files(pattern='_base_')){
  # decide whether or not to read in file
  # if yes,
  if(read_file_or_not(f)){
    # read in file
    measures <- read.csv(f)
    
    # clean it
    measures_clean <- clean_up(measures)
    
    # double check it looks good compared to your sublist, defined above
    print('comparing data to sublist_base...')
    print(df_check(sublist_base,measures_clean))
    print('comparing data to demos...')
    demos$sub_id <- demos$famlbnr
    print(df_check(demos,measures_clean))
    readline(prompt = "press enter to continue")
    
    # and selected the columns to look at for outliers
    # desired columns should be:
    # 15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48
    selected_cols <- select_columns_console(measures_clean)
    print(selected_cols)
    
    for(m in selected_cols){
      
      # regress out age, sex, and site
      
      # make analysis df
      vari <- measures_clean[[m]]
      df <- cbind(vari, demos[,c('mixed_age','mixed_site','sex_code')])
      # run model
      model <- lm(vari ~ mixed_age + mixed_site + sex_code, data = df)
      # get residual
      vari_resid <- rstandard(model) # standardized residuals
      # add to outlier df
      outlier_flags[[m]] <- get_outs(vari_resid)
    }
  }
}

# ==== post base ====

outlier_flags # generated above

outlier_flags_mixed_vars <- outlier_flags

outlier_means <- rowSums(outlier_flags[, -1])

hist(outlier_means,
     breaks=20,
     xlab = 'Number of outliers',
     main = 'Outlier distribution',
     col = 'lightblue')

# create an outlier df with:
# - sub_id
# - number of total outliers, per sub
if(all(sublist_base$sub_id==outlier_flags$sub_id)){
  print('subs are consistent; outlier_df created')
  outlier_df <- as.data.frame(cbind(outlier_flags$sub_id,outlier_means))
  colnames(outlier_df) <- c('sub_id','outlier_means')
}

# determine subjects with outliers
i <- 0
for(r in 1:length(outlier_means)){
  if(outlier_df$outlier_means[r]>0){
    print(paste('Subject ',outlier_df$sub_id[r],' has ',outlier_df$outlier_means[r],' outliers',sep=''))
    i <- i + 1
  }
}
print(paste('Subjects with at least 1 outlier: ',i,sep=""))

## temporary workspace here. Save outlier data from alternative approaches

#outlier_flags_mixed_vars <- outlier_flags
#outlier_means_mixed_vars <- outlier_means

r_value <- cor(outlier_means_mixed_vars, outlier_means_sep_models)
plot(outlier_means_mixed_vars,outlier_means_sep_models,
     xlab = 'Mixed age and site variables',
     ylab = 'Separate models',
     col = ifelse(outlier_df_mixed_vars$highlight_4 == 1, 'red', 'black'),) # change highlight here
abline(lm(outlier_means_sep_models ~ outlier_means_mixed_vars), col = 'blue', lwd = 1)
text(x = 5, y = 65, 
     labels = paste0('r = ', round(r_value, 2)), 
     pos = 4, cex = 1)
text(x = 5, y = 45, 
     #labels = 'red = at least 5% (10.2) of parcels are outliers', # for highlight_2
     #labels = 'red = at least 5% (10.2) of parcels are outliers\nAND any rating >1', # for highlight_3
     labels = 'red = at least 5% (10.2) of parcels are outliers\nAND any rating >1\nAND data looks bad', # for highlight_4
     pos = 4, cex = 1, col = 'red')

# determine top X%
outlier_df_mixed_vars <- as.data.frame(cbind(outlier_flags_mixed_vars$sub_id,
                               rowSums(outlier_flags_mixed_vars[, -1])))
colnames(outlier_df_mixed_vars) <- c('sub_ids','outliers')

outlier_df_mixed_vars$highlight <- 0
outlier_df_mixed_vars$highlight[which(outlier_df_mixed_vars$outliers >= quantile(outlier_df_mixed_vars$outliers, 0.90))] <- 1
outlier_df_mixed_vars$sub_ids[which(outlier_df_mixed_vars$outliers >= quantile(outlier_df_mixed_vars$outliers, 0.90))]

# determine if X% are outliers
outlier_df_mixed_vars$highlight_2 <- 0
outlier_df_mixed_vars$highlight_2[which(outlier_df_mixed_vars$outliers >= 204*0.05)] <- 1
outlier_df_mixed_vars$sub_ids[which(outlier_df_mixed_vars$outliers >= 204*0.05)]

setwd('~/Desktop/VIA11_VIA15_project/')
orig_ratings <- read.csv('freesurfer_paths_quality.csv')
orig_ratings <- orig_ratings[,c('famlbnr','time','Rating_USE')]
orig_ratings <- pivot_wider(orig_ratings, 
                       names_from = time, 
                       values_from = Rating_USE)
orig_ratings <- orig_ratings[orig_ratings$famlbnr %in% outlier_df_mixed_vars$sub_ids,]
orig_ratings <- orig_ratings[order(orig_ratings$famlbnr),]
all(orig_ratings$famlbnr==outlier_df_mixed_vars$sub_ids)

setwd('/mnt/projects/VIA_longitudin/adam/qc/')
cortical_ratings <- read_ods('Cortical_QC.ods')
cortical_ratings <- cortical_ratings[!is.na(cortical_ratings$Subject),]
cortical_ratings$Subject <- gsub("sub-via|_base", "", cortical_ratings$Subject)
cortical_ratings <- cortical_ratings[!grepl("_OLD", cortical_ratings$Subject),]
cortical_ratings$Subject <- as.numeric(cortical_ratings$Subject)
cortical_ratings <- cortical_ratings[cortical_ratings$Subject %in% outlier_df_mixed_vars$sub_ids,]
cortical_ratings <- cortical_ratings[order(cortical_ratings$Subject),]
all(cortical_ratings$Subject==outlier_df_mixed_vars$sub_ids)

outlier_master <- as.data.frame(cbind(outlier_df_mixed_vars$sub_ids,
                        outlier_df_mixed_vars$highlight_2,
                        orig_ratings$VIA11,
                        orig_ratings$VIA15,
                        cortical_ratings$Internal_QC,
                        cortical_ratings$External_QC))
colnames(outlier_master) <- c('sub_id','fivepct_outliers','ses01_rating',
                              'ses02_rating','base_internal','base_external')

outlier_master <- outlier_master[outlier_master$fivepct_outliers==1,]
outlier_master <- outlier_master[is.na(outlier_master$ses01_rating) | outlier_master$ses01_rating>1 |
                                 is.na(outlier_master$ses02_rating) | outlier_master$ses02_rating>1 |
                                 outlier_master$base_internal>1 | outlier_master$base_external >1,]
                
outlier_master$keep <- 1            
for(i in 1:dim(outlier_master)[1]){
  temp <- c(outlier_master$fivepct_outliers[i],
            outlier_master$ses01_rating[i],
            outlier_master$ses02_rating[i],
            outlier_master$base_internal[i],
            outlier_master$base_external[i])
  if(all(is.na(temp) | temp == 1)){
    outlier_master$keep[i] <- 0
  }
}    

outlier_master <- outlier_master[outlier_master$keep==1,]

outlier_df_mixed_vars$highlight_3 <- 0
outlier_df_mixed_vars$highlight_3[outlier_df_mixed_vars$sub_ids %in% outlier_master$sub_id] <- 1

# visual inspection, bad subjects determined and listed below

outlier_df_mixed_vars$highlight_4 <- 0
outlier_df_mixed_vars$highlight_4[outlier_df_mixed_vars$sub_ids %in% c(5, 14, 135, 220, 248, 285, 297, 429, 466, 513)] <- 1

## end workspace


# ==== pre longitudinal ====

# purpose here is to create violin plots to show change in residuals from 
# VIA11 to VIA15 - how much variability in slopes is there? What will our 
# best approach be to determine outliers?
# Since we're dealing >200 ROIs, to simplify we will plot whole brain measures
# only

setwd('/mnt/projects/VIA_BIDS/BIDS_VIA11-15_derivatives/stat_tables/ANALYSIS_CORTICAL/')
# determine who has longitudinal data
sublist_longi <- sublist[duplicated(sublist$sub_id),]
rownames(sublist_longi) <- NULL
# match demo subs to these subs
demos_longi <- demos[demos$famlbnr %in% sublist_longi$sub_id,]
demos_longi$sub_id <- demos_longi$famlbnr
for(f in list.files(pattern='_long_')){
  # decide whether or not to read in file
  # if yes,
  if(read_file_or_not(f)){
    # read in file
    measures <- read.csv(f)
    
    # clean it
    measures_clean <- clean_up(measures)
    
    # longitudinal data only
    measures_clean_longi <- measures_clean[duplicated(measures_clean$famlbnr) | duplicated(measures_clean$famlbnr, fromLast = TRUE), ]
    rownames(measures_clean_longi) <- NULL
    measures_clean_longi_1 <- measures_clean_longi[measures_clean_longi$session_id=='ses01',]
    measures_clean_longi_2 <- measures_clean_longi[measures_clean_longi$session_id=='ses02',]
    
    # double check it looks good compared to your sublist, defined above
    print('comparing ses01 and ses02 data to sublist_base...')
    print(df_check(sublist_longi,measures_clean_longi_1))
    print(df_check(sublist_longi,measures_clean_longi_2))
    print('comparing ses01 and ses02 data to demos...')
    print(df_check(demos_longi,measures_clean_longi_1))
    print(df_check(demos_longi,measures_clean_longi_2))
    readline(prompt = "press enter to continue")
    
    # and select the columns to look at for outliers
    selected_cols <- select_columns_console(measures_clean_longi)
    print(selected_cols)
    
    for(m in selected_cols){
      
      # regress out age, sex, and site
      
      # make analysis dfs
      vari_1 <- measures_clean_longi_1[[m]]
      vari_2 <- measures_clean_longi_2[[m]]
      df_1 <- cbind(vari_1, demos_longi[,c('via11_mri_age','via11_mri_site','sex_code')])
      df_2 <- cbind(vari_2, demos_longi[,c('via15_mri_age','via15_mri_site','sex_code')])
      # run models
      model_1 <- lm(vari_1 ~ via11_mri_age + via11_mri_site + sex_code, data = df_1)
      model_2 <- lm(vari_2 ~ via15_mri_age + via15_mri_site + sex_code, data = df_2)
      # get residual
      vari_1_resid <- residuals(model_1)
      vari_2_resid <- residuals(model_2)
      
      dif <- vari_2_resid - vari_1_resid
      outs <- get_outs(dif)
      
      # create violin plots and save to your directory
      df_plot <- as.data.frame(cbind(measures_clean_longi_1$famlbnr,
                                     vari_1_resid,
                                     vari_2_resid,
                                     outs))
      
      colnames(df_plot) <- c('sub_id','sess01','sess02','outs')
      
      df_plot <- df_plot %>%
        pivot_longer(cols = c(sess01, sess02), 
                     names_to = "Time", 
                     values_to = "Value")
      
      viol_plot <- ggplot(df_plot, aes(x = Time, y = Value)) +
        geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +  # Violin plot
        geom_jitter(aes(color = as.factor(sub_id)), width = 0.1, size = 2) +  # Individual points
        geom_line(aes(group = sub_id, color = as.factor(outs)), alpha = 0.5) +  # Connecting lines
        scale_color_manual(values = c("0" = "darkgrey", "1" = "red")) +
        theme_minimal() +
        labs(x = "Time", y = paste(m," residuals adjusting for age, site, sex",sep=""), title = m) +
        theme(legend.position = "none")
      
      output_dir <- '/mnt/projects/VIA_longitudin/adam/qc/violin_plots/'
      output_file <- paste0(output_dir, m,"_viol.png")
      ggsave(filename = output_file, plot = viol_plot, width = 8, height = 6, dpi = 300, bg = 'white')
      print(paste("Plot saved to:", output_file))
      
    }
  }
}


# ==== FOR LONGITUDINAL ====

# set directory for data
setwd('/mnt/projects/VIA_BIDS/BIDS_VIA11-15_derivatives/stat_tables/ANALYSIS_CORTICAL/')
# determine who has longitudinal data
sublist_longi <- sublist[duplicated(sublist$sub_id),]
rownames(sublist_longi) <- NULL
# match demo subs to these subs
demos_longi <- demos[demos$famlbnr %in% sublist_longi$sub_id,]
demos_longi$sub_id <- demos_longi$famlbnr
# loop through files, initialize df to keep outliers
outlier_flags_longi <- data.frame(sub_id = sublist_longi$sub_id)
for(f in list.files(pattern='_long_')){
  # decide whether or not to read in file
  # if yes,
  if(read_file_or_not(f)){
    # read in file
    measures <- read.csv(f)
    
    # clean it
    measures_clean <- clean_up(measures)
    
    # longitudinal data only
    measures_clean_longi <- measures_clean[duplicated(measures_clean$famlbnr) | duplicated(measures_clean$famlbnr, fromLast = TRUE), ]
    rownames(measures_clean_longi) <- NULL
    measures_clean_longi_1 <- measures_clean_longi[measures_clean_longi$session_id=='ses01',]
    measures_clean_longi_2 <- measures_clean_longi[measures_clean_longi$session_id=='ses02',]
    
    # double check it looks good compared to your sublist, defined above
    print('comparing ses01 and ses02 data to sublist_base...')
    print(df_check(sublist_longi,measures_clean_longi_1))
    print(df_check(sublist_longi,measures_clean_longi_2))
    print('comparing ses01 and ses02 data to demos...')
    print(df_check(demos_longi,measures_clean_longi_1))
    print(df_check(demos_longi,measures_clean_longi_2))
    readline(prompt = "press enter to continue")
    
    # and select the columns to look at for outliers
    selected_cols <- select_columns_console(measures_clean_longi)
    print(selected_cols)
    
    for(m in selected_cols){
      
      # regress out age, sex, and site
      
      # make analysis dfs
      vari_1 <- measures_clean_longi_1[[m]]
      vari_2 <- measures_clean_longi_2[[m]]
      df_1 <- cbind(vari_1, demos_longi[,c('via11_mri_age','via11_mri_site','sex_code')])
      df_2 <- cbind(vari_2, demos_longi[,c('via15_mri_age','via15_mri_site','sex_code')])
      # run models
      model_1 <- lm(vari_1 ~ via11_mri_age + via11_mri_site + sex_code, data = df_1)
      model_2 <- lm(vari_2 ~ via15_mri_age + via15_mri_site + sex_code, data = df_2)
      # get residual
      vari_1_resid <- rstandard(model_1) # standardized residuals
      vari_2_resid <- rstandard(model_2) # standardized residuals
      # get difference
      dif <- (vari_2_resid - vari_1_resid) # /((vari_2_resid + vari_1_resid)/2) Kathrine and I decided NOT to divide by the mean
      # add to outlier df
      outlier_flags_longi[[m]] <- get_outs(dif)
    }
  }
}

# ==== post longitudinal ====

outlier_flags_longi # generated above

#outlier_flags_longi_1 <- outlier_flags_longi ... without dividing by mean
#outlier_flags_longi_2 <- outlier_flags_longi ... WITH dividing by mean

outlier_means_longi <- rowSums(outlier_flags_longi[, -1])

hist(outlier_means_longi,
     breaks=10,
     xlab = 'Number of outliers',
     main = 'Outlier distribution',
     col = 'lightblue')

# create an outlier df with:
# - sub_id
# - number of total outliers, per sub
if(all(sublist_longi$sub_id==outlier_flags_longi$sub_id)){
  print('subs are consistent; outlier_df_longi created')
  outlier_df_longi <- as.data.frame(cbind(outlier_flags_longi$sub_id,outlier_means_longi))
  colnames(outlier_df_longi) <- c('sub_id','outlier_means')
}

# determine subjects with outliers
i <- 0
for(r in 1:length(outlier_means_longi)){
  if(outlier_df_longi$outlier_means[r]>0){
    print(paste('Subject ',outlier_df_longi$sub_id[r],' has ',outlier_df_longi$outlier_means[r],' outliers',sep=''))
    i <- i + 1
  }
}
print(paste('Subjects with at least 1 outlier: ',i,sep=""))

## temporary workspace here. Save outlier data from alternative approaches

# just thickness
colnames(outlier_flags_longi)
outlier_flags_longi_area <- outlier_flags_longi[,c(1:69)]
outlier_flags_longi_thickness <- outlier_flags_longi[,c(1,70:137)]
outlier_flags_longi_volume <- outlier_flags_longi[,c(1,138:205)]

# subs to check for thickness: 244 ses01, 449 ses01
rowSums(outlier_flags_longi_thickness[,-1])[outlier_flags_longi_thickness$sub_id==244] # only 2
rowSums(outlier_flags_longi_thickness[,-1])[outlier_flags_longi_thickness$sub_id==244] # 


# determine top X%
outlier_df_longi$highlight <- 0
outlier_df_longi$highlight[which(outlier_df_longi$outlier_means >= quantile(outlier_df_longi$outlier_means, 0.90))] <- 1
outlier_df_longi$sub_id[which(outlier_df_longi$outlier_means >= quantile(outlier_df_longi$outlier_means, 0.90))]

# determine if X% are outliers
outlier_df_longi$highlight_2 <- 0
outlier_df_longi$highlight_2[which(outlier_df_longi$outlier_means >= 204*0.05)] <- 1
outlier_df_longi$sub_id[which(outlier_df_longi$outlier_means >= 204*0.05)]

setwd('~/Desktop/VIA11_VIA15_project/')
orig_ratings <- read.csv('freesurfer_paths_quality.csv')
orig_ratings <- orig_ratings[,c('famlbnr','time','Rating_USE')]
orig_ratings <- pivot_wider(orig_ratings, 
                            names_from = time, 
                            values_from = Rating_USE)
orig_ratings <- orig_ratings[orig_ratings$famlbnr %in% outlier_df_longi$sub_id,]
orig_ratings <- orig_ratings[order(orig_ratings$famlbnr),]
all(orig_ratings$famlbnr==outlier_df_longi$sub_id)

setwd('/mnt/projects/VIA_longitudin/adam/qc/')
cortical_ratings <- read_ods('Cortical_QC.ods')
cortical_ratings <- cortical_ratings[!is.na(cortical_ratings$Subject),]
cortical_ratings$Subject <- gsub("sub-via|_base", "", cortical_ratings$Subject)
cortical_ratings <- cortical_ratings[!grepl("_OLD", cortical_ratings$Subject),]
cortical_ratings$Subject <- as.numeric(cortical_ratings$Subject)
cortical_ratings <- cortical_ratings[cortical_ratings$Subject %in% outlier_df_longi$sub_id,]
cortical_ratings <- cortical_ratings[order(cortical_ratings$Subject),]
all(cortical_ratings$Subject==outlier_df_longi$sub_id)

outlier_master <- as.data.frame(cbind(outlier_df_longi$sub_id,
                                      outlier_df_longi$highlight_2,
                                      orig_ratings$VIA11,
                                      orig_ratings$VIA15,
                                      cortical_ratings$Internal_QC,
                                      cortical_ratings$External_QC))
colnames(outlier_master) <- c('sub_id','fivepct_outliers','ses01_rating',
                              'ses02_rating','base_internal','base_external')

outlier_master <- outlier_master[outlier_master$fivepct_outliers==1,]
outlier_master <- outlier_master[is.na(outlier_master$ses01_rating) | outlier_master$ses01_rating>1 |
                                   is.na(outlier_master$ses02_rating) | outlier_master$ses02_rating>1 |
                                   outlier_master$base_internal>1 | outlier_master$base_external >1,]

outlier_master$keep <- 1            
for(i in 1:dim(outlier_master)[1]){
  temp <- c(outlier_master$fivepct_outliers[i],
            outlier_master$ses01_rating[i],
            outlier_master$ses02_rating[i],
            outlier_master$base_internal[i],
            outlier_master$base_external[i])
  if(all(is.na(temp) | temp == 1)){
    outlier_master$keep[i] <- 0
  }
}    

outlier_master <- outlier_master[outlier_master$keep==1,]

outlier_df_longi$highlight_3 <- 0
outlier_df_longi$highlight_3[outlier_df_longi$sub_id %in% outlier_master$sub_id] <- 1

# visual inspection, bad subjects determined and listed below

outlier_df_longi$highlight_4 <- 0
outlier_df_longi$highlight_4[outlier_df_longi$sub_id %in% c()] <- 1

## end workspace


# ==== end









# ==== euler ====

# read in euler data
setwd('/mnt/projects/VIA_longitudin/adam/tasks/analysis_1/')
fs_qa_dat <- read.csv('VIA11-15_fs741_long_qa_measures_restructured.csv')
fs_qa_dat <- fs_qa_dat[fs_qa_dat$famlbnr %in% sublist$sub_id,]
fs_qa_dat <- fs_qa_dat[order(fs_qa_dat$famlbnr),]
length(unique(fs_qa_dat$famlbnr))==length(unique(demos$famlbnr))

# as of now, excluding for orig image quality and my personal qc
# NOT excluding for outliers + enhanced inspection of internal/external qc
# SO.... we NEED TO EXCLUDE THE STATISTICAL OUTLIERS
# subs which are outliers to be excluded:
# - 5 (ses01), 135 (ses01), 285 (ses01), 449 (ses01)

euler_dat <- merge(fs_qa_dat, demos, by='famlbnr')

subs_to_remove <- c(5,135,285,449) # 449 was the only longitudinal outlier - only ses01 was bad; it's been removed
euler_dat <- euler_dat[!(euler_dat$session_id == "ses01" & euler_dat$famlbnr %in% subs_to_remove),]

# regress out age, site, and sex, as before
# select only necessary variables for this
# and also correctly convert to long
euler_dat <- euler_dat[,c(1,8,19,49,51,55,56,57)]
colnames(euler_dat) <- c('sub_ID','time','eul_num','VIA11_site','VIA15_site',
                         'sex','VIA11_age','VIA15_age')
euler_dat <- euler_dat %>%
  mutate(mri_age = case_when(
    time == "ses01" ~ VIA11_age,
    time == "ses02" ~ VIA15_age
  ), mri_site = case_when(
    time == "ses01" ~ VIA11_site,
    time == "ses02" ~ VIA15_site
  )) %>%
  select(-VIA11_age, -VIA15_age, -VIA11_site, -VIA15_site)

euler_dat$sex <- as.factor(euler_dat$sex)
euler_dat$mri_site <- as.factor(euler_dat$mri_site)

# loop through time - find residuals for euler number separately by time
for(t in c('ses01','ses02')){
  print(paste('TIME = ',t,sep=''))
  euler_dat_temp <- euler_dat[euler_dat$time==t,]
  # set up model
  model_eul <- lm(eul_num ~ mri_age + mri_site + sex, data = euler_dat_temp)
  # get residual
  euler_dat_temp$eul_resid <- rstandard(model_eul)
  # get outliers -- THREE SDs, different than before
  cutoff <- function(data, x) {
    mean(data, na.rm = TRUE) + c(-1, 1) * x * sd(data, na.rm = TRUE) }
  x <- 3 # x is number of SDs
  low <- cutoff(euler_dat_temp$eul_resid,x)[1]
  high <- cutoff(euler_dat_temp$eul_resid,x)[2]
  euler_dat_temp$eul_outs <- as.numeric(euler_dat_temp$eul_resid < low | euler_dat_temp$eul_resid > high)
  # look at the outliers (what's their euler number? ... what's their sub ID and ses #?)
  print(euler_dat_temp$eul_num[euler_dat_temp$eul_outs==1])
  print(euler_dat_temp$sub_ID[euler_dat_temp$eul_outs==1])
  print(euler_dat_temp$time[euler_dat_temp$eul_outs==1])
}










