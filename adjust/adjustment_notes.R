## stashing this here to do age adjustments of the overall SCRs for 1-9 year olds

##########################################################

#setup data frames for plots and export
sp_bins_df <- sp_bins %>%
  mutate(age = as.factor(row.names(sp_bins)))

age_bins_mid_df <- as.data.frame(age_bins_mid) %>%
  filter(age_bins_mid != 9.5) %>%
  mutate(age = as.factor(sp_bins_df$age))

age_sp_bins_df <- left_join(sp_bins_df, age_bins_mid_df, by = "age") %>%
  mutate(age = as.numeric(age))

# age weights obtained from Tropical Data###################

age	<- as.numeric(c(1,2,3,4,5,6,7,8,9))

weights <- as.numeric(c(0.10691605, 0.10691605, 0.10691605, 0.10691605, 
                        0.11446716, 0.11446716, 0.11446716, 0.11446716, 
                        0.11446716))

age_weights_drc <- data.frame(age, weights) %>%
  write_excel_csv(files$ageweights_drc_clean)

##########################################################

# adjust data by age weights
# age weight * med/med est
# age weight * lcl/lcl est 
# age weight * lcl/ucl est

age_adj_obs_data <- age_sp_bins_df %>%
  mutate(aa_med = med*age_weights_drc$weights, 
         aa_low95 = low_95*age_weights_drc$weights, 
         aa_high95 = high_95*age_weights_drc$weights)

write_excel_csv(age_adj_obs_data, files$drc1_ct694_aa_obsdf)

##########################################################