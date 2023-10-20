#############################
### 0. Create Directories ###
#############################

# Create data directory
data_dir <- dir.create("data")

# Designate subdirectories
data_subdirectories <- c("a_raw_data",
                         "b_intermediate_data",
                         "c_submodel_data",
                         "d_suitability_data",
                         "e_rank_data",
                         "f_sensitivity_data",
                         "g_uncertainty_data",
                         "zz_miscellaneous")

# Create sub-directories within data directory
for (i in 1:length(data_subdirectories)){
  subdirectories <- dir.create(paste0("data/", data_subdirectories[i]))
}

#####################################

# Create code directory
code_dir <- dir.create("code")

#####################################

# Create figure directory
figure_dir <- dir.create("figure")

#####################################
#####################################

# Delete directory (if necessary)
### ***Note: change directory path for desired directory
#unlink("data/a_raw_data", recursive = T)
