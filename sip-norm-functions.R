setwd("/Users/DVyshenska/R_wd")
require(dplyr)
require("mgsub")
require("reshape2")
require("ggplot2")

# FUNCTIONS
parse_bbmap_cov <- function(coverage_path, sp_meta){
  
  # get all files in the directory
  files <- list.files(coverage_path)
  
  # import first file to create a draft of an output dataframe
  df <- read.csv(paste(coverage_path, files[1], sep = ""), sep = "\t", check.names = F,
                 stringsAsFactors = F)
  full_df <- df[, c("#ID", "Avg_fold")] # strip all except id and coverage info
  
  # generate new column names for the output dataframe
  new_colnames <- c("ID", mgsub(files, c("^.*_", "[.txt]"), c("", ""))) 
  colnames(full_df) <- new_colnames[c(1,2)] # give new column names to the first input
  
  # looping through the rest of the files
  for(i in 2:length(files)){
    # get full file
    tempdf <- read.csv(paste(coverage_path, files[i], sep = ""), sep = "\t", check.names = F,
                       stringsAsFactors = F)
    # strip not needed columns
    tempdf <- tempdf[,c("#ID", "Avg_fold")]
    # assign new names
    colnames(tempdf) <- c("ID", new_colnames[i+1])
    # merge with the rest
    full_df <- merge(full_df, tempdf, by = "ID")
  }
  
  # separating all data into:
  # only sequins data
  sequin_df <- full_df[full_df$ID %in% sp_meta[,1],]
  sequin_df <- sequin_df[order(sequin_df$ID),]
  rownames(sequin_df) <- sequin_df$ID
  sequin_df <- sequin_df[,-1]
  # the rest of the data (scaffolds/bacteria/etc.)
  bact_df <- full_df[!(full_df$ID %in% sp_meta[,1]),]
  bact_df <- bact_df[order(bact_df$ID),]
  rownames(bact_df) <- bact_df$ID
  bact_df <- bact_df[,-1]
  return(list(mags=bact_df, spequins=sequin_df))
}
parse_bbmap_read <- function(raw_files_path, sp_meta){# impoaring and parsing coverage data from mags and spikins
  files <- list.files(raw_files_path)
  
  # import first file to create a draft of an output dataframe
  df <- read.csv(paste(raw_files_path, files[1], sep = ""), sep = "\t",
                 stringsAsFactors = F, header = F)
  df <- df[-(nrow(df)),] # remove last row (contains unmapped reads #)
  #raw_colnames <- c("scaffold_id", "scaffold_length", 
  #                  "mapped_reads", "placed_reads") # see header explanation here: https://toolshed.g2.bx.psu.edu/repository/display_tool?repository_id=23812865694c57b7&render_repository_actions_for=tool_shed&tool_config=%2Fsrv%2Ftoolshed%2Fmain%2Fvar%2Fdata%2Frepos%2F001%2Frepo_1520%2Fsamtools_idxstats.xml&changeset_revision=88b8c2916784
  #colnames(df) <- raw_colnames
  
  full_df <- df[, c(1, 3)] # strip all except id and read info
  
  # generate new column names for the output dataframe
  new_colnames <- c("ID", mgsub(files, c("^.*_", "[.txt]"), c("", ""))) 
  colnames(full_df) <- new_colnames[c(1,2)] # give new column names to the first input
  
  # looping through the rest of the files
  for(i in 2:length(files)){
    # get full file
    tempdf <- read.csv(paste(raw_files_path, files[i], sep = ""), sep = "\t", 
                       header = F,
                       stringsAsFactors = F)
    # strip not needed columns
    tempdf <- tempdf[,c(1,3)]
    # assign new names
    colnames(tempdf) <- c("ID", new_colnames[i+1])
    # merge with the rest
    full_df <- merge(full_df, tempdf, by = "ID")
  }
  
  # separating all data into:
  # only sequins data
  sequin_df <- full_df[full_df$ID %in% sp_meta[,1],]
  sequin_df <- sequin_df[order(sequin_df$ID),]
  rownames(sequin_df) <- sequin_df$ID
  sequin_df <- sequin_df[,-1]
  # the rest of the data (scaffolds/bacteria/etc.)
  bact_df <- full_df[!(full_df$ID %in% sp_meta[,1]),]
  bact_df <- bact_df[order(bact_df$ID),]
  rownames(bact_df) <- bact_df$ID
  bact_df <- bact_df[,-1]
  return(list(mags=bact_df, spequins=sequin_df))
}

get_molecules <- function(input_mass, sp_meta){
  # cunction that calculates atommoles of dsDNA fragments & approximate number
  # of molecules in the given volume of DNA mix
  #
  # input: 
  #   input_mass: integer, total mass of spike-in DNA mix, in ng;
  #   sp_meta: table with first column spike-in ids, second - spike-in relative 
  # abundance and third column - length of the fragments in bp
  #
  # output:
  #   table with extra columns - calculated for each spike-in mass in grams, 
  # atommoles and molecules.
  
  g_inmass <- input_mass*10^-9
  sp_meta$g <- ((sp_meta[,2] * g_inmass)/ sum(sp_meta[,2]))
  sp_meta$moles <- sp_meta$g/((sp_meta[,3] * 607.4) + 157.9) # average molecular weight is from here https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html
  sp_meta$atommoles <- sp_meta$moles*10^18
  sp_meta$molecules <- (sp_meta$moles * 6.022*10^23)
  return(sp_meta)
}

get_scale_factors_sd <- function(sp_test_data, thresholdPerMean=100000){
  # This function takes spike-in data, calculates scaling factors for each 
  # sample based on spike-in data and returns calculated scaling factors.
  # it also filters out sets of spike-ins if at that concentration
  # not all spike-ins were detected (removed all spike-ins at that concentration)
  # it can also remove sets of spike-ins if their standar deviation 
  # as % from the mean was larger than thresholdPerMean value
  # this function also outputs graphs of spike-in behavior for each sample.
  
  # Function input:
  #   sp_test_data: data frame with spike-in data that we use for calculating
  # scaling factors. Column names are sample names, row names are spike-in names,
  # MUST contain additional column "rel_conc" that contains real relative
  # concentration numbers for each spike-in (taken from spike-in metadata file)
  #   thresholdPerMean: spike-in variation threshold (in %, ex. 20, 50, 100). If
  # for any given spike-in concentration measured values will deviate from the
  # mean more than the given threshold - this spike in will be disregarded from
  # the analysis (if located in the lower or upper ranges of the standard curve)
  # or will be rised as a warning (if located within a standard curve).
  
  # Function output is a list consisting from 4 data frames:
  #   scaleVal_df: dataframe that contains scaling factors for each sample
  # (intersept and slope). Used to normalize metagenomics data.
  #   pass_table: dataframe of logical values for each spike-in concentration in
  # each sample indicating if it passed the pre-set variation threshold or not.
  #   sdPercent_df: dataframe containing percent of standard deviation from 
  # the mean for each spike-in concentration for each sample.
  #   detection_pass: dataframe of logical values for each spike-in concentration
  # in each sample indicating if all spike-ins of each concentration were
  # detected as expected (FALSE - some or all were not detected, TRUE - all
  # expected spike-ins were detected)
  
  #
  #
  
  # checking if any spike-ins were not detected at all
  if(any(sp_test_data == 0)){
    cat("\nWarning! Some spike-ins have zero counts! Check your data!\n\n")
  }
  
  # getting means, standard deviations and percentages of sd from mean for 
  # each spike in. 
  mean_sp <- as.data.frame(sp_test_data %>% group_by(rel_conc) %>% summarise_all(mean))
  sd_sp <- abs(as.data.frame(sp_test_data %>% group_by(rel_conc) %>% summarise_all(sd)))
  per_sp <- (sd_sp[,-1]/mean_sp[,-1])*100
  
  # Getting a table of logical values for each spikein concentration based on 
  # the pre-set filtering threshold for percent of deviation.
  logic_per_sp <- per_sp <= thresholdPerMean
  rownames(logic_per_sp) <- sd_sp[,1]
  
  # setting percent deviation from mean values for undetected spikins to FALSE
  per_sp[mean_sp[,-1] == 0] <- NA

  # setting logical values FALSE where number of detected 
  # spikins less then expected
  sp_counts <- as.data.frame(sp_test_data %>% 
                               group_by(rel_conc) %>% 
                               summarise_all(function(x){sum(x!=0)}))
  logic_detected <- sp_counts[,-1] == as.vector(table(sp_test_data$rel_conc))
  
  # printing out warnings about underdetected samples and spike-in concentrations
  for(i in 1:ncol(logic_detected)){
    if(any(logic_detected[,i] == FALSE)){
      cat("\nIn sample ", colnames(logic_detected)[i], 
          "following concentrations were underdetected:",
          sp_counts$rel_conc[logic_detected[,i] == FALSE], "\n")
    }
  }
  rownames(logic_detected) <- sp_counts[,1]
  logic_per_sp[logic_detected == FALSE] <- FALSE

  # relocating "real concentration" column from spike-in dataframe to a vector
  sample_df <- sp_test_data[ , -which(names(sp_test_data) %in% c("rel_conc"))]
  rel_conc <- sp_test_data$rel_conc
  
  # creating new dataframe for scaling factors
  scale_val <- data.frame(matrix(nrow=2, ncol=ncol(sample_df)))
  colnames(scale_val) <- colnames(sample_df)
  rownames(scale_val) <- c("intersept", "slope")
  
  #looping through each sample separately
  for(i in 1:ncol(sample_df))  {
    sample <- sample_df[,i]
    logical <- logic_per_sp[,i]
    #cat(sample, "\n")
    #cat(logical, "\n")
    # getting the spike-in concentration that pass thresholds
    conc_to_use <- names(logical)[logical == TRUE]
    
    # calculating scaling factors with stripped spike-ins that are located on 
    # the fare ends of the concentration range and did not pass threshold
    
    x <- rel_conc[rel_conc %in% conc_to_use]
    y <- sample[rel_conc %in% conc_to_use]
    lm_res <- lm(y~x)
    scale_val[, i] <- as.vector(lm_res$coefficients)
    r_adj <-summary(lm_res)$adj.r.squared # also getting R-sq. adj for plot
    
    # # colorcoding spike-in concentrations
    dot_color <- vector()
    dot_color[!(rel_conc %in% conc_to_use)] <- "red" # did not pass threshold
    dot_color[rel_conc %in% conc_to_use] <- "green" # pass threshold
    
    # plot the data into the pdf file in working directory
    pdf(paste(colnames(sample_df)[i],  "_PerDevSD", thresholdPerMean,".pdf", sep = ""))
    plot(rel_conc, sample, xlab = "real relative concentration", ylab="measured value", 
         col=dot_color)
    # cat("Processing i = ", i, "\n")
    # cat("intercept is ", scale_val[1,i], " and slope is ", scale_val[2,i], "\n\n")
    abline(scale_val[1,i], scale_val[2,i])
    legend("topleft", legend=c("variation OUTSIDE the threshold", 
                               "variation WITHIN threshold",
                               paste("Adjusted R-squared:  ", round(r_adj, 4), sep = "")),
           col=c("red", "green", "white"), lty=1, cex=0.8)
    dev.off()
  }
  rownames(per_sp) <- rownames(logic_detected)

  return(list(scaleVal_df=scale_val, pass_table=logic_per_sp, 
              sdPercent_df=per_sp, detection_pass=logic_detected))
}

lm_normalize <- function(sampl_data, scaleVal_df) {
  # this function scales samples data using input scaling factors (intersept and
  # slope, calculated using function get_scale_factors())
  
  # Function input:
  #   sampl_data - a dataframe where row names are MAGs' names, column names
  # are samples' names; the data stored - data that needs to be scaled
  #   scaleVal_df: dataframe that contains scaling factors for each sample
  # (intersept- first row,  and slope - second row). Used to normalize 
  # metagenomics data. Column names must be the same as in sampl_data
  
  # Function output:
  # norm_data - contains scaled samples' data (columns are samples, rows are
  # MAGs)
  
  # making sure that column names for input data frames are identical
  if(all(colnames(sampl_data) == colnames(scaleVal_df))){
    # normalizing input sample data using intersept and slope from spikins:
    # creating empty dataframe to store normalized samples data
    norm_sample_data <- data.frame(matrix(nrow=nrow(sampl_data), ncol = ncol(sampl_data)))
    colnames(norm_sample_data) <- colnames(sampl_data)
    rownames(norm_sample_data) <- rownames(sampl_data)
    
    # normalizing each sample using calculated values and populating norm_sample_data
    for(i in 1:ncol(sampl_data)){
      norm_sample_data[,i] <- (sampl_data[,i] - scaleVal_df[1,i])/scaleVal_df[2,i]
    }
    return(norm_data = norm_sample_data)
  } else{
    cat("Check column names in input data frames! Exiting analysis.\n")
    break
  }
}

reorder_df <- function(input_table){
  # function tailored for sequins paper Mixes experiment
  # renames row names and column names, reorders table.
  # input must be a df with names of MAGs as rownames, sample names (1-6)
  # as column names. Must include
  # row "Sequins" (the row can have NA instead of numbers)
  plot_names <- c("Synechocystis sp. PCC 6803",
                  "Leptolyngbya sp. PCC 7376",
                  "Synechococcus elongatus PCC 7942",
                  "Nostoc punctiforme ATCC 29133",
                  "Sequins")
  mix_names <- c("Mix 1",	"Mix 2",	"Mix 3",
                 "Mix 4",	"Mix 5",	"Mix 6")
  
  rownames(input_table) <- mgsub(rownames(input_table), 
                                 c("NC_000911.1 Synechocystis sp. PCC 6803, complete sequence",
                                   "NC_019683.1 Leptolyngbya sp. PCC 7376, complete sequence",
                                   "NC_007604.1 Synechococcus elongatus PCC 7942, complete genome",
                                   "NC_010628.1 Nostoc punctiforme PCC 73102, complete sequence"),
                                 plot_names[-5])
  
  input_table <- input_table[match(plot_names, rownames(input_table)),] 
  

  
  colnames(input_table) <- mgsub(colnames(input_table), 
                                 c("2",	"1",	"4",	"3",	"6",	"5"), 
                                 c("Mix 1",	"Mix 2",	"Mix 3", "Mix 4",	"Mix 5",	"Mix 6"))
  reordered <- input_table[,match(mix_names, colnames(input_table))] 
  return(reordered)
}
plot_stackbar <- function(plot_df, file_name){
  # function creates a stacked barplot for mags as pdf file
  # input must be prepared by reorder_df. can handle NA row of "Sequins"
  
  if("Sequins" %in% rownames(plot_df)){
    plot_df <- plot_df[-which(rownames(plot_df) == "Sequins"),]
  }
  plot_df$mags <- rownames(plot_df)
  plot_df <- melt(plot_df, id.vars = "mags")
  plot_df$mags <- factor(plot_df$mags, levels = unique(plot_df$mags))

  
  ggplot(plot_df, aes(x = variable, y = value, fill=mags)) +
    geom_bar(stat='identity', width = 0.6, colour="black") +
    scale_fill_manual(values=c("lightskyblue1", "dodgerblue4", "palegreen3", "palegreen"))+
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    theme_classic()
  ggsave(file_name)

}







# ADDITIONAL FUNCTIONS
# get_scale_factors_dev <- function(sp_test_data, thresholdPerMean=100){
#   # This function takes spike-in data, calculates scaling factors for each 
#   # sample based on spike-in data and returns calculated scaling factors.
#   # befor calculating scaling factors it filters out sets of spike-ins if 
#   # at that concentration not all spike-ins were detected 
#   # it can also remove sets of spike-ins if their deviation [abs(x - mean)]
#   # as % from the mean was larger than thresholdPerMean value
#   # this function also outputs graphs of spike-in behavior for each sample.
#   
#   # Function input:
#   #   sp_test_data: data frame with spike-in data that we use for calculating
#   # scaling factors. Column names are sample names, row names are spike-in names,
#   # MUST contain additional column "rel_conc" that contains real relative
#   # concentration numbers for each spike-in (taken from spike-in metadata file)
#   #   thresholdPerMean: spike-in variation threshold (in %, ex. 20, 50, 100). If
#   # for any given spike-in concentration measured values will deviate from the
#   # mean more than the given threshold - this spike in will be disregarded from
#   # the analysis
#   
#   # Function output is a list consisting from 4 data frames:
#   #   scaleVal_df: dataframe that contains scaling factors for each sample
#   # (intersept and slope). Used to normalize metagenomics data.
#   #   PercentDev_df: dataframe containing percent of deviation from 
#   # the mean for each spike-in concentration for each sample.
#   #   varPass: dataframe of logical values for each spike-in in each sample
#   # indicating if it passed % deviation from the mean threshold (thresholdPerMean)
#   #   detection_pass: dataframe of logical values for each spike-in concentration
#   # in each sample indicating if all spike-ins of each concentration were
#   # detected as expected (FALSE - some or all were not detected, TRUE - all
#   # expected spike-ins were detected)
#   #   stDev_df: dataframe of standard deviations for each concentration set of
#   # spike-ins
#   
#   #
#   #
#   
#   if(any(sp_test_data == 0)){
#     cat("\nWarning! Some spike-ins have zero counts! Check your data!\n\n")
#   }
#   
#   # getting means, standard deviations and percentages of sd from mean for 
#   # each spike in concentration. 
#   mean_sp <- as.data.frame(sp_test_data %>% group_by(rel_conc) %>% summarise_all(mean))
#   sd_sp <- abs(as.data.frame(sp_test_data %>% group_by(rel_conc) %>% summarise_all(sd)))
#   
#   #### THIS should be done better - tranfromations using tydivese?
#   
#   # expanding mean table to fit same size as samples data
#   n.times <- as.vector(table(sp_test_data$rel_conc))
#   mean_expand <- mean_sp %>% slice(rep(1:n(), n.times)) # NEEDS DOUBLE CHECKING!
#   
#   # reordering both tables - mean table and data table
#   sp_test_data <- sp_test_data[order(sp_test_data$rel_conc),]
#   mean_expand <- mean_expand[order(mean_expand$rel_conc),]
#   rownames(mean_expand) <- rownames(sp_test_data)
#   
#   # calculating how far away each measurement deviates from the mean of 
#   # each spike-in concentration set (in %)
#   perDel <- (abs(mean_expand[,-1] - sp_test_data[,-7]) * 100)/mean_expand[,-1]
#   
#   #####
#   
#   # for concentrations where none spike-ins per sample were detected - set 
#   # deviation % to be NA
#   perDel[mean_expand[,-1] == 0] <- NA
#   
#   # getting the table (same size as data table) of logical values where TRUE
#   # means that for this measurement the deviation was below pre-set threshold
#   logical_perDel <- perDel <= thresholdPerMean
#   logical_perDel[mean_expand[,-1] == 0] <- FALSE
#   
#   # preparing perDel for output (adding concentrations column)
#   perDel$rel_conc <- mean_expand$rel_conc
#   
#   # setting logical values FALSE where number of detected 
#   # spikins less then expected
#   sp_counts <- as.data.frame(sp_test_data %>% 
#                                group_by(rel_conc) %>% 
#                                summarise_all(function(x){sum(x!=0)}))
#   logic_detected <- sp_counts[,-1] == as.vector(table(sp_test_data$rel_conc))
#   rownames(logic_detected) <- sp_counts[,1]
#   
#   # printing out warnings about underdetected samples and spike-in concentrations
#   for(i in 1:ncol(logic_detected)){
#     if(any(logic_detected[,i] == FALSE)){
#       cat("\nIn sample ", colnames(logic_detected)[i], 
#           "following concentrations were underdetected:",
#           sp_counts$rel_conc[logic_detected[,i] == FALSE], "\n")
#     }
#   }
#   
#   
#   # relocating "real concentration" column from spike-in dataframe to a vector
#   sample_df <- sp_test_data[ , -which(names(sp_test_data) %in% c("rel_conc"))]
#   rel_conc <- sp_test_data$rel_conc
#   
#   
#   #############
#   # creating new dataframe for scaling factors
#   scale_val <- data.frame(matrix(nrow=2, ncol=ncol(sample_df)))
#   colnames(scale_val) <- colnames(sample_df)
#   rownames(scale_val) <- c("intersept", "slope")
#   
#   #looping through each sample separately
#   for(i in 1:ncol(sample_df))  {
#     #i <- 3
#     sample <- sample_df[,i]
#     
#     # getting the spike-in concentration that pass thresholds
#     detected <- rel_conc %in% rownames(logic_detected)[logic_detected[,i] == TRUE]
#     sp_varnorm <- rownames(sample_df) %in% rownames(logical_perDel)[logical_perDel[,i] == TRUE]
#     select_logical <- detected & sp_varnorm
#     # calculating scaling factors with stripped spike-ins that are located on 
#     # the fare ends of the concentration range and did not pass threshold
#     
#     x <- rel_conc[select_logical]
#     y <- sample[select_logical]
#     lm_res <- lm(y~x)
#     scale_val[, i] <- as.vector(lm_res$coefficients)
#     r_adj <-summary(lm_res)$adj.r.squared # also getting R-sq. adj for plot
#     
#     # colorcoding spike-in concentrations
#     dot_color <- vector()
#     dot_color[!select_logical] <- "red" # did not pass threshold
#     dot_color[select_logical] <- "green" # pass threshold
#     
#     # plot the data into the pdf file in working directory
#     pdf(paste(colnames(sample_df)[i], "_PerDevMean", thresholdPerMean,".pdf", sep = ""))
#     plot(rel_conc, sample, xlab = "real relative concentration", ylab="measured value", 
#          col=dot_color)
#     # cat("Processing i = ", i, "\n")
#     # cat("intercept is ", scale_val[1,i], " and slope is ", scale_val[2,i], "\n\n")
#     abline(scale_val[1,i], scale_val[2,i])
#     legend("topleft", legend=c("variation OUTSIDE the threshold", 
#                                "variation WITHIN threshold",
#                                paste("Adjusted R-squared:  ", round(r_adj, 4), sep = "")),
#            col=c("red", "green", "white"), lty=1, cex=0.8)
#     dev.off()
#   }
#   
#   return(list(scaleVal_df=scale_val, PercentDev_df=perDel, varPass=logical_perDel,
#               detection_pass=logic_detected, stDev_df = sd_sp))
# }
