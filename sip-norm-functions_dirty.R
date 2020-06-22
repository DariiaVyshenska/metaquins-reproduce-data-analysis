setwd("/Users/DVyshenska/R_wd")
library(dplyr)

######
# temp helping code
# important training objects:
# sampl_data = simulation of coverage data from samples
# sample_spikedf = simulation of coverage data from spikeins
# spikes_meta = metadata (concentrations) of spikeins

# control tables for checking the functions:
# real_lm = data for the intersept (first row) and slope (second row)
# real_samples = coverage data before adding intersept and slope
spikein_df <- read.csv("spikein_conc.csv", check.names = F, stringsAsFactors = F)
sampl_data <- data.frame(matrix(nrow=100, ncol = 10))
sample_names <- paste("sample", rep(1:10), sep="")
colnames(sampl_data) <- sample_names
rownames(sampl_data) <- paste("mag", rep(1:100), sep = "")



spikein_conc <- spikein_df$`Mix A (relative abundance)`
spikes_meta <- spikein_df[, c(1,4)]
real_lm <- data.frame(matrix(nrow=2, ncol=10))
colnames(real_lm) <- sample_names
real_lm[1, ] <- abs(rnorm(10)*10)
real_lm[2, ] <- abs(rnorm(10)*10)

sample_spikedf <- data.frame(matrix(nrow = 86, ncol = 10))
colnames(sample_spikedf) <- sample_names
rownames(sample_spikedf) <- spikein_df$`Metagenome sequin ID`
for(i in 1:10){
  sample_spikedf[, i] <- real_lm[1, i] + real_lm[2, i]*spikein_conc
}
for(i in 1:100){
  sampl_data[i, ] <- abs(rnorm(10))
}
real_samples <- sampl_data

for(i in 1:10){
  sampl_data[,i] <- real_lm[1, i] + real_lm[2, i]*sampl_data[,i]
}


#######
# FUNCTIONS
lm_normalize <- function(sampl_data, sample_spikedf, spikes_meta){
  # this function scales samples data using linear model (intersept and
  # slope are calculated on spike-in data and then applied to samples data)
  # INPUT
  #
  # sampl_data - a dataframe where row names are MAGs' names, column names
  # are samples' names; the data stored - data that needs to be scaled
  #
  # sample_spikedf - a dataframe where row names are spike-ins' names, column names
  # are samples' names; the data stored - spike-ins data that is used for
  # samles' data scaling
  #
  # spikes_meta - a dataframe with two columns: first contains names of 
  # spikin used in each sample, second contains real concentration of that 
  # spike-in in all samples
  # 
  # OUTPUT is a list of two dataframes:
  # norm_data - contains scaled samples' data (columns are samples, rows are
  # MAGs)
  # scaling_val contains intersept (row 1) and slope (row 2) for each sample
  
  # checking and sorting spike-in data for both files: metadata and reads
  sp_meta <- spikes_meta[order(spikes_meta[,1]),]
  sp_raw <- sample_spikedf[order(rownames(sample_spikedf)),]
  
  if(all(sp_meta[,1] == rownames(sp_raw))){
    print("The spikins names are correct!")
    
    # getting intersept and slope for each sample from spike-in data:
    # creating empty dataframe to store intersept and slope for each sample
    n_val <- data.frame(matrix(nrow=2, ncol=ncol(sampl_data)))
    colnames(n_val) <- colnames(sp_raw)
    # looping through spike-in data for each sample to populate n_val
    for(i in 1:ncol(sp_raw)){
      n_val[,i] <- as.vector(lm(sp_raw[,i]~sp_meta[,2])$coefficients)
    }
  } else{break}
  
  # normalizing input sample data using intersept and slope from spikins:
  # creating empty dataframe to store normalized samples data
  norm_sample_data <- data.frame(matrix(nrow=nrow(sampl_data), ncol = ncol(sampl_data)))
  colnames(norm_sample_data) <- colnames(sampl_data)
  rownames(norm_sample_data) <- rownames(sampl_data)
  # normalizing each sample using calculated values and populating norm_sample_data
  for(i in 1:ncol(sampl_data)){
    norm_sample_data[,i] <- (sampl_data[,i] - n_val[1,i])/n_val[2,i]
  }
  return(list(norm_data = norm_sample_data, scaling_val = n_val))
}
# example:
test_out <- lm_normalize(sampl_data, sample_spikedf, spikes_meta)





#############
# working on limit of quantification
# essential functions for generating test data
get_sample <- function(r_mean, r_sd, n){
  x <- rnorm(n)
  new_x <- r_sd*(x-mean(x))/sd(x)+r_mean
  return(new_x)
}
##
sp_meta_ord <- spikes_meta[order(spikes_meta[,2]),]
sp_meta_ord$sdPer <- NA
new_sp_sam <- data.frame(matrix(nrow=nrow(sp_meta_ord), ncol=ncol(real_samples)))
rownames(new_sp_sam) <- sp_meta_ord[,1]
colnames(new_sp_sam) <- colnames(real_samples)
#e_df <- new_sp_sam

conc <- unique(sp_meta_ord$`Mix A (relative abundance)`)
sp_meta_ord[sp_meta_ord[,2] %in% conc[1:2],3] <- 35
sp_meta_ord[sp_meta_ord[,2] %in% conc[3:5],3] <- 25
sp_meta_ord[sp_meta_ord[,2] %in% conc[6:8],3] <- 15
sp_meta_ord[sp_meta_ord[,2] %in% conc[9:16],3] <- 5

for(i in 1:ncol(new_sp_sam)){
  new_sp_sam[, i] <- real_lm[1, i] + real_lm[2, i]*sp_meta_ord[,2]
}

sd_per <- unique(sp_meta_ord$`Mix A (relative abundance)`)
sim_spdata <- new_sp_sam
sim_spdata[,] <- NA

s <- sd_per[1]

for(s in sd_per){
  sim_sdP <- unique(sp_meta_ord[sp_meta_ord$`Mix A (relative abundance)`==s,3])
  n <- nrow(new_sp_sam[sp_meta_ord$`Mix A (relative abundance)`==s,])
  for(i in 1:ncol(sim_spdata)){
    sim_mean <- mean(new_sp_sam[sp_meta_ord$`Mix A (relative abundance)`==s,i])
    sim_sd <- (sim_sdP*sim_mean)/100
    sim_spdata[sp_meta_ord$`Mix A (relative abundance)`==s,i] <- get_sample(sim_mean, sim_sd, n)
  }
}


# getting simulated mean and standard deviation to make sure our thresholds
# will work on the data
e_dfTemp <- sim_spdata
e_dfTemp$conc <- sp_meta_ord$`Mix A (relative abundance)`
sim_errmeanT <- e_dfTemp %>% group_by(conc) %>% summarise_all(funs(mean))
sim_errsdT <- e_dfTemp %>% group_by(conc) %>% summarise_all(funs(sd))
m <- as.data.frame(sim_errmeanT[,-1])
e <- as.data.frame(sim_errsdT[,-1])
#write.csv(m, "err_m.csv")
#write.csv(sim_mean, "siminit_m.csv")

p <-(e*100)/m
write.csv(p, "pers.csv")

threshCheck <- abs(p) < 10
apply(threshCheck, 1, function(x) all(x))
sum(apply(threshCheck, 1, function(x) all(x)))
apply(threshCheck, 1, function(x) all(!x))
sum(apply(threshCheck, 1, function(x) all(!x)))
write.csv(sim_spdata, "simulated_spikinsdata.csv", row.names = F)


# working on the final function
# I will need to use df sim_spdata as spike-in data to calculate sd % and 
# and move from there

#input test data
# function input
sp_test_data <- read.csv("simulated_spikinsdata.csv", check.names = F, stringsAsFactors = F)
sp_meta <- read.csv("spikes_meta.csv", check.names = F, stringsAsFactors = F)
#samp_test_data <- 
sp_test_data <- sp_test_data[order(sp_test_data[,1]),]
sp_meta <- sp_meta[order(sp_meta[,1]),]
if(all(sp_test_data[,1] == sp_meta[,1])){
  print("ordered correct")
  sp_test_data$rel_conc <- sp_meta[,2]
  rownames(sp_test_data) <- sp_test_data[,1]
  sp_test_data <- sp_test_data[,-1]
  sp_test_data <- sp_test_data[order(sp_test_data$rel_conc),]
}else{
  print("spike-in names in input data tables do not match!")
  break}


# testing functions
sample_test_df <- read.csv("sampl_data.csv", stringsAsFactors = F, check.names = F)
rownames(sample_test_df) <- sample_test_df[,1]
sample_test_df <- sample_test_df[,-1]
source("sip-norm-functions.R")




step1_ls <- get_scale_factors(sp_test_data = sp_test_data, thresholdPerMean = 20)
step2_df <- lm_normalize(sampl_data = sample_test_df, scaleVal_df = step1_ls$scaleVal_df)

step1_lsR <- get_scale_factorsR(sp_test_data = sp_test_data, thresholdPerMean = 20)
step2_dfR <- lm_normalizeR(sampl_data = sample_test_df, scaleVal_df = step1_ls$scaleVal_df)




x <- rnorm(20) + 20
y <- 5 + 2*x

x_lm <- sp_test_data$sample2 #x[1:14]
y_lm <- sp_test_data$rel_conc #y[1:14]
x_test <- x[15:20]
y_test <- y[15:20]
model1 <- lm(y_lm~x_lm)
model2 <- lm(x_lm~y_lm)
pred1 <- model1$coefficients[1] +model1$coefficients[2]*sample_test_df$sample2
pred2 <- (sample_test_df$sample2 - model2$coefficients[1])/model2$coefficients[2]




sp_reads <- sp_test_data$sample2
real_conc <- sp_test_data$rel_conc

new_reads <- sample_test_df$sample2

model1 <- lm(sp_reads ~ real_conc)
model2 <- lm(real_conc ~ sp_reads)

res1 <- (new_reads - model1$coefficients[1]) / model1$coefficients[2]
res2 <- model2$coefficients[1] + model2$coefficients[2]*new_reads


model1$coefficients[1] + model1$coefficients[2]*res1


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
  
  # getting the spike-in concentration that pass thresholds
  conc_to_use <- names(logical)[logical == TRUE]
  
  # calculating scaling factors with stripped spike-ins that are located on 
  # the fare ends of the concentration range and did not pass threshold
  
  x <- sample[rel_conc %in% conc_to_use]
  y <- rel_conc[rel_conc %in% conc_to_use]
  lm_res <- lm(y~x)
  print(colnames(scale_val)[i])
  print(summary(lm_res))
  scale_val[, i] <- as.vector(lm_res$coefficients)
  print(as.vector(lm_res$coefficients))
  print("\n")
  r_adj <-summary(lm_res)$adj.r.squared # also getting R-sq. adj for plot
  
  # colorcoding spike-in concentrations
  dot_color <- vector()
  dot_color[!(rel_conc %in% conc_to_use)] <- "red" # did not pass threshold
  dot_color[rel_conc %in% conc_to_use] <- "green" # pass threshold
  
  # plot the data into the pdf file in working directory
  pdf(paste(colnames(sample_df)[i], ".pdf", sep = ""))
  plot(x, y, xlab = "measured value", ylab="real relative concentration", 
       col=dot_color)
  abline(scale_val[1,i], scale_val[2,i])
  legend("topleft", legend=c("variation OUTSIDE the threshold", 
                             "variation WITHIN threshold",
                             paste("Adjusted R-squared:  ", round(r_adj, 4), sep = "")),
         col=c("red", "green", "white"), lty=1, cex=0.8)
  dev.off()













plot(x,y)
plot(y,x)

############### 
thresholdPerMean <- 20


if(any(sp_test_data == 0)){
  cat("Warning! Some spike-ins have zero counts! Check your data!")
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
# setting logical values for undetected spikins to FALSE
logic_per_sp[mean_sp[,-1] == 0] <- FALSE
per_sp[mean_sp[,-1] == 0] <- NA
# setting logical values and percentages to F and NA where number of detected 
# spikins less then desired
sp_counts <- as.data.frame(sp_test_data %>% 
                             group_by(rel_conc) %>% 
                             summarise_all(function(x){sum(x!=0)}))


logic_detected <- sp_counts[,-1] == as.vector(table(sp_test_data$rel_conc))

for(i in 1:ncol(logic_detected)){
  if(any(logic_detected[,i] == FALSE)){
    cat("In sample ", colnames(logic_detected)[i], 
        "following concentrations were underdetected:",
        sp_counts$rel_conc[logic_detected[,i] == FALSE], "\n\n")
  }
}

logic_per_sp[logic_detected == FALSE] <- FALSE
per_sp[logic_detected == FALSE] <- NA

# relocating "real concentration" column from spike-in dataframe to a vector
sample_df <- sp_test_data[ , -which(names(sp_test_data) %in% c("rel_conc"))]
rel_conc <- sp_test_data$rel_conc

# creating new dataframe for scaling factors
scale_val <- data.frame(matrix(nrow=2, ncol=ncol(sample_df)))
colnames(scale_val) <- colnames(sample_df)
rownames(scale_val) <- c("intersept", "slope")



######

# Poisson regression
colors <- c("Red", "Blue", "Gold", "Black", "Pink", "Green")
poisson.dist <- list()
a < - c(1, 2, 3, 4, 5, 6) # A vector for values of u
for (i in 1:6) {
  poisson.dist[[i]] <- c(dpois(0:20, i)) # Store distribution vector for each corresponding value of u
}
plot(unlist(poisson.dist[1]), type = "o", xlab="y", ylab = "P(y)",
     col = colors[i])
for (i in 1:6) {
  lines(unlist(poisson.dist[i]), type = "o", col = colors[i])
}
# Adds legend to the graph plotted
legend("topright", legend = a, inset = 0.08, cex = 1.0, fill = colors, title = "Values of u")


library(datasets)
data <- warpbreaks
columns <-  names(data) # Extract column names from dataframe
columns # show columns
ls.str(warpbreaks)
poisson.model <-  glm(breaks ~ wool + tension, data, family = quasipoisson(link = "log"))
summary(poisson.model)


#######
sp_test_data <- read.csv("simulated_spikinsdata.csv", check.names = F, stringsAsFactors = F)
sp_meta <- read.csv("spikes_meta.csv", check.names = F, stringsAsFactors = F)
#samp_test_data <- 
sp_test_data <- sp_test_data[order(sp_test_data[,1]),]
sp_meta <- sp_meta[order(sp_meta[,1]),]
if(all(sp_test_data[,1] == sp_meta[,1])){
  print("ordered correct")
  sp_test_data$rel_conc <- sp_meta[,2]
  rownames(sp_test_data) <- sp_test_data[,1]
  sp_test_data <- sp_test_data[,-1]
  sp_test_data <- sp_test_data[order(sp_test_data$rel_conc),]
}else{
  print("spike-in names in input data tables do not match!")
  break}


# testing functions
sample_test_df <- read.csv("sampl_data.csv", stringsAsFactors = F, check.names = F)
rownames(sample_test_df) <- sample_test_df[,1]
sample_test_df <- sample_test_df[,-1]
sp_meta <- sp_meta[order(sp_meta$`Mix A (relative abundance)`),]



r_t <- sample_test_df$sample1
r <- sp_test_data$sample1
c <- sp_test_data$rel_conc

lm(r~c)$coefficients
lm(c~r)$coefficients

8.974323+0.242596*r
plot(r,c)


x <- rnorm(10) + 10
y <- 5+2*x
coef1 <- lm(y~x)$coefficients
coef1[1] + coef1[2]*x

coef2 <- lm(x~y)$coefficients
(x-coef2[1])/coef2[2]

#########
### Function to convert total input ng & sequin rel amounts into atommoles of 
# spikeins

# intputs: spike-in metadata & total input of DNA spikein
# expected inputs (IN ng!!!) for JGI for spikins are "0.1ng of sequins to each well"
# sp_meta - second column must be relative units, third coumn must be spike-in length (bp)
sp_meta <- read.csv("spikes_meta.csv", check.names = F, stringsAsFactors = F)
sp_meta <- sp_meta[order(sp_meta[,2]),]
input_mass <- 0.666

test <- get_molecules(input_mass = input_mass, sp_meta = sp_meta)


##########
getwd()
setwd("/Users/DVyshenska/R_wd")

df <- read.csv("./stat_files/constats_1.txt", sep = "\t", check.names = F)
sp_meta <- read.csv("spikes_meta.csv")
sp_data <- df[df[,1] %in% sp_meta[,1],c(1,2)]
sp_meta$log_conc <- log2(sp_meta$Mix.A..relative.abundance.)

w_data <- merge(sp_meta, sp_data, by.x = "Metagenome.sequin.ID", by.y = "#ID")
w_data[w_data$Avg_fold == 0,5] <- NA
w_data$log_cov <- log2(w_data$Avg_fold)
plot(w_data$log_conc, w_data$log_cov)
summary(lm(log_cov~log_conc, data=w_data))
