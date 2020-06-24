setwd("/Users/DVyshenska/R_wd")
source("sip-norm-functions.R")
# Section one =================================================================

# DATA IMPORT AND PARSING
# importing metadata and coverage files
raw_files_path <- "./read_stat_files/"
sp_meta <- read.csv("spikes_meta.csv", check.names = F, stringsAsFactors = F)
df_list <- parse_bbmap_read(raw_files_path, sp_meta)

# preparing spike-ins data to be an input for get_scale_factors
sequins_df <- df_list$spequins # retreaving only spike-ins data
#######
# step only for sequin article data - filtering out only these sequins
# that are same concentrations across ALL samples
mixes_sequins_meta <- read.csv("sequins_mixes_meta.csv", check.names = F,
                               stringsAsFactors = F)
same_sequins <- mixes_sequins_meta$`Metagenome sequin ID`[mixes_sequins_meta$`Mix A (relative abundance)` == mixes_sequins_meta$`Mix B (relative abundance)`]
sequins_df <- sequins_df[rownames(sequins_df) %in% same_sequins,]
######
# addint metadata (rel_conc) to sequins coverage table
sp_meta <- sp_meta[order(sp_meta[,1]),] 
rownames(sp_meta) <- sp_meta[,1]
sequins_df <- merge(sequins_df, sp_meta["Mix A (relative abundance)"], by=0)
colnames(sequins_df)[which(colnames(sequins_df) == "Mix A (relative abundance)")] <- "rel_conc"
rownames(sequins_df) <- sequins_df$Row.names
sequins_df <- sequins_df[,-which(colnames(sequins_df) == "Row.names")]



# Section two =================================================================

# DATA (reads) INITIAL NORMALIZATION TRANSFORMATION
# normalizing to size of the bin (in sequin article case - genome size)
# getting sizes of scaffolds:
l <- read.csv("./read_stat_files/readstat_1.txt", sep = "\t", header = F,
                 stringsAsFactors = F)
l <- l[-nrow(l), c(1,2)]
l <- l[order(l[,1]),]

# getting mapped raw read numbers for all scaffolds
all_reads <- rbind(df_list$mags, df_list$spequins)
all_reads <- all_reads[order(rownames(all_reads)),]
total_reads <- colSums(all_reads)

# getting RPKM values for each scaffold
rpkm <- data.frame(matrix(nrow = nrow(all_reads), ncol = ncol(all_reads)))
rownames(rpkm) <- rownames(all_reads)
colnames(rpkm) <- colnames(all_reads)
if(all(l[,1] == rownames(all_reads))){
  for(i in 1:ncol(rpkm)){
  rpkm[,i] <- ((all_reads[,i] * 10^6)/total_reads[i])/l[,2]
  }
}

# splitting mags and sequins
mags <- rpkm[!(rownames(rpkm) %in% sp_meta$`Metagenome sequin ID`), ]
sequins_df <- rpkm[(rownames(rpkm) %in% sp_meta$`Metagenome sequin ID`), ]
sequins_df <- merge(sequins_df, sp_meta["Mix A (relative abundance)"], by=0)
colnames(sequins_df)[which(colnames(sequins_df) == "Mix A (relative abundance)")] <- "rel_conc"
rownames(sequins_df) <- sequins_df$Row.names
sequins_df <- sequins_df[,-which(colnames(sequins_df) == "Row.names")]

# doing log2(x+1) tranformation for the counts
log2rpkm_mags <- log2(mags + 1)
log2rpkm_seq <- log2(sequins_df + 1)

# Section three =================================================================
# DATA NORMALIZATIONS


# LOG LM SCALING
# getting linear regression scaling factors
#source("sip-norm-functions.R")
loglm_res <- get_scale_factors_sd(log2rpkm_seq,thresholdPerMean = 100)
# normalizing using linear regression scaling factors
data_norm1 <- lm_normalize(sampl_data = log2rpkm_mags, scaleVal_df = loglm_res$scaleVal_df)
write.csv(data_norm1, "normalized_tables/rpkm_loglmScaling.csv")

input_table <- reorder_df(rbind(data_norm1, c(0,0,0,0,0,0)))
plot_stackbar(plot_df = input_table, file_name = "rpkm_log_lm.pdf")

# NON-LOG LM SCALING
# getting linear regression scaling factors
#source("sip-norm-functions.R")
lm_res <- get_scale_factors_sd(sequins_df,thresholdPerMean = 100)
# normalizing using linear regression scaling factors
data_norm2 <- lm_normalize(sampl_data = mags, scaleVal_df = lm_res$scaleVal_df)
write.csv(data_norm2, "normalized_tables/rpkm_lmScaling.csv")

input_table <- reorder_df(rbind(data_norm2, c(0,0,0,0,0,0)))
plot_stackbar(plot_df = input_table, file_name = "rpkm_lm.pdf")

# RATIO NORMALIZED

#---------
# reads(scaffold)/reads(all_equal_sequins) 

# taking the table with # of all reads and separating it into mag table and 
# sequin table
reads_mags <- all_reads[!(rownames(all_reads) %in% sp_meta$`Metagenome sequin ID`), ]
reads_sequins_df <- all_reads[(rownames(all_reads) %in% sp_meta$`Metagenome sequin ID`), ]
# selecting only sequins added at equal concentrations to all samples:
eq_sp_names <- mixes_sequins_meta$`Metagenome sequin ID`[mixes_sequins_meta$`Mix A (relative abundance)` == mixes_sequins_meta$`Mix B (relative abundance)`]
reads_sequins_df <- reads_sequins_df[rownames(reads_sequins_df) %in% eq_sp_names, ]
# getting the sum of spikein reads per sample - normalization factors
sum_reads_s <- as.data.frame(t(colSums(reads_sequins_df)))
# building normalization table for easier normalization with dataframes
sp_eq_sums <- sum_reads_s[rep(seq_len(nrow(sum_reads_s)), nrow(reads_mags)), ]

# normalizing the data
read_ratio_df <- reads_mags/sp_eq_sums

write.csv(read_ratio_df, "normalized_tables/reads_ratio.csv")
# plotting normalized data
input_table <- reorder_df(read_ratio_df)
plot_stackbar(plot_df = input_table, file_name = "reads_ratio.pdf")

#----------
# (reads(scaffold)/length(scaffold)) / reads(all_equal_sequins)

# taking the table with # of all reads and separating it into mag table and 
# sequin table
reads_mags <- all_reads[!(rownames(all_reads) %in% sp_meta$`Metagenome sequin ID`), ]
reads_sequins_df <- all_reads[(rownames(all_reads) %in% sp_meta$`Metagenome sequin ID`), ]
# dividing read numbers of mags by their lengths
len_mags <- l[l$V1 %in% rownames(reads_mags),]
all(len_mags$V1 == rownames(reads_mags))
for(i in 1:ncol(reads_mags)){
  reads_mags[,i] <- (reads_mags[,i]/len_mags$V2)
}

# selecting only sequins added at equal concentrations to all samples:
eq_sp_names <- mixes_sequins_meta$`Metagenome sequin ID`[mixes_sequins_meta$`Mix A (relative abundance)` == mixes_sequins_meta$`Mix B (relative abundance)`]
reads_sequins_df <- reads_sequins_df[rownames(reads_sequins_df) %in% eq_sp_names, ]
# getting the sum of spikein reads per sample - normalization factors
sum_reads_s <- as.data.frame(t(colSums(reads_sequins_df)))
# building normalization table for easier normalization with dataframes
sp_eq_sums <- sum_reads_s[rep(seq_len(nrow(sum_reads_s)), nrow(reads_mags)), ]

# normalizing the data
read_ratio_df <- reads_mags/sp_eq_sums

write.csv(read_ratio_df, "normalized_tables/readsbylength_ratio.csv")
# plotting normalized data
input_table <- reorder_df(read_ratio_df)
plot_stackbar(plot_df = input_table, file_name = "readsbylength_ratio.pdf")


# RUVg NORMALIZATION OF READS
# using untransformed data!
require("RUVSeq")
raw_files_path <- "./read_stat_files/"
parsed <- parse_bbmap_read(raw_files_path, sp_meta)
set <- as.matrix(rbind(parsed$mags,parsed$spequins))
is.numeric(set)
spikes <- rownames(set) %in% same_sequins

setRUVg <- RUVg(set, spikes, k=1)
write.csv(setRUVg$normalizedCounts, "normalized_tables/reads_RUVg.csv")

# plotting normalized data
bac_logic <- rownames(setRUVg$normalizedCounts) %in% rownames(parsed$mags)
ruvg_bac <- as.data.frame(setRUVg$normalizedCounts[bac_logic,])
input_table <- reorder_df(ruvg_bac)
plot_stackbar(plot_df = input_table, file_name = "ruvg_read_bact.pdf")


