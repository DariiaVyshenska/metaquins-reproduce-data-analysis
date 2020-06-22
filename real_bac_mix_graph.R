setwd("/Users/DVyshenska/R_wd")
source("sip-norm-functions.R")


bac_table <- read.csv("cyano_mixes_metaquins.csv", row.names = 1, check.names = F)
#plot_stackbar(bac_table, "article_stoch.pdf")

# have to multiply first three mixes - otherwise I can't recreate the
# graph from the article
bac_table_mult[,1:3] <- bac_table[,1:3]*2
plot_stackbar(bac_table_mult, "article_stoch.pdf")
