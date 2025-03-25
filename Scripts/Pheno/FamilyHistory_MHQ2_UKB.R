library(ukbkings)
library(dplyr)
library(stringr)

project_dir <- "/datasets/ukbiobank/ukb82087"

f <- bio_field(project_dir)

head(f)
glimpse(f)

f %>%
    distinct(basket)
    
f %>%
    select(field, name) %>%
    filter(str_detect(field, "29001")) %>%
    bio_field_add("FH_field_subset.txt")

bio_phen(
    project_dir,
    field = "FH_field_subset.txt",
    out = "FH_field_subset"
)

# read in data 

tep<-readRDS("FH_field_subset.rds")

fh.none<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "0" %in% row)

# family history of depression
fh.depression<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "1" %in% row)
tep$fh.depression[c(fh.depression)]<-1
tep$fh.depression[is.na(tep$fh.depression)&c(fh.none)]<-0

# family history of mania
fh.mania<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "2" %in% row)
tep$fh.mania[c(fh.mania)]<-1

# family history of schizophrenia
fh.scz<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "3" %in% row)
tep$fh.scz[c(fh.scz)]<-1

# family history of psychosis
fh.psychosis<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "4" %in% row)
tep$fh.psychosis[c(fh.psychosis)]<-1

# family history of personality disorder
fh.personality<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "5" %in% row)
tep$fh.personality[c(fh.personality)]<-1

# family history of autism
fh.autism<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "6" %in% row)
tep$fh.autism[c(fh.autism)]<-1

# family history of ADHD
fh.adhd<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "7" %in% row)
tep$fh.adhd[c(fh.adhd)]<-1

# family history of Anxiety
fh.anxiety<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "8" %in% row)
tep$fh.anxiety[c(fh.anxiety)]<-1

# family history of eating disorder
fh.eating_disorder<-apply(tep[,grep("29001", colnames(tep))], 1, function(row) "9" %in% row)
tep$fh.eating_disorder[c(fh.eating_disorder)]<-1

saveRDS(tep,file="FH_field_subset.rds")

# set the 0 and NA for family history

tep$fh.depression.all<-NA
tep$fh.depression.all[tep$fh.depression==1]<-1
tep$fh.depression.all[tep$"29001-0.0">=0&is.na(tep$fh.depression.all)]<-0

tep$fh.mania.all<-NA
tep$fh.mania.all[tep$fh.mania==1]<-1
tep$fh.mania.all[tep$"29001-0.0">=0&is.na(tep$fh.mania.all)]<-0

tep$fh.scz.all<-NA
tep$fh.scz.all[tep$fh.scz==1]<-1
tep$fh.scz.all[tep$"29001-0.0">=0&is.na(tep$fh.scz.all)]<-0

tep$fh.psychosis.all<-NA
tep$fh.psychosis.all[tep$fh.psychosis==1]<-1
tep$fh.psychosis.all[tep$"29001-0.0">=0&is.na(tep$fh.psychosis.all)]<-0

tep$fh.personality.all<-NA
tep$fh.personality.all[tep$fh.personality==1]<-1
tep$fh.personality.all[tep$"29001-0.0">=0&is.na(tep$fh.personality.all)]<-0

tep$fh.autism.all<-NA
tep$fh.autism.all[tep$fh.autism==1]<-1
tep$fh.autism.all[tep$"29001-0.0">=0&is.na(tep$fh.autism.all)]<-0

tep$fh.adhd.all<-NA
tep$fh.adhd.all[tep$fh.adhd==1]<-1
tep$fh.adhd.all[tep$"29001-0.0">=0&is.na(tep$fh.adhd.all)]<-0

tep$fh.anxiety.all<-NA
tep$fh.anxiety.all[tep$fh.anxiety==1]<-1
tep$fh.anxiety.all[tep$"29001-0.0">=0&is.na(tep$fh.anxiety.all)]<-0

tep$fh.eating_disorder.all<-NA
tep$fh.eating_disorder.all[tep$fh.eating_disorder==1]<-1
tep$fh.eating_disorder.all[tep$"29001-0.0">=0&is.na(tep$fh.eating_disorder.all)]<-0

saveRDS(tep,file="FH_field_subset.rds")

# subset family history cleaned datasets

fh<-tep[,c(1,20:28)]
names(fh)<-c("ID","fh.depression","fh.mania","fh.scz","fh.psychosis","fh.personality","fh.autism","fh.adhd","fh.anxiety","fh.eating_disorder")

# generate FH predictors

tep1<-subset(fh,!is.na(tep$fh.depression))
tep1$FID<-tep1$ID
tep1$IID<-tep1$ID

predictors<-c("fh.depression","fh.mania","fh.scz","fh.psychosis","fh.personality","fh.autism","fh.adhd","fh.anxiety","fh.eating_disorder")
predictors_base_path <- "/FH-PRS-MDD/predictors/"

combine_and_write_subset <- function(input_data, predictors_base_path) {
  for (predictor in predictors) {
    # Create a subset data frame with selected predictor
    subset_data <- input_data[, c("FID", "IID", predictor)]

    # Construct the output file path for the current predictor
    output_file <- paste0(predictors_base_path, predictor, ".txt")

    # Write the subset data frame to the output file
    write.table(subset_data, file = output_file, sep = " ", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}

combine_and_write_subset(tep1, predictors_base_path)

q()


