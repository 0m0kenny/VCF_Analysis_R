#Download necessary packages


library("vcfR")
library("ggplot2")
library("tidyr")
library("dplyr") 
library("RColorBrewer")
library("tidyverse")
library("tibble")


# Reading the file and getting appropriate data


x <- read.vcfR("DellyVariation.vcf") 

#check column names to find column no of relevant samples
colnames(x@gt) 


#Extract relevant sample columns in vcf file gt region


# Subset genotype matrix to exclude the unwanted sample columns
x@gt <- x@gt[, -c(2:7)]

# Save the filtered VCF file
write.vcf(x, "Dellyfiltered.vcf")


# Load the new file


#load the vcf file
y <- read.vcfR("Dellyfiltered.vcf") 

#check column names to confirm presence of only relevant samples
colnames(y@gt) 




# A- Summarise the various combinations of variants, their location, type and function

### Summarise the vcf file



y[is.biallelic(y)] # gives summary of the total number of biallelic variants

y[is.polymorphic(y, na.omit = TRUE),] #gives summary of the total no of polymorphic variants- omit data with na in them


Results show that all variants are biallelic and only 10,197 are polymorphic

### Used vcfRtidy() function to extract each element in the info and format field so it is easier to extract specific information outside of the vcfR package


#separate the individual elements in the info and format fields into separate rows
Z <- vcfR2tidy(y, single_frame = TRUE) 

#gives info on the abbreviations used in the format and info fields stored in meta
col_abrv <-Z$meta 

#extract just the elements from the data region
data_out <-Z$dat 

# separate each sample into different columns
z_wide <- data_out %>%
  pivot_wider(
    names_from = Indiv, # Use the 'Indiv' column to create new columns
    values_from = c(
      gt_GT, gt_GL, gt_FT, gt_RC, gt_RCL, gt_RCR, gt_CN, gt_DR, 
      gt_DV, gt_RR, gt_RV, gt_GQ, gt_GT_alleles  # Columns to pivot
    ),
    names_glue = "{.value}_{Indiv}"  #Format the new column names to differentiate (e.g., gt_GT_Sample1)
  )



## Summarise type of variants

### Extract all the variants types and plot the frequency of each type.


#check for all possible variant types in the vcf file
queryMETA(y)


#specifies all possible variant types using the abbreviations specified in the ID column
var_types<- c("DEL", "INS", "DUP", "INV", "TRA") 

#find rows and extract rows with specified variants in the ID column then count total no of rows 
del_var <- nrow(z_wide[grep("DEL", z_wide$ID), ]) 
ins_var <- nrow(z_wide[grep("INS", z_wide$ID), ])
dup_var <- nrow(z_wide[grep("DUP", z_wide$ID), ])
inv_var <- nrow(z_wide[grep("INV", z_wide$ID), ])
tra_var <- nrow(z_wide[grep("TRA", z_wide$ID), ]) 

#store the total no of each variant type in a list
var_count <- c(del_var, ins_var, dup_var, inv_var, tra_var)

#create a data frame of no of variants 
var_count_df<- cbind(var_types,var_count)

#plot the total number of each variant type 
ggplot(var_count_df, aes(x = var_types, y = as.numeric(var_count), fill = var_types)) + 
  geom_bar(stat = "identity") + 
  labs(title = "Frequency of Variant Types", 
       x = "Types of Variants", 
       y = "Frequency") + 
  theme(plot.title =element_text(hjust=0.5))
  
        


#There are no inversions or translocations but there are 11,823 deletions, 173 duplications and 5,316 insertions.

#This is a lot of data so filtering is required.

### Perform QC on the data

#Remove duplicates and filter for variants that pass in the Filter and the FT column for each sample, that have a Genotype Quality \>= 20, Read Coverage \>= 40 and variants where Precise = TRUE.


#Filter data based on specific conditions
qc_passz <- z_wide %>% distinct() %>%
  filter(if_all(c(gt_FT_NA19238, gt_FT_NA19239, gt_FT_NA19240, FILTER), ~ . == "PASS")) %>%
  filter(if_all(c(gt_GQ_NA19238, gt_GQ_NA19239, gt_GQ_NA19240), ~ . >= 20)) %>%
  filter(if_all(c(gt_RC_NA19238, gt_RC_NA19239, gt_RC_NA19240), ~ . >= 40)) %>%
  filter(PRECISE == TRUE)



### Re-extract all the variant types and plot frequency of each type


#specifies all possible variant types using the abbreviations specified in the ID field
var_types_qc<- c("DEL", "INS", "DUP") 

#count total no of rows with specified variants
del_var_qc <- qc_passz[grep("DEL", qc_passz$ID), ]
ins_var_qc <- qc_passz[grep("INS", qc_passz$ID), ]
dup_var_qc <- qc_passz[grep("DUP", qc_passz$ID), ]

#store the total no of each variant type in a list
var_count_qc <- c(nrow(del_var_qc), nrow(ins_var_qc), nrow(dup_var_qc))

#create a data frame of no of variants 
vcount_df_qc <- cbind(var_types_qc,var_count_qc)

#plot the total number of each variant type 
ggplot(vcount_df_qc, aes(x = var_types_qc, y = as.numeric(var_count_qc), fill = var_types_qc)) + 
  geom_bar(stat = "identity") + 
  labs(title = "Frequency of Variant Types After Additional QC", 
       x = "Types of Variants", 
       y = "Frequency") + 
  theme(plot.title =element_text(hjust=0.5), legend.position = "none")
  
        


#Now we see 11,295 deletions (528 less than before filtering), 169 duplications (4 less than before) and 0 insertions.

### Find chromosomes with most variants


#initialise an empty matrix to store counts
chr_counts <- matrix(0, nrow = 23, ncol =3)

#change the column names in the matrix
colnames(chr_counts) <- c("Chromosome", "Deletions", "Duplications")

#loop 23 times reflecting the 22 different chromosomes
for (chrom in 1:22){
      chrom_pos <- paste("chr", chrom, sep = "") #creates a string chr1, chr2 etc according to the CHROM column
      
      #counts the total no of chr1 in the deletion variant df and the duplication variant df and stores it
      chr_counts[chrom,] <- c(chrom_pos, sum(del_var_qc[,1] == chrom_pos), sum(dup_var_qc[,1] == chrom_pos))
      
      #for the last chromosome, search for chrX 
      chr_counts[23, ] <- c("chrX", sum(del_var_qc[,1]  == "chrX"), sum(dup_var_qc[,1] == "chrX"))
}

#turn the counts matrix into a dataframe
chr_counts_df  <- as.data.frame(chr_counts)

#group the counts according to the variant type
chr_counts_df <- pivot_longer(chr_counts_df, 
                             cols = c("Deletions", "Duplications"), 
                             names_to = "Variant Type", 
                             values_to = "Total number")

#forces ggplot to plot the graph according to row number
chr_counts_df$Chromosome <- factor(chr_counts_df$Chromosome, 
                                     levels = c(paste("chr", 1:22, sep = ""), "chrX"))


#plot the total number of each variant type 
ggplot(chr_counts_df, aes(x = Chromosome, y = as.numeric(`Total number`) , fill = `Variant Type`)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Frequency of Variant Types in different chromosomes", 
       x = "Chromosome", 
       y = "Total no of variant") + 
  theme(plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 90, hjust = 0.5))



#Chromosome 1 has the most deletions but chromosome 2 has the most duplications.

# B- Identify De novo mutations

### Determine Parents and Child

#Using the incompatible allele combination, if we get a response we know its not the child. We expect to see responses only if the sample filtered for 0/1 genotype is not the child because as a child this genotype is impossible if the parents are both 1/1 unless a de novo mutation occurred. The genotype 1/0 is the same as 0/1 so we did not include this in the code.



#filter using genotype to find variants where all conditions are true
is_notchild_38 <- filter(qc_passz, 
                             (gt_GT_NA19238 == "0/1") & 
                             (gt_GT_NA19239 == "1/1") & 
                             (gt_GT_NA19240 == "1/1"))

is_notchild_39 <- filter(qc_passz, 
                             (gt_GT_NA19238 == "1/1") &
                             (gt_GT_NA19239 == "0/1") & 
                             (gt_GT_NA19240 == "1/1"))

is_notchild_40 <- filter(qc_passz, 
                             (gt_GT_NA19238 == "1/1") & 
                             (gt_GT_NA19239 == "1/1" ) & 
                             (gt_GT_NA19240 == "0/1"))





#Based on inheritance pattern we see that NA19238 and NA19239 are the parents and NA19240 is the child.

### Identify de novo mutation(s)


# check for de novo mutations where both parents are homozygous for ref or alt allele
de_novo_hom <- filter(qc_passz,
                             ((gt_GT_NA19238 == "1/1") & 
                             (gt_GT_NA19239 == "1/1" ) & 
                             (gt_GT_NA19240 == "0/1"| gt_GT_NA19240 == "0/0")) |
                             ((gt_GT_NA19238 == "0/0") & 
                             (gt_GT_NA19239 == "0/0") & 
                             (gt_GT_NA19240 == "0/1"| gt_GT_NA19240 == "1/1")) ) 

# check for de novo mutations when one parent is homozygous for ref
de_novo_38 <- filter(qc_passz,
                             ((gt_GT_NA19238 == "1/1") & 
                             (gt_GT_NA19239 == "0/1") & 
                             (gt_GT_NA19240 == "0/0")) |
                             ((gt_GT_NA19238 == "1/1") & 
                             (gt_GT_NA19239 == "0/0" ) & 
                             (gt_GT_NA19240 == "0/0" | gt_GT_NA19240 == "1/1")))

de_novo_39 <- filter(qc_passz,
                             ((gt_GT_NA19238 == "0/1") & 
                             (gt_GT_NA19239 == "1/1" ) & 
                             (gt_GT_NA19240 == "0/0")) | 
                             ((gt_GT_NA19238 == "0/0") & 
                             (gt_GT_NA19239 == "1/1" ) & 
                             (gt_GT_NA19240 == "0/0" | gt_GT_NA19240 == "1/1")))
                     


#combine all possible de novo into a single dataframe
all_de_novo <- rbind(de_novo_hom, de_novo_38, de_novo_39)



#There are only 3 de novo mutations occurring in the child.

# C- Annotate the variants

#### Variant Effect Predictor (VEP) web interface was used to annotate the vcf file and the results were downloaded as a txt file and uploaded on Rstudio.


#load the annotated file and convert to a data frame
annotated_df <-  read.csv("VEPannotated.txt", header = TRUE, sep = "\t" ) %>% as.data.frame() 

#check column names
colnames(annotated_df)




#### Filter the annotated variants to extract only the variants that passed the previous QC filtering in raw VCF.



#get total number of variants in qc filtered dataframe
n_vep <- nrow(qc_passz)

#initialise an empty list to store filtered variants
filtered_results <- list()

#for loop to pick each variant that passed qc and extract the corresponding variants in the annotated file.
for (i in 1:n_vep){
  #select variant ID
  var_id <- qc_passz[i,3]
  
  #extract variants in the annotated file that have the same ID
  filtered_var <- annotated_df %>%
  filter(X.Uploaded_variation %in% var_id)
  
  #store each filtered data into the list
  filtered_results[[i]] <- filtered_var
}

#convert list into a dataframe
filtered_annotated_df <- bind_rows(filtered_results)



### Summary of Annotated Variants

##### Percentage counts of all consequences


conseq_all <- c(
  "intron_variant",
  "non_coding_transcript_variant",
  "downstream_gene_variant",
  "upstream_gene_variant",
  "NMD_transcript_variant",
  "non_coding_transcript_exon_variant",
  "intergenic_variant",
  "feature_truncation",
  "regulatory_region_variant")

n2 <- length(conseq_all)

#initialise an empty matrix to store counts
conseq_all_counts <- matrix(0, nrow = n2, ncol =2)

#change the row and column names in the matrix
rownames(conseq_all_counts) <- conseq_all
colnames(conseq_all_counts) <- c("Counts", "%Counts")


 for (con in 1:n2){
   #count total no of times the coding consequence appears in the consequence column
       conseq_all_counts[con,1] <-  sum(stringr::str_detect(filtered_annotated_df$Consequence, conseq_all[con]))
 }

for (con in 1:n2){
   conseq_all_counts[con,2] <- round((conseq_all_counts[con,1]/sum(conseq_all_counts[,1])) * 100, digits = 0)
}



###### Pie chart to summarise the coding sequences



#convert the counts matrix into a dataframe
conseq_all_counts_df  <- as.data.frame(conseq_all_counts)

#detect row names
has_rownames(conseq_all_counts_df)

#give the rows a header
conseq_all_counts_df <- rownames_to_column(conseq_all_counts_df, var = "Consequences")

#create a new column for labels combining the consequence and %counts column and adding% eg downstream_variant: 6%
conseq_all_counts_df$labels <- paste0(conseq_all_counts_df$Consequence, ": ", conseq_all_counts_df$`%Counts`, "%")

#create a list from the values in the label column
ordered_labels <- as.list(conseq_all_counts_df$labels)

# Reorder the labels based on the list to force ggplot to plot according to order
conseq_all_counts_df$labels <- factor(conseq_all_counts_df$labels, levels = ordered_labels)


ggplot(conseq_all_counts_df, aes(x = "", y = `%Counts`, fill = labels)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Coding Sequences (all)") +
  theme_void() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust=0.5)) +
  scale_fill_manual(
    values = setNames(RColorBrewer::brewer.pal(n = n2, name = "Set3"), conseq_all_counts_df$labels)) # Customizing the legend labels to include percentages
   



##### Percentage of variants that are in the coding regions



Coding_conseq<- c(
"coding_sequence_variant",
"stop_lost",
"inframe_deletion",
"frameshift_variant",
"inframe_insertion",
"stop_gained",
"start_lost",
"start_retained_variant")

n1 <- length(Coding_conseq)

#initialise an empty matrix to store counts
coding_conseq_counts <- matrix(0, nrow = n1, ncol =2)

#change the column names in the matrix
rownames(coding_conseq_counts) <- Coding_conseq
colnames(coding_conseq_counts) <- c("Counts", "%Counts")

#for loop to search for consequence in each variant and count total no
for (con in 1:n1){
      coding_conseq_counts[con,1] <-  sum(stringr::str_detect(filtered_annotated_df$Consequence,  Coding_conseq[con]))

}

#for loop to get counts in percetange
for (con in 1:n1){
  coding_conseq_counts[con,2] <- round((coding_conseq_counts[con,1]/sum(coding_conseq_counts[,1])) *100, digits = 0)
}



###### Pie chart to summarise the percentage count of the coding sequences


#convert the counts matrix into a dataframe
coding_conseq_counts_df  <- as.data.frame(coding_conseq_counts)

#detect row names
has_rownames(coding_conseq_counts_df)

#give the rows a header
coding_conseq_counts_df <- rownames_to_column(coding_conseq_counts_df, var = "Consequences")


#create a new column for labels combining the consequence and %counts column and adding% eg downstream_variant: 6%
coding_conseq_counts_df$labels <- paste0(coding_conseq_counts_df$Consequence, ": ", coding_conseq_counts_df$`%Counts`, "%")

#create a list from the values in the label column
ordered_labels <- as.list(coding_conseq_counts_df$labels)

# Reorder the labels based on the list to force ggplot to plot according to order specified in list
coding_conseq_counts_df$labels <- factor(coding_conseq_counts_df$labels, levels = ordered_labels)

#plot the pie chart
ggplot(coding_conseq_counts_df, aes(x = "", y = `%Counts`, fill = labels)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "Coding consequences") +
  theme_void() +
  theme(legend.title = element_blank(),plot.title = element_text(hjust=0.5)) +
  scale_fill_manual(
    values = setNames(RColorBrewer::brewer.pal(n = n1, name = "Set2"), coding_conseq_counts_df$labels)) # Customizing the legend labels to include percentages
   



### Annotations of the De novo mutations

##### Extract the annotations for the de novo mutations.


#filter each de novo mutation by ID
de_novo1<- filter(filtered_annotated_df, 
                  (X.Uploaded_variation == "DEL00062209"))

de_novo2<- filter(filtered_annotated_df, 
                  (X.Uploaded_variation == "DEL00064459"))

de_novo3<- filter(filtered_annotated_df, 
                  (X.Uploaded_variation == "DEL00056339"))




## Using IGV, we looked at a list of proneural development genes to identify variants in these genes. 

# Found 2 genes with variants: DCX and TCF4. 

# The insertion variant affecting DCX was filtered out using our QC conditions. 

# However, the variant affecting TCF4 passed QC. 

#extract the annotations for the variant
tcf<- filter(filtered_annotated_df, 
                  (X.Uploaded_variation == "DEL00073588"))



##We then looked at the data using different functions to find specific information about the variants.

#EXCTRACT ALL THE PHENOTYPES ASSOCIATED WITH THE VARIANT
unique(tcf$PHENOTYPES)
decay <- filter(tcf, (BIOTYPE == "nonsense_mediated_decay"))


