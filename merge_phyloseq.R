merge_data <- merge_phyloseq(EMP500_104940_phyloseq, EMP500_104945_phyloseq)
merge_data <- merge_phyloseq(merge_data, EMP500_104949_phyloseq)
merge_data <- merge_phyloseq(merge_data, EMP500_104960_phyloseq)
merge_data <- merge_phyloseq(merge_data, EMP500_104961_phyloseq)
merge_data <- merge_phyloseq(merge_data, EMP500_88564_phyloseq)
merge_data <- merge_phyloseq(merge_data, EMP500_88574_phyloseq)

merge_EMP500 <- merge_data 

merge_data <- merge_phyloseq(merge_data, ps_EMP)

merge_data <- merge_phyloseq(merge_data, Tara_phyloseq)

merge_data <- merge_phyloseq(merge_data, Qiita_10178_phyloseq)

merge_data <- merge_phyloseq(merge_data, Qiita_2318_phyloseq)

merge_data <- merge_phyloseq(merge_data, Qiita_1552_phyloseq)


merge_data <- merge_phyloseq(merge_data, ps_TaraFL)

merge_data <- merge_phyloseq(merge_data, ps_TaraLPA)

merge_data <- merge_phyloseq(merge_data, ps_TaraSPA)

merge_atlantic <- merge_phyloseq(ps_TaraFL, ps_TaraLPA)
merge_atlantic <- merge_phyloseq(merge_atlantic, ps_TaraSPA)

metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/MetaData/Filtered_Metadata.csv", row.names=1)
sample_data_object <- sample_data(metadata)
combined_data <- merge_phyloseq(merge_data, sample_data_object)

##EMP

metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/MetaData/EMP.csv", row.names=1)
sample_data_object <- sample_data(metadata)
ps_EMP <- merge_phyloseq(ps_EMP, sample_data_object)
#ps_EMP_TSS <- tss_norm(ps_EMP)

##EMP500
metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/MetaData/EMP500.csv", row.names=1)
sample_data_object <- sample_data(metadata)
merge_EMP500 <- merge_phyloseq(merge_EMP500, sample_data_object)
#merge_EMP500_TSS <- tss_norm(merge_EMP500)

##Qiita 2318
metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/MetaData/Qiita_2318.csv", row.names=1)
sample_data_object <- sample_data(metadata)
Qiita_2318_phyloseq <- merge_phyloseq(Qiita_2318_phyloseq, sample_data_object)
#Qiita2318_TSS <- tss_norm(Qiita_2318_phyloseq)

##Qiita 1552
metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/MetaData/Qiita_1552.csv", row.names=1)
sample_data_object <- sample_data(metadata)
Qiita_1552_phyloseq <- merge_phyloseq(Qiita_1552_phyloseq, sample_data_object)
#Qiita1552_TSS <- tss_norm(Qiita_1552_phyloseq)

##Qiita 10178
metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/MetaData/Qiita_10178.csv", row.names=1)
sample_data_object <- sample_data(metadata)
Qiita_10178_phyloseq <- merge_phyloseq(Qiita_10178_phyloseq, sample_data_object)
#Qiita10178_TSS <- tss_norm(Qiita_10178_phyloseq)

## Atlantic Oceans Project
metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/MetaData/Atlantic_Ocean.csv", row.names=1)
sample_data_object <- sample_data(metadata)
merge_atlantic <- merge_phyloseq(merge_atlantic, sample_data_object)
#mergeAtlantic_TSS <- tss_norm(merge_atlantic)

##Tara Oceans Project
metadata <- read.csv("/Users/pranathir/Desktop/PranathiR/IITM/ComputationalSystemsBiologyLab/Data/Tara5/Metadata.csv", row.names=1)
sample_data_object <- sample_data(metadata)
Tara_phyloseq <- merge_phyloseq(Tara_phyloseq, sample_data_object)
#TaraTSS <- tss_norm(Tara_phyloseq)


variables <- ls()
rm(list = variables[!variables %in% c("combined_data", "ps_EMP", "Qiita_10178_phyloseq", "Qiita_1552_phyloseq", "Qiita_2318_phyloseq", "merge_EMP500", "Tara_phyloseq", "merge_atlantic", "merged_complete")])

