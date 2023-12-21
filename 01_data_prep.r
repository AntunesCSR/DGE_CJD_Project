# Script that compiles the medical data from the .gz files into a single data table

######################### Step 1: Set Working Directory and Load Libraries

# Load required libraries
library(data.table)  # for fread function
library(dplyr)       # for data manipulation
library(DESeq2)      # for DESeq2 analysis

# Set working directory to the location where the .gz files are stored
setwd("../Omics_project_codes")



######################### Step 2: Compile Data

# Get a list of all .gz files in the directory
directory_path <- "../Omics_project_codes/GSE160208_RAW"
file_list <- list.files(directory_path, pattern = ".*\\.gz$", full.names = TRUE)

# Initialize an empty data table to store the compiled data
compiled_data <- data.table()

# Loop through each file and compile the data
for (file_path in file_list) {
  # Read the expression data from the file
  expression_data <- fread(cmd = shQuote(file_path), header = TRUE, sep = "\t", quote = "")
  

  # Extract the gene counts, names, and class
  counts <- expression_data$Count
  gene_names <- expression_data$Name
  gene_class <- expression_data$CodeClass

  # Extract information about diagnosis, patient number, and sample type from the file name
  parts <- unlist(strsplit(basename(file_path), "_|\\-|\\.", perl = TRUE))
  diagnosis <- parts[5]
  patient_number <- as.numeric(parts[6])
  sample_type <- parts[7]

  # Create a data.table with gene counts and sample information
  sample_data <- data.table(
    Diagnosis = diagnosis,
    PatientNumber = patient_number,
    SampleType = sample_type,
    CountColumn = counts,
    Gene = gene_names,
    Class = gene_class
  )

  # Append the data to the compiled_data data table
  compiled_data <- rbindlist(list(compiled_data, sample_data), use.names = TRUE)
}

# View the compiled data table
View(compiled_data)

# check dimensions of compiled data
dim(compiled_data)

# check dimensions of sample_data
dim(sample_data)

# Export as txt file as sample_data.txt 
write.table(compiled_data, file = "sample_data.txt", sep = "\t", quote = FALSE, row.names = TRUE)

