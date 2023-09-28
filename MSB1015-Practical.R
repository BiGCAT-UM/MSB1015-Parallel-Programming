install.packages("data.table")
install.packages("microbenchmark")
install.packages("doParallel")
install.packages("httr")
install.packages("jsonlite")
install.packages("randomForest")
install.packages("caret")
install.packages("BiocManager")
BiocManager::install("EBImage")
#apt update && apt-get install -y libfftw3-dev

library(doParallel)
library(microbenchmark)

detectCores()
detectCores(logical = FALSE)

#In order to measure the time needed to execute code, we will use the "microbenchmark" package
#which tests the code execution time through multiple runs and reports the average time.
# Documentation: https://cran.r-project.org/web/packages/microbenchmark/microbenchmark.pdf

#However, if you want to have a similar functionality in your R scripts aside from "microbenchmark", you can use the following code snippet.

time_runs <- c()

n_cores <- c(1,2,3,4,5,6,7,8)

set.seed(123)

for(n in n_cores){
  
  print(paste0("Starting run ",n,": ",0," seconds"))
  
  start_time <- Sys.time()
  
  ### Start of code
  
  
  ### End of code
  
  end_time <- Sys.time()
  
  time_taken <- round(end_time - start_time, digits = 2)
  
  print(paste0("Time taken for run ",n,": ",time_taken," seconds"))
  
  time_runs <- c(time_runs, time_taken)
}

#####################################################################

library(doParallel)
library(data.table)

setwd('/home/rstudio/workspace/parallel-practical')

setDTthreads(1)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)

setDTthreads(2)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)


setDTthreads(4)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)


setDTthreads(8)

tm <- microbenchmark(fread("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)


tm <- microbenchmark(read.csv("data/augmented_gene_expression.csv"),
                     times = 10L
)
print(tm)


#https://github.com/Rdatatable/data.table/blob/master/src/openmp-utils.c

#################################

library(doParallel)
library(httr)
library(jsonlite)

process_df_row <- function(df){
  
  for(i in 1:nrow(df)) {
    row <- df[i,]
    
    uniprot_id <- row[["uniprot_id"]]
    
    uniprot_url <- paste0("https://rest.uniprot.org/uniprotkb/",uniprot_id,".json")
    
    response <- GET(url = uniprot_url)
    
    response_text <- content(response, "text", encoding = "UTF-8")
    # Parsing data in JSON
    response_json <- fromJSON(response_text)
    
    protein_name <- response_json$proteinDescription$recommendedName$fullName$value
    
    if(is.null(protein_name)){
      protein_name <- ""
    }
    
    protein_length <- response_json$sequence$length
    
    if(is.null(protein_length)){
      protein_length <- ""
    }

    protein_organism <- response_json$organism$scientificName
    
    if(is.null(protein_organism)){
      protein_organism <- ""
    } 
    
    protein_sequence <- response_json$sequence$value
    
    if(is.null(protein_sequence)){
      protein_sequence <- ""
    } 
    
    df[i, 'protein_name'] <- protein_name
    df[i, 'protein_length'] <- protein_length
    df[i, 'protein_organism'] <- protein_organism
    df[i, 'protein_sequence'] <- protein_sequence
    
  }
  
  return(df)
}

uniprot_df <- read.csv('data/uniprot_ids_df.csv')

uniprot_df_enriched <- process_df_row(uniprot_df)


microbenchmark(process_df_row(uniprot_df), times=10L)

microbenchmark( "exp" = {
  
    NUM_OF_CORES <- 4
    
    cl <- makeCluster(NUM_OF_CORES)
    registerDoParallel(cl)
    
    uniprot_df <- read.csv('data/uniprot_ids_df.csv')
    
    df_splits <- split(uniprot_df, (as.numeric(rownames(uniprot_df))-1) %/% 25)
    
    uniprot_df_enriched <- foreach (i = 1:length(df_splits), .combine = rbind, .packages = c("httr", "jsonlite")) %dopar% {
    
        tmp_df = df_splits[[i]]
        
        tmp_df = process_df_row(tmp_df)
    }
    
    stopCluster(cl)
  
  }, times = 10L
)

#################################

library(doParallel)
library(EBImage)

list_of_files <- list.files(path="data/original_images", pattern=".jpg", all.files=FALSE, full.names=FALSE)

length(list_of_files)

process_image <- function(image_file){
  
  original_folder <- "data/original_images"
  processed_folder <- "data/processed_images"
  
  img <- readImage(paste0(original_folder,'/',image_file))
  
  #display(img)
  
  gray_img = channel(img, "gray")
  thresh_img = gray_img < 0.7
  binary_img = bwlabel(thresh_img)
  
  #display(binary_img)

  writeImage(binary_img, paste0(processed_folder,'/',image_file))
}

microbenchmark( "exp" = {
  
    for(image_file in list_of_files){
      process_image(image_file)
    }
    
  }, times = 3L
)


microbenchmark( "exp" = {
  
    NUM_OF_CORES <- 4
    
    cl <- makeCluster(NUM_OF_CORES)
    registerDoParallel(cl)
    
    foreach (i = 1:length(list_of_files), .combine = rbind, .packages = c("EBImage")) %dopar% {
      
      image_file = list_of_files[i]
      
      process_image(image_file)
    }
    
    stopCluster(cl)
    
  }, times = 3L
)

#################################

library(doParallel)
library(caret)

dataset <- matrix(rnorm(13000),nrow=500)

time_runs <- c()

n_cores <- c(2,3,4,5,6,7,8)

set.seed(825)

for(n in n_cores){

  print(paste0("Starting run ",n,": ",0," seconds"))
  
  start_time <- Sys.time()
  
  cl <- makeCluster(n)
  registerDoParallel(cl)
  
  model <- train(V26 ~ ., data = as.data.frame(dataset), method = "rf", verbose = FALSE)  
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  
  time_taken <- round(end_time - start_time, digits = 2)
  
  print(paste0("Time taken for run ",n,": ",time_taken," seconds"))
  
  time_runs <- c(time_runs, time_taken)
}

plot(n_cores, time_runs)

#################################################################


## Exercise 5: Solve the following problem using a parallel programming approach

# 1. Load the CSV file in the "data" folder, named "pdb_ids.csv" using Pandas
# 2. Create a function to process each row of the dataframe. The function should take one argument of type DataFrame and return a dataframe object.
# The function should iterate through the dataframe rows and perform the following steps:
#   * Get the PDB ID from the relevant column and make an HTTP call to download the protein image from PDB (use the following URL template: http://cdn.rcsb.org/images/structures/dl/{}/{}_assembly-1.jpeg).
# * Save the content of the response (binary content) to an image file stored in the folder "data/pdb_images" named with the PDB id and the extension "jpeg"
# * read the image file from the folder using OpenCV and extract the size of the image (i.e. width and height)
# * store the width and the height of the image in the relevant columns in the dataframe
# 3. Save the dataframe to a file
# 4. Use the timing template provided at the beginning of this notebook to time your code and test it using 2, 4 and 8 cores


