

fileLines <- readLines("eQTL_QTLtools-18759968.out")


submittedJobs <- numeric(0)
completedJobs <- numeric(0)
pattern <- "\\b\\d+\\b"
for (line in fileLines){
  # Extracted submitted jobs
  if (grepl("^Submitted job", line)) {
    # Pull out first snakemake job id
    jobNumber <- as.numeric(unlist(regmatches(line, gregexpr(pattern, line)))[1])
    submittedJobs <- c(submittedJobs, jobNumber)
  }
  # Extract completed jobs
  if (grepl("^Finished job", line)){
    jobNumber <- as.numeric(unlist(regmatches(line, gregexpr(pattern, line)))[1])
    completedJobs <- c(completedJobs, jobNumber)
  }
  
}

incompleteJobs <- submittedJobs[which(!submittedJobs %in% completedJobs)]
