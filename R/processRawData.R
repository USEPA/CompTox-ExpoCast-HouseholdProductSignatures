

#' processRawData
#'
#' @param data_files Character vector of the raw data files complete with the path
#' @param myPath Path to the directory in which you want to save files, plots, etc.
#'
#' @return
#' @export
#'
#' @examples
processRawData <- function(data_files, myPath){
  
  # Each NTA result is in an Excel workbook.  Each workbook has a sheet called
  # List of Samples that describes all samples that make up the other sheets in
  # the workbook.  Work with this sheet to populate a list of sample results.
  
  save_names <- substring(data_files, 58)
  save_names <- gsub(".xlsx", "", save_names)
  
  allData <- list()
  
  for (i in 1:length(data_files)){
    
    metadata <- read_xlsx(data_files[i], sheet = "List of Samples")
    
    # The System ID column gives the name of the sample, Type identifies blanks and 
    # duplicates, and the 2 dilution factor columns tell the dilution value.  All
    # of these make up the sheet names.
    if (save_names[i] == "Shampoo"){
      # Has only 1 dilution factor but there are 2 samples for the blanks
      temp <- metadata[, colnames(metadata) %in% c("System ID", "Type", "Dilution Factor")]
      colnames(temp)[3] <- "value"
      temp$Type[is.na(temp$Type)] <- ""
      temp$Type[temp$Type == "Duplicate"] <- "DUP"
      temp$Type[temp$Type == "Solvent Blank"] <- "SB"
      sheets <- paste(temp$`System ID`, temp$Type, temp$value, sep = "_")
      sheets <- gsub("__", "_", sheets)
      sheets <- c(paste(sheets[1], "_R1", sep = ""), paste(sheets[1], "_R2", sep = ""), paste(sheets[2], "_R1", sep = ""),
                  paste(sheets[2], "_R2", sep = ""), sheets[-c(1,2)])
    } else {
      temp <- metadata[, colnames(metadata) %in% c("System ID", "Type", "Dilution Factor (A)", "Dilution Factor (B)")]
      temp <- melt(temp, id.vars = c("System ID", "Type"))
      temp <- temp[!temp$value == "--",]
      temp$Type[is.na(temp$Type)] <- ""
      temp$Type[temp$Type == "Duplicate"] <- "DUP"
      temp$Type[temp$Type == "Solvent Blank"] <- "SB"
      sheets <- paste(temp$`System ID`, temp$Type, temp$value, sep = "_")
      sheets <- gsub("__", "_", sheets)
    }
    
    # Now read in each sheet to create a list of samples
    data <- list()
    data[[1]] <- metadata
    for (j in 1:length(sheets)){
      data[[j+1]] <- read_xlsx(data_files[i], sheet = sheets[j]) 
    }
    names(data) <- c("List_of_Samples", sheets)
    
    # Save data in this minimally processed form
    #save(data, file = paste(myPath, "data", paste(save_names[i], ".RData", sep = ""), sep = "/"))
    
    allData[[i]] <- data
  }
  names(allData) <- save_names
  
  
  # Make a master list of the List of Samples
  allSamps <- c()
  for (i in 1:length(allData)){
    if (length(allSamps) == 0){
      allSamps <- cbind("Product" = save_names[i], allData[[i]][[1]])
    } else {
      if (save_names[i] == "Shampoo"){
        # Only has 1 dilution factor, so I need to alter this one a bit to be able to 
        # concatenate it with the rest of the product tables
        temp <- allData[[i]][[1]]
        temp2 <- cbind(temp[,1:7], "Dilution Factor (B)" = NA, temp[,8], "Reported Peaks (B)" = NA, temp[,8])
        colnames(temp2)[c(7,9,11)] <- c("Dilution Factor (A)", "Reported Peaks (A)", "Total Reported Peaks (A+B)")
        allSamps <- rbind(allSamps, cbind("Product" = save_names[i], temp2))
      } else {
        allSamps <- rbind(allSamps, cbind("Product" = save_names[i], allData[[i]][[1]]))
      }
    }
  }
  # Now add the full sample names as a column (e.g. B29896_SB_20)
  #allSamps <- cbind("SampName" = NA, allSamps)
  #for (i in 1:save_names){
    
  #}
  # Save this result for future reference
  #save(allSamps, file = paste(myPath, "data", "masterSampleList.RData", sep = "/"))
  
  
  # Further process the data to build chemical x sample matrices
  # Want one that is a binary occurrence matrix for chemicals and the other
  # to give the sample concentrations.
  chemList <- c()
  for (i in 1:length(allData)){
    for (j in 2:length(allData[[i]])){  # 2 b/c 1 is always the list of samples table
      allData[[i]][[j]] <- allData[[i]][[j]][!allData[[i]][[j]]$Group %in% c("is", "s", "ns", "xt;ui", "xu"),]
      if (length(chemList) == 0){
        chemList <- allData[[i]][[j]][,c("Name", "CAS", "Formula", "Group")]
      } else {
        chemList <- rbind(chemList, allData[[i]][[j]][,c("Name", "CAS", "Formula", "Group")])
      }
    }
  }
  
  # Remove duplicated chemicals
  chemList <- chemList[!duplicated(chemList$Name),]
  chemList <- as.data.frame(chemList)
  # Save this set for querying the dashboard
  #write.table(chemList, file = paste(myPath, "data", "uniqueChems_xc_and_xt_all_samples.txt", sep = "/"), quote = FALSE,
   #           sep = "\t", na = "NA", row.names = FALSE)
  
  
  # Now build 3 data frames in the form of chemicals x samples
  # One will be binary (was chemical measured in sample), one will contain the sample concentration,
  # and one will contain the group code for each chemical
  sampNames <- c()
  for (i in 1:length(save_names)){
    sampNames <- c(sampNames, names(allData[[i]])[-1])
  }
  chemXsamp_binary <- as.data.frame(matrix(data = NA, nrow = length(chemList$Name), ncol = length(sampNames),
                                          dimnames = list(chemList$Name, sampNames)))
  chemXsamp_group <- chemXsamp_binary
  chemXsamp_conc <- chemXsamp_binary  
  
  # Populate the data frames
  for (i in 1:length(allData)){
    for (j in 2:length(allData[[i]])){
      ind <- match(chemList$Name, allData[[i]][[j]]$Name)
      indBinary <- ifelse(!is.na(ind), 1, 0)
      chemXsamp_binary[,names(allData[[i]])[j]] <- indBinary
      chemXsamp_group[,names(allData[[i]])[j]] <- allData[[i]][[j]]$Group[ind]
      chemXsamp_conc[,names(allData[[i]])[j]] <- allData[[i]][[j]]$`Sample Concentration`[ind]
    }
  }
  # Save them
  #save(chemXsamp_binary, file = paste(myPath, "data", "chemXsamp_binary_full.RData", sep = "/"))
  #save(chemXsamp_group, file = paste(myPath, "data", "chemXsamp_group_full.RData", sep = "/"))
  #save(chemXsamp_conc, file = paste(myPath, "data", "chemXsamp_conc_full.RData", sep = "/"))
  
}
