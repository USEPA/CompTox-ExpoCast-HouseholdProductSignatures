

#' processSSAfiles
#'
#' @param data_files Character vector of the raw data files complete with the path
#' @param myPath Path to the directory in which you want to save files, plots, etc.
#' @param prodGroups Full, easily discernible product category names to be used in
#'                   downstream analyses/plotting/etc.
#'
#' @return
#' @export
#'
#' @examples
processSSAfiles <- function(data_files, myPath, prodGroups){
  
  # All data files are in the form path/to/directory/Sequence# prod cat.xlsx.
  # Use the curated prodGroups object to save files and populate data objects.
  
  # Each SSA result is in an Excel workbook.  Each workbook has a sheet called
  # List of Samples that describes all samples that make up the other sheets in
  # the workbook.  Work with this sheet to populate a list of sample results.
  
  allData <- list()
  for (i in 1:length(data_files)){
    
    metadata <- read_xlsx(data_files[i], sheet = "List of Samples")
    
    # The System ID column gives the name of the sample, Type identifies blanks and 
    # duplicates, and the 2 dilution factor columns tell the dilution value.  All
    # of these make up the sheet names.
    if (prodGroups[i] == "Shampoo"){
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
      temp <- reshape2::melt(temp, id.vars = c("System ID", "Type"))
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
    #save(data, file = paste(myPath, "data", paste(prodGroups[i], ".RData", sep = ""), sep = "/"))
    
    allData[[i]] <- data
  }
  names(allData) <- prodGroups
  
  # allData is a list (1 for each product category) of lists (1 object for each sheet in the 
  # raw Excel file, which include the detected peaks in each sample and a master table for 
  # each category, which is the first object).
  
  # Make a master list of the List of Samples
  allSamps <- c()
  for (i in 1:length(allData)){
    if (length(allSamps) == 0){
      allSamps <- cbind("Product" = prodGroups[i], allData[[i]][[1]])
    } else {
      if (prodGroups[i] == "Shampoo"){
        # Only has 1 dilution factor, so I need to alter this one a bit to be able to 
        # concatenate it with the rest of the product tables
        temp <- allData[[i]][[1]]
        temp2 <- cbind(temp[,1:7], "Dilution Factor (B)" = NA, temp[,8], "Reported Peaks (B)" = NA, temp[,8])
        colnames(temp2)[c(7,9,11)] <- c("Dilution Factor (A)", "Reported Peaks (A)", "Total Reported Peaks (A+B)")
        allSamps <- rbind(allSamps, cbind("Product" = prodGroups[i], temp2))
      } else {
        allSamps <- rbind(allSamps, cbind("Product" = prodGroups[i], allData[[i]][[1]]))
      }
    }
  }
  # Save this result for future reference
  #save(allSamps, file = paste(myPath, "data", "masterSampleList.RData", sep = "/"))
  
  
  # Save another table for the manuscript supplement. Concatenate all raw output files for
  # each sample (contained in allData) to have a complete list of all identified chemicals
  # across all samples in the study.
  fullSampChemInfo <- c()
  for (i in 1:length(allData)){
    for (j in 2:length(allData[[i]])){
      temp <- allData[[i]][[j]][,2:25]
      #dil <- strsplit(names(allData[[i]]), "_(?!.*_)", perl = TRUE)[[j]][2]
      temp <- cbind(rep(prodGroups[i], dim(temp)[1]), temp)
      fullSampChemInfo <- rbind(fullSampChemInfo, temp)
    }
  }
  # Change some column names
  colnames(fullSampChemInfo)[1:3] <- c("Category", "MySysID", "Customer ID")
  save(fullSampChemInfo, file = paste(myPath, "data", "fullSetOfSamplesWithAllChems.RData", sep = "/"))
  
  
  ## Further process the data to build chemical x sample matrices
  ## Want one that is a binary occurrence matrix for chemicals, one for the group codes,
  ## and one to give the sample concentrations.
  ## Do for in 2 sets.  1. xc and xt chemicals and 2. ns, xt;ui, and xu chemicals
  #-----------------------------------------------------------------------------
  # Set 1 (chemicals we will carry out our analysis pipeline on)
  ind <- fullSampChemInfo$Group %in% c("xc", "xt")
  chemList <- fullSampChemInfo[ind, c("Name", "CAS", "Formula", "Group")]
  # Remove duplicated chemicals
  chemList <- chemList[!duplicated(chemList$Name),]
  # Save this set for querying the dashboard
  #write.table(chemList, file = paste(myPath, "data", "uniqueChems_xc_and_xt_all_samples.txt", sep = "/"), quote = FALSE,
   #           sep = "\t", na = "NA", row.names = FALSE)
  
  # Now build 3 data frames in the form of chemicals x samples
  sampNames <- c()
  sampCats <- c()
  for (i in 1:length(prodGroups)){
    sampNames <- c(sampNames, names(allData[[i]])[-1])
    sampCats <- c(sampCats, rep(prodGroups[i], length(allData[[i]][-1])))
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
  chemXsamp_matrices_xcxt <- list("binary" = chemXsamp_binary, "group" = chemXsamp_group,
                                  "conc" = chemXsamp_conc, "category" = sampCats)
  #save(chemXsamp_matrices_xcxt, file = paste(myPath, "data", "chemXsamp_matrices_xcxt_full.RData", sep = "/"))

  #-----------------------------------------------------------------------------
  # Set 2 (chemicals we will check for prevalence as prioritization)
  ind <- fullSampChemInfo$Group %in% c("ns", "xt;ui", "xt;ui;PCB", "xu")
  chemList <- fullSampChemInfo[ind, c("Name", "CAS", "Formula", "Group")]
  # Remove duplicated chemicals
  chemList <- chemList[!duplicated(chemList$Name),]
  # Save this set for querying the dashboard
  #write.table(chemList, file = paste(myPath, "data", "uniqueChems_ns-xtui-xtuiPCB-xu_all_samples.txt", sep = "/"), quote = FALSE,
    #         sep = "\t", na = "NA", row.names = FALSE)
  
  # Now build 3 data frames in the form of chemicals x samples
  sampNames <- c()
  sampCats <- c()
  for (i in 1:length(prodGroups)){
    sampNames <- c(sampNames, names(allData[[i]])[-1])
    sampCats <- c(sampCats, rep(prodGroups[i], length(allData[[i]][-1])))
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
  chemXsamp_matrices_other <- list("binary" = chemXsamp_binary, "group" = chemXsamp_group,
                                  "conc" = chemXsamp_conc, "category" = sampCats)
  #save(chemXsamp_matrices_other, file = paste(myPath, "data", "chemXsamp_matrices_ns-xtui-xtuiPCB-xu_full.RData", sep = "/"))
  
}




#' cleanForUPCs
#'
#' @param upc_file 
#' @param myPath 
#' @param prodGroups
#'
#' @return
#' @export
#'
#' @examples
cleanForUPCs <- function(upc_file, myPath, prodGroups){
  
  myData <- list()
  for (i in 1:length(prodGroups)){
    myData[[i]] <- as.data.frame(read_xlsx(paste(myPath, "raw_data", upc_file, sep = "/"), sheet = i, col_names = TRUE, skip = 1))
    
    # An empty column separates the 2 task orders
    indNA <- which(is.na(myData[[i]][1,]))
    temp <- myData[[i]][,(indNA+1):dim(myData[[i]])[2]]
    myData[[i]] <- myData[[i]][,1:(indNA-1)]
    
    # Some sheets have a notes column at the end to provide more details about the samples. Create one if not there.
    if (length(grep("Notes", colnames(temp))) == 0){
      colnames(temp) <- sapply(colnames(temp), function(x) strsplit(x, "[.]")[[1]][1])
      colnames(myData[[i]]) <- colnames(temp)
      myData[[i]] <- rbind(myData[[i]], temp)
      myData[[i]]$Notes <- NA
    } else {
      colnames(temp) <- sapply(colnames(temp), function(x) strsplit(x, "[.]")[[1]][1])
      myData[[i]]$Notes <- NA
      colnames(myData[[i]]) <- colnames(temp)
      myData[[i]] <- rbind(myData[[i]], temp)
    }
    
    # Add a column for product category
    myData[[i]] <- cbind("Category" = rep(prodGroups[i], dim(myData[[i]])[1]), myData[[i]])
  }

  # Make the final table
  myProdUPCs <- as.data.frame(rbindlist(myData, fill = TRUE))
  # There are some extra rows due to unneeded info at the bottom of the sheets
  myProdUPCs <- myProdUPCs[!is.na(myProdUPCs$`SwRI ID`),]
  
  # Save result
  save(myProdUPCs, file = paste(myPath, "data", "prodPurchaseDetails_withUPCs.RData", sep = "/"))
  
}












