#' A function to generate an FPWM class object.
#'
#' This function creates the FPWM object and performs the fork in the motif.
#' @param mainTF [character] character with the name of the main TF.
#' @param partners [character list or vector] List or character vector with the names of the partner TFs.
#' @param cell [character] with the name of the main TF.
#' @param mainTF_MMID [character] with the name of the main TF.
#' @param partners_MMID [character vector] with the name of the main TF.
#' @param forkPosition [numeric] defines the postion in the matrix where the motif will be forked.
#' @param probabilityMatrix [logical] whether the function should return a frequency matrix or probability matrix (Default FALSE). 
#' @param scaleFrequencyCounts [logical] whether the count matrix should have equal rowSums across all the rows (Default FALSE). 
#' @param flipMatrix [logical] whether to apply reverse complement in case the core motif is after the forkPosition (Default FALSE). 
#' @examples
#' fpwm <- createFPWM(mainTF ="CEBPB", partners = c("ATF4","ATF7","JUND","CEBPD"), cell = "K562", forkPosition = 5)
#' @return returns a FPWM class object that can be used to plot or write in transfact format.
#' @export
createFPWM <- function( mainTF = NULL,
						partners = NULL,
						cell = NULL,
						mainTF_MMID = NULL,
						partners_MMID = NULL,
						forkPosition = NULL,
            probabilityMatrix = FALSE,
            scaleFrequencyCounts = FALSE,
            flipMatrix = FALSE)
						{
  
  tfnames <- FALSE
  tfIDs <- FALSE
  
  # Check if inputs are compatible
  if(!is.null(mainTF) & !is.null(partners) & !is.null(cell) & !is.null(forkPosition) ) { tfnames <- TRUE }
  if(!is.null(mainTF_MMID) & !is.null(partners_MMID) & !is.null(forkPosition) ) { tfIDs <- TRUE }
  
  if(tfnames & tfIDs) { stop("Incompatible input. Provide 'mainTF', 'partners' and 'cell'. !!OR!! 'mainTF_MMID' and 'partners_MMID'. ") }

  if(!(tfnames | tfIDs) ) { stop("Missing input.\nProvide: 'mainTF', 'partners', 'cell' and 'forkPosition' \nOR\n 'mainTF_MMID', 'partners_MMID' and 'forkPosition'") }

  if( probabilityMatrix == TRUE & scaleFrequencyCounts == TRUE ) { stop("scaleFrequencyCounts only works for count matrices, disable probabilityMatrix.") }

  if (tfnames){
        xTF_cell_tissue_name <- suppressMessages(dataBrowser(tf = mainTF, cell_tissue_name = cell)$ID[1]) 
  
        partners <- as.list(unlist(partners)) ; 
        MMpartners <- partners
        for(i in 1:length(partners)){ suppressMessages( MMpartners[[i]] <- dataBrowser(tf = partners[[i]], cell_tissue_name = cell)$ID[1] ) }
        
	    if(is.null(xTF_cell_tissue_name)){ stop("Please check the spelling of your mainTF or cell. \nOr there might no be information for this TF-cell combination.") }
	
	    if( length(xTF_cell_tissue_name) != 1 ){ stop("More than one record for the combination of TF and cell tissue") ; } 
  
  		if( sum(MMpartners %in% partners) ){
 			 message( paste("The following partners were not found in combination with the cell:", unlist(partners[which(MMpartners %in% partners)] ) ) )
 			 stop("\n")
 			 }
 	    peak_id_y_list <- MMpartners
 	    peak_id_x <- xTF_cell_tissue_name
   }
   
  if (tfIDs){
      partners_MMID = as.list(unlist(partners_MMID))
      peak_id_y_list = partners_MMID ; peak_id_x = mainTF_MMID
    }

    if( length(peak_id_y_list) < 2) { stop("The list of partner TF should contain 2 or more elements.") }

    Motif <- TFregulomeR::intersectPeakMatrix(peak_id_x = peak_id_x, motif_only_for_id_x = TRUE, peak_id_y = peak_id_y_list, motif_only_for_id_y = TRUE)

    motif_length <- dim(Motif[1,1][[1]]@MethMotif_x@MMBetaScore)[2] # get number of positions in motif

    if(as.numeric(forkPosition)>=motif_length){ stop("The forkPosition is larger than the motif length.") }

    if(flipMatrix==TRUE){
      for( i in 1:length(peak_id_y_list) ){
        Motif[1,i][[1]]@MethMotif_x@MMBetaScore <- Motif[1,i][[1]]@MethMotif_x@MMBetaScore[,motif_length:1] # reverse meth info
        colnames(Motif[1,i][[1]]@MethMotif_x@MMBetaScore) <- rev(colnames(Motif[1,i][[1]]@MethMotif_x@MMBetaScore)) # reverse meth info colnames
        Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[motif_length:1,] # reverse
        tmp_matrix <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix
        tmp_matrix[,'A'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'T'] # complementary
        tmp_matrix[,'T'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'A'] # complementary
        tmp_matrix[,'C'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'G'] # complementary
        tmp_matrix[,'G'] <- Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix[,'C'] # complementary
        Motif[1,i][[1]]@MethMotif_x@MMmotif@motif_matrix <- tmp_matrix

      }
      forkPosition <- which(motif_length:1 %in% forkPosition)-1 # adjust fork position
    }

    TFregulomeDataovlaplist <- list()
    TFregulomeDataBetalist <- list()
    TFregulomeDatamatrixlist <- list()
    TFregulomeDatanPeaks <- vector()
    TFregulomeDatanSites <- vector()

    for( i in 1:length(peak_id_y_list) ){
      y <- Motif[1,i]
      TFregulomeDataovlaplist[[i]] = y[[1]]@overlap_percentage_x
      TFregulomeDataBetalist[[i]] = y[[1]]@MethMotif_x@MMBetaScore
      TFregulomeDatamatrixlist[[i]] = y[[1]]@MethMotif_x@MMmotif@motif_matrix
      TFregulomeDatanPeaks[i] <- y[[1]]@MethMotif_x@MMmotif@nPeaks
      TFregulomeDatanSites[i] <- y[[1]]@MethMotif_x@MMmotif@nsites
      }
  
  FPWM <- new("FPWMClassObj")
  FPWM <- updateFPWMClassObj(FPWM, id = peak_id_y_list,
                                 nSites = TFregulomeDatanSites,
                                 nPeaks = TFregulomeDatanPeaks,
                                 matrix=TFregulomeDatamatrixlist,
                                 betalevel = TFregulomeDataBetalist,
                                 score= TFregulomeDataovlaplist,
                                 forkPosition = forkPosition)
  FPWM <- MatrixAdder( FPWM, forkPosition, probabilityMatrix)
  FPWM <- BetaAdder(FPWM, forkPosition)
  FPWM <- ConvertToFTRANSFAC(FPWM, probabilityMatrix, scaleFrequencyCounts)
  
  for (i in c(1:length(FPWM@betalevel))) {
    X<-FPWM@betalevel[[i]]
    FPWM@betalevel[[i]] <- X[,(forkPosition+1):ncol(FPWM@betalevel[[i]])]
  }
  FPWM@xid <- peak_id_x
  FPWM <- ModifyBetaFormat(FPWM)


  # order object based on overlapping score
  FPWM_tmp <- FPWM
  objOrder <- order(unlist(FPWM@score),decreasing=TRUE)
  FPWM_tmp@betalevel <- FPWM@betalevel[objOrder]
  FPWM_tmp@id <- FPWM@id[objOrder]
  FPWM_tmp@matrix <- FPWM@matrix[objOrder]
  FPWM_tmp@score <- FPWM@score[objOrder]

  FPWMPO <- FPWM@forked$PO
  from <- min(FPWMPO[duplicated(FPWMPO)])
  to <- max(FPWMPO)
  ix <- cbind(which(FPWMPO %in% from) , which(FPWMPO %in% to))

  for ( jx in 1:length(objOrder) ){
  from_row_ix <- ix[objOrder[jx],1] : ix[objOrder[jx],2]
  to_row_ix <- ix[jx,1] : ix[jx,2]
  FPWM_tmp@forked[to_row_ix ,2:5] <-  FPWM@forked[ from_row_ix, 2:5]
  }
  
  # add colnames
   FPWM <- FPWM_tmp
   colnames(FPWM@forked) <- c("PO","A","C","G","T")
   colnames(FPWM@parentbeta) <- c("PO","number","meth")
   for( i in 1:length(FPWM@betalevel) ) { colnames(FPWM@betalevel[[i]]) <- c("PO","number","meth")  }
   
return(FPWM)
}

MatrixAdder <- function( fpwmObject, forkPosition, probabilityMatrix)
{ 
  nSites <- fpwmObject@nSites
  if(probabilityMatrix==TRUE){
    total_nSites <- sum(nSites)
    S <- (fpwmObject@matrix[[1]][1:forkPosition,]) * (nSites[1]/total_nSites)
    for ( i in 2:length(fpwmObject@matrix) ) { S <- S + (fpwmObject@matrix[[i]][1:forkPosition,] * (nSites[i]/total_nSites) ) }
  }

  if(probabilityMatrix==FALSE){
    S <- round(fpwmObject@matrix[[1]][1:forkPosition,] * nSites[1])
    for ( i in 2:length(fpwmObject@matrix) ) { S <- S + round(fpwmObject@matrix[[i]][1:forkPosition,] * nSites[i])  }
  }

  fpwmObject@parentmatrix <- S
  return(fpwmObject)
}

BetaAdder <- function( fpwmObject, forkPosition)
{
  S <- fpwmObject@betalevel[[1]][,1:forkPosition]
  for ( i in 2:length(fpwmObject@id) ) {
      S <- S + fpwmObject@betalevel[[i]][,1:forkPosition]
      }
  fpwmObject@parentbeta <- S
  return(fpwmObject)
}

ConvertToFTRANSFAC <- function(fpwmObject, probabilityMatrix, scaleFrequencyCounts)
{ 
  Cnumber = length(fpwmObject@matrix)
  RowNum = nrow(fpwmObject@matrix[[1]])
  forkPosition = fpwmObject@forkPosition
  R = rep(c((forkPosition + 1):RowNum), times = Cnumber)
  Step = (RowNum - forkPosition)
  DF <- data.frame(colnames(c("PO","A","C","G","T")))
  DF[1:forkPosition,"PO"]=c(1:forkPosition)
  DF[(forkPosition+1):(length(R)+forkPosition),"PO"] <- R
  DF[1:forkPosition,2:5]=fpwmObject@parentmatrix
  c = 1

  if(probabilityMatrix == TRUE){
    for(i in seq(forkPosition+1, dim.data.frame(x = DF)[1], Step)){
      DF[i:((Step+i)-1),2:5] <- fpwmObject@matrix[[c]][(forkPosition+1):RowNum,]
      c <- c+1
    }
  }

  if(probabilityMatrix == FALSE){
    for(i in seq(forkPosition+1, dim.data.frame(x = DF)[1], Step)){

      if(scaleFrequencyCounts == FALSE){
        DF[i:((Step+i)-1),2:5] <- round(fpwmObject@matrix[[c]][(forkPosition+1):RowNum,] * fpwmObject@nSites[c])
      }

      if(scaleFrequencyCounts == TRUE){
        scaleNumber <- sum(DF[1,2:5])
        tmp_scaled_mat <- round(fpwmObject@matrix[[c]][(forkPosition+1):RowNum,] * scaleNumber)
        # check the round didn't add an extra count
        if( sum(rowSums(tmp_scaled_mat) != scaleNumber) ){
          for(ifix in which(rowSums(tmp_scaled_mat) != scaleNumber) ){
            diff_val <- sum(tmp_scaled_mat[ifix,]) - scaleNumber
            tmp_scaled_mat[ ifix, which.max(tmp_scaled_mat[ifix,]) ] <- tmp_scaled_mat[ ifix, which.max(tmp_scaled_mat[ifix,]) ] - diff_val
          }


        }

        DF[i:((Step+i)-1),2:5] <- tmp_scaled_mat
      }
      c <- c+1
    }
  }

  fpwmObject@forked <- DF
  return(fpwmObject)
}

ModifyBetaFormat <- function(fpwmObject)
{
  BS1 <- fpwmObject@parentbeta
  BS1 <-
    cbind(c("beta score<10%", "beta score 10-90%", "beta score>90%"),
          BS1)
  BS1 <- rbind(c("position", c(1:(ncol(BS1) - 1))), BS1)


  M1 <- matrix(nrow = ((ncol(BS1) - 1) * 3), ncol = 3)


  pos1 <- rep(t(BS1[1, 2:ncol(BS1)]), times = 3)
  M1[, 1] <- pos1

  pos1 <- t(BS1[2, 2:ncol(BS1)])
  M1[1:(ncol(BS1) - 1), 2] <- pos1
  pos1 <- t(BS1[3, 2:ncol(BS1)])
  M1[ncol(BS1):((ncol(BS1) - 1) * 2), 2] <- pos1
  pos1 <- t(BS1[4, 2:ncol(BS1)])
  M1[(((ncol(BS1) - 1) * 2) + 1):((ncol(BS1) - 1) * 3), 2] <- pos1

  pos1 <- rep(t(BS1[2, 1]), times = (ncol(BS1) - 1))
  M1[1:(ncol(BS1) - 1), 3] <- pos1
  pos1 <- rep(t(BS1[3, 1]), times = (ncol(BS1) - 1))
  M1[ncol(BS1):((ncol(BS1) - 1) * 2), 3] <- pos1
  pos1 <- rep(t(BS1[4, 1]), times = (ncol(BS1) - 1))
  M1[(((ncol(BS1) - 1) * 2) + 1):((ncol(BS1) - 1) * 3), 3] <- pos1

  fpwmObject@parentbeta <- M1
  
  for (i in c(1:length(fpwmObject@betalevel))) {
    
    BS2 <- as.matrix(fpwmObject@betalevel[[i]])
    BS2 <-
      cbind(c("beta score<10%", "beta score 10-90%", "beta score>90%"),
            BS2)
    BS2 <- rbind(c("position", c(fpwmObject@forkPosition + 1:(ncol(
      BS2
    ) - 1))), BS2)
    M2 <- matrix(nrow = ((ncol(BS2) - 1) * 3), ncol = 3)
  
  
    pos2 <- rep(t(BS2[1, 2:ncol(BS2)]), times = 3)
    M2[, 1] <- pos2
    
    pos2 <- t(BS2[2, 2:ncol(BS2)])
    M2[1:(ncol(BS2) - 1), 2] <- pos2
    pos2 <- t(BS2[3, 2:ncol(BS2)])
    M2[ncol(BS2):((ncol(BS2) - 1) * 2), 2] <- pos2
    pos2 <- t(BS2[4, 2:ncol(BS2)])
    M2[(((ncol(BS2) - 1) * 2) + 1):((ncol(BS2) - 1) * 3), 2] <- pos2
    
    pos2 <- rep(t(BS2[2, 1]), times = (ncol(BS2) - 1))
    M2[1:(ncol(BS2) - 1), 3] <- pos2
    pos2 <- rep(t(BS2[3, 1]), times = (ncol(BS2) - 1))
    M2[ncol(BS2):((ncol(BS2) - 1) * 2), 3] <- pos2
    pos2 <- rep(t(BS2[4, 1]), times = (ncol(BS2) - 1))
    M2[(((ncol(BS2) - 1) * 2) + 1):((ncol(BS2) - 1) * 3), 3] <- pos2
    
    fpwmObject@betalevel[[i]] <- M2}
  
  return(fpwmObject)
}