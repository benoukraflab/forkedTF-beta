#' Writes the FPWM in transfac or FPWMtransfac format
#'
#' This function saves the FPWM in transfac format.
#' @param FPWM FPWM object.
#' @param format  [character] the "transfac" option will output a standard transfact matrix per binding partner. The "FPWMtransfac" prints a single matrix with all the binding partners in FPWM format.
#' @param fileName [character] name of the file where the FPWM is going to be written.
#' @return A plain text file with the FPWM in transfac or FPWMtransfac format.
#' @keywords write.FPWM
#' @export
#' @examples
#' fpwm <- createFPWM(mainTF ="CEBPB",partners = c("ATF4","ATF7","ATF3"), cell = "K562", forkPosition = 5)
#' write.FPWM(FPWM = fpwm, format = "transfac", fileName = "FPWM.transfact" )
#' write.FPWM(FPWM = fpwm, format = "FPWMtransfac", fileName = "FPWM.FPWMtransfac" )

write.FPWM <- function( FPWM = NULL,
							format = "transfac",
							fileName = NULL
                           )
{

	if( is.null(FPWM) ){
        stop("Please provide an FPWM object: write.TRANSFAC( FPWM = yourFPWM, ...)  ") ; }

	if( format!="transfac" & format!="FPWMtransfac" ){
        stop("formats supported : transfac or FPWMtransfac ") ; }

    # Get format from forked
    diff_rowSums <- length(unique(rowSums(FPWM@forked[,2:5])))
    first_rowsum <- rowSums(FPWM@forked[,2:5])[1]
    matrix_format <- "CC Matrix"
    if(diff_rowSums != 1 ){ matrix_format <- "CC Count matrix" }else{
    	if(first_rowsum > 1){ matrix_format <- "CC Scale Count matrix" }
    	if(first_rowsum == 1){ matrix_format <- "CC Probability matrix" }
    }

	if( format == "FPWMtransfac"){
		ynames <- paste(unlist(FPWM@id),collapse="_&_")

		if( is.null(fileName) ){
			fileName <- paste0(FPWM@xid,"_+_",ynames,".FPWMtransfac")
		}

		fileConn <- file(fileName)

		transfac_vector <- c()
		transfac_vector <- c( transfac_vector, paste0("AC ",FPWM@xid,"_+_",ynames) )
		transfac_vector <- c( transfac_vector, "XX" )
		transfac_vector <- c( transfac_vector, paste0("parentLogo : ",FPWM@xid) )
		transfac_vector <- c( transfac_vector, paste0("leafLogos : ", paste(unlist(FPWM@id),collapse=",")) )
		transfac_vector <- c( transfac_vector, paste0("overlappingScore : ", paste(unlist(FPWM@score),collapse=",") ) )
		transfac_vector <- c( transfac_vector, paste0("numberOfMotifs : ", paste(unlist(FPWM@nSites),collapse=",") ) )
		transfac_vector <- c( transfac_vector, paste0("numberOfOverlappingPeaks : ", paste(unlist(FPWM@nPeaks),collapse=",") ) )
		transfac_vector <- c( transfac_vector, paste0("forkPosition : ",FPWM@forkPosition) )
		transfac_vector <- c( transfac_vector, "XX", paste("P0","A","C","G","T",sep="\t")  )

		for ( jx in 1:dim(FPWM@forked)[1] ){
				transfac_vector <- c( transfac_vector, paste(FPWM@forked[jx,],collapse="\t") )
			}
		transfac_vector <- c( transfac_vector, "XX", "CC FPWMtransfac format from FPWM",matrix_format, "XX" , "//"  )

		writeLines(transfac_vector, fileConn)
		close(fileConn)
		message(paste("FPWM saved as FPWMtransfac file:",fileName))
	}


	if( format == "transfac"){
		if( is.null(fileName) ){
			ynames <- paste(unlist(FPWM@id),collapse="_&_")
			fileName <- paste0(FPWM@xid,"_+_",ynames,".transfac")
		}

		fileConn <- file(fileName)

		FPWMPO <- FPWM@forked$PO
 	 	from <- min(FPWMPO[duplicated(FPWMPO)])
 		to <- max(FPWMPO)
 		ix_table <- cbind(which(FPWMPO %in% from) , which(FPWMPO %in% to))

		transfac_vector <- c()

		for( ix in 1:length(FPWM@id) ){
			transfac_vector <- c( transfac_vector, paste0("AC ",FPWM@xid,"_+_",FPWM@id[ix]), "XX" )
			transfac_vector <- c( transfac_vector, paste0("ID ",FPWM@xid,"_+_",FPWM@id[ix]), "XX" )
			transfac_vector <- c( transfac_vector, paste0("DE ",FPWM@xid,"_+_",FPWM@id[ix], " ; from MethMotif")  )
			transfac_vector <- c( transfac_vector, paste("P0","A","C","G","T",sep="\t")  )

			for ( jx in 1:as.numeric(FPWM@forkPosition) ){
				transfac_vector <- c( transfac_vector, paste(FPWM@forked[jx,],collapse="\t") )
			}

			ppForked <- FPWM@forked[ ix_table[ix,1]:ix_table[ix,2] , ]
			for ( jx in 1:dim(ppForked)[1] ){
				transfac_vector <- c( transfac_vector, paste(ppForked[jx,],collapse="\t") )
			}

			transfac_vector <- c( transfac_vector, c("XX","CC program: forkedTF")  )
			transfac_vector <- c( transfac_vector, paste0("CC numberOfMotifs: ",FPWM@nSites[ix])  )
			transfac_vector <- c( transfac_vector, paste0("CC numberOfOverlappingPeaks: ",FPWM@nPeaks[ix])  )
			transfac_vector <- c( transfac_vector, matrix_format  )
			transfac_vector <- c( transfac_vector, c("XX","//")  )
		}

		writeLines(transfac_vector, fileConn)
		close(fileConn)
		message(paste("FPWM saved as transfac file:",fileName))
	}
}
