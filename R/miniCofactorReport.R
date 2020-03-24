#' Produces a miniCofactorReport
#'
#' This function allows you to get a PDF report of top cofactors along with DNA methylation information for a motif of a TF.
#' @param TF [character] Main TF of interest.
#' @param cell [character] Cell of interest.
#' @param cobinding_threshold [numeric] Only the co-factors with co-binding percentages more than this threshold value will be reported. By default the threshold is 0.05.
#' @param Methylation [logical] TRUE to retrieve cytosine methylation information.
#' @param includeMotifOnly [logical] TRUE if you wish to include only peaks that contain the known binding motif.
#' @param height [numeric] Height in inch for the plot.
#' @param width [numeric] Width in inch for the plot.
#' @param pdfName [character] Name of the pdf to be saved.
#' @return A PDF file with the cofactorReport.
#' @keywords cofactorReport
#' @export
#' @examples
#' miniCofactorReport(TF = "CEBPB",cell = "K562")

miniCofactorReport <- function(TF,
              cell,
              cobinding_threshold=0.05,
              Methylation = TRUE,
              includeMotifOnly = TRUE,
              height = 12,
              width = 7,
              pdfName = NULL) {

   # Setting the max number of binding partners to 10
   NumberofTop = 10

   # check other input arguments
  if (cobinding_threshold <= 0 || cobinding_threshold > 1){
        stop("'cobinding_threshold' should be >0 but <=1.")
    }

              
  if( !(class(TF) == "character" | class(TF) == "data.frame") ){ 
        stop("Invalid input for 'TF'. Either write the name of a TF or provide a bedfile in a data.frame format") ; } 

  if( class(TF) == "character" ){
      TF_cell_tissue_name <- dataBrowser(tf = TF, cell_tissue_name = cell) ;
      
      if(is.null(TF_cell_tissue_name)){ stop("Please check the spelling of your TF or cell. This is case-sensitive.") }
  
      if( dim(TF_cell_tissue_name)[1] != 1 ){ stop("More than one record for the combination of TF and cell tissue") ; } 

        TF_cell_peaks <- loadPeaks(id = TF_cell_tissue_name$ID[1], includeMotifOnly = includeMotifOnly) ;
    
        cell_TFBS <- dataBrowser(cell_tissue_name = cell) ;
        message("#################################")
        message("This might take a few minutes!!!")
        message("#################################")

        intersectMatrix_forCofactorReport <- intersectPeakMatrix(user_peak_list_x = list(TF_cell_peaks),
                                              user_peak_x_id = TF_cell_tissue_name$ID[1], 
                                              peak_id_y = cell_TFBS$ID, 
                                              motif_only_for_id_y = includeMotifOnly, 
                                              methylation_profile_in_narrow_region = Methylation) ;
                                              
        FPWMcofactorReport_ui(intersectPeakMatrix = intersectMatrix_forCofactorReport,
                    top_num = NumberofTop,
                    cobinding_threshold = cobinding_threshold,
                    height = height,
                    width = width,
                    pdfName = pdfName
                    ) ;
    }

  if( class(TF) == "data.frame" ){
        cell_TFBS <- dataBrowser(cell_tissue_name = cell) ;
        
        used_TF_list <- deparse(substitute(TF))
        
        message("#################################")
        message("This might take a few minutes!!!")
        message("#################################")
        
        intersectMatrix_forCofactorReport <- intersectPeakMatrix(user_peak_list_x = list(TF_cell_peaks),
                                              user_peak_x_id = used_TF_list , 
                                              peak_id_y = cell_TFBS$ID, 
                                              motif_only_for_id_y = includeMotifOnly, 
                                              methylation_profile_in_narrow_region = TRUE) ;
                                              
        FPWMcofactorReport_ui(intersectPeakMatrix = intersectMatrix_forCofactorReport,
                    top_num = NumberofTop,
                    cobinding_threshold = cobinding_threshold,
                    height = height,
                    width = width,
                    pdfName = pdfName
                    ) ;
  }
}


FPWMcofactorReport_ui <- function(intersectPeakMatrix,
                           top_num = 10,
                           cobinding_threshold=0.05,
        				   height = height,
	                       width = width,
        				   pdfName = pdfName)
{
    # check input arguments
    if (missing(intersectPeakMatrix))
    {
        stop("Please provide output of 'intersectPeakMatrix()' using 'intersectPeakMatrix ='!")
    }
    # check the validity of input intersectPeakMatrix
    if (class(intersectPeakMatrix[1,1][[1]])[1] != "IntersectPeakMatrix")
    {
        stop("The input 'intersectPeakMatrix' is not valid. Please use the output of function 'intersectPeakMatrix()'")
    }


    # start reporting
    message("Start miniCofactorReport ...")
    message(paste0("... The maximum number of cofactors to be reported is ", top_num))
    message(paste0("... The minimum percent of co-binding peaks for a cofactor is ", cobinding_threshold*100,"%"))

    intersectPeakMatrix_i <- intersectPeakMatrix
    id_i <- rownames(intersectPeakMatrix_i) # mainTF id
    is_from_TFregulomeR <- intersectPeakMatrix_i[1,1][[1]]@isxTFregulomeID # check if the mm id is correct
        
		message(paste0("... ... Start reporting peak id '",id_i,"' ..."))
  	suppressMessages(intersectPeakMatrix_res_i <- intersectPeakMatrixResult(intersectPeakMatrix_i,
                                                                                return_intersection_matrix = TRUE,
                                                                                angle_of_matrix = "x"))
    intersect_matrix_i <- intersectPeakMatrix_res_i$intersection_matrix
    intersect_matrix_t_i <- as.data.frame(t(intersect_matrix_i))
    intersect_matrix_filter_i <- intersect_matrix_t_i[(intersect_matrix_t_i[,1] >= cobinding_threshold*100),,drop=FALSE]

    if (nrow(intersect_matrix_filter_i) < 1)
    {
        stop("No overlap between your TF of interest and other TFs was found with the specified paramaters.")
    }

    intersect_matrix_order_i <- intersect_matrix_filter_i[order(intersect_matrix_filter_i[,1]),,drop = FALSE]

    # cobinding barplot
    intersect_matrix_heatmap_i <- as.data.frame(matrix(nrow=nrow(intersect_matrix_order_i), ncol = 5))

    colnames(intersect_matrix_heatmap_i) <- c("x","new_x","y","new_y","value")
    intersect_matrix_heatmap_i$y <- rownames(intersect_matrix_order_i)
    intersect_matrix_heatmap_i$new_y <- paste(unlist(lapply(rownames(intersect_matrix_order_i),
                                                          function(x) tail(unlist(strsplit(x,split = "_")),1))),
                                                   rev(rownames(intersect_matrix_heatmap_i)), sep = "-")

    intersect_matrix_heatmap_i$x <- colnames(intersect_matrix_order_i)
    intersect_matrix_heatmap_i$new_x <- unlist(lapply(colnames(intersect_matrix_order_i),
                                                          function(x) tail(unlist(strsplit(x,split = "_")),1)))

    intersect_matrix_heatmap_i$value <- intersect_matrix_order_i[,1]
    intersect_matrix_heatmap_i$new_y <- factor(intersect_matrix_heatmap_i$new_y,
                                                   levels = as.character(intersect_matrix_heatmap_i$new_y))

    cobinding_ylabel <- as.character(intersect_matrix_heatmap_i$new_y[rev(seq(1,nrow(intersect_matrix_heatmap_i),1))])
    cobinding_ylabel_new <- paste0(as.character(intersect_matrix_heatmap_i$new_x),
                                       "\n+\n",cobinding_ylabel)

    colors_cobinding <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))(11)
    cobinding_ylabel_new <- rev(gsub("\\-\\d+$","",cobinding_ylabel_new,perl=TRUE))

    topIX <- length(cobinding_ylabel_new)
    cobinding_ylabel_new[ topIX ] <-  paste0("All\n",as.character(intersect_matrix_heatmap_i$new_x)[topIX]) # adding the label all

    p1 <- ggplot(intersect_matrix_heatmap_i, aes(x=new_y, 
                                       y=value,
                                       fill = value)) +
                     geom_bar(stat="identity")+ scale_fill_gradientn(colours=colors_cobinding,
                                             breaks=c(seq(0, 100, length=11)),
                                             limits=c(0,100)) + scale_y_reverse(position = "right") +  coord_flip() +
                     geom_text(aes(label=paste0(round(value,digits=1)," %")), vjust=-.2, angle =90) +
                     theme(panel.background = element_blank(),
                           plot.background = element_blank(),
                           axis.title=element_blank(),
                           axis.ticks = element_blank(),
                           plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
                           axis.text.x = element_text(size=10),axis.text.y = element_text(size=12),
                           legend.position = "none") + 
                     scale_x_discrete(labels=cobinding_ylabel_new) 

    # motifs
    motif_plot_list_p <- list()
    count <- 1
    for (j in order(-seq(1,nrow(intersect_matrix_heatmap_i),1)))
    {
       x_j <- intersect_matrix_heatmap_i$x[j]
       y_j <- intersect_matrix_heatmap_i$y[j]
       motif_j <- t(intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_x@MMmotif@motif_matrix)
       top <- 10
       bottom <- 0
       if (j==1)
       {
              bottom <- 10
              top <- 0
       }

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
 #  Edited from plotLogo.R in TFregulomeR-dev 2019-09-18 :
 
       if (is_from_TFregulomeR) { MM_object = intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_x } 
       if (!is_from_TFregulomeR) { MM_object = intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_y } 
   
       y_max <- 2
       logo_method <- "bits"

       # get beta score matrix and motif matrix
       MMBetaScore <- MM_object@MMBetaScore
       MMmotif <- MM_object@MMmotif
       motif_matrix <- t(MMmotif@motif_matrix)

       motif_length <- ncol(motif_matrix)

       if (!(is.na(MMBetaScore[1,1])))
       {
            # generate a dataframe for beta score plotting
              plot_beta_score <- matrix(rep(0,length(MMBetaScore)*3), ncol = 3)
              colnames(plot_beta_score) <- c("number","pos","meth")
              plot_beta_score[seq(1,motif_length,1),1] <- as.vector(MMBetaScore[3,])
              plot_beta_score[(motif_length+1):(2*motif_length),1] <- as.vector(MMBetaScore[2,])
              plot_beta_score[(2*motif_length+1):(3*motif_length),1] <- as.vector(MMBetaScore[1,])
              plot_beta_score[seq(1,motif_length,1),2] <- seq(1, motif_length, 1)
              plot_beta_score[(motif_length+1):(2*motif_length),2] <- seq(1, motif_length, 1)
              plot_beta_score[(2*motif_length+1):(3*motif_length),2] <- seq(1, motif_length,1)
              plot_beta_score[seq(1,motif_length,1), 3] <- "beta score>90%"
              plot_beta_score[(motif_length+1):(2*motif_length),3] <- "beta score 10-90%"
              plot_beta_score[(2*motif_length+1):(3*motif_length),3] <- "beta score<10%"
              plot_beta_score <- as.data.frame(plot_beta_score)

              # plot all, methylated or unmethylated
              # make levels in beta score plotting matrix
              plot_beta_score$meth <- factor(plot_beta_score$meth,levels = c("beta score>90%",  "beta score 10-90%","beta score<10%"))
              plot_beta_score$pos <- factor(plot_beta_score$pos, levels = seq(1,motif_length,1))
              ylim <- round(max(as.vector(apply(MMBetaScore,2,sum)))/1000+1)*1000+500
              barplot_color <- c("darkorange1","darkgreen", "dodgerblue1")
      
              sum_of_pos <- aggregate(as.numeric(as.character(plot_beta_score$number)),
                                      by=list(pos=plot_beta_score$pos),
                                      FUN=sum)
              colnames(sum_of_pos) <- c("pos", "sum")
              #plot beta score
              p1j <- ggplot(data = plot_beta_score[order(plot_beta_score$meth, decreasing = FALSE),],
                            aes(x=pos,y=as.numeric(as.character(number)),
                            fill=meth)) +
                            geom_bar(colour="black", stat="identity") +
                            scale_fill_manual(values = barplot_color) + ylim(0, ylim) +
                            theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
                            axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
                            plot.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "pt"),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(),legend.position="none") +
                            stat_summary(fun = sum, aes(label = stat(sum_of_pos$sum), group = pos), geom = "text",vjust = -0.5)
    }
    else
    {
              p1j <- ggplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank())
    }

              # motif logo position size
              #size xlab
    if (motif_length>40)
    {
              xlab_size <- 4
    }
    else
    {
              xlab_size <- -0.5*motif_length+24
    }
              #plot motif logo
    p2j <- ggplot() + geom_logo(data = motif_matrix, method = logo_method) +
    theme(axis.title.y=element_blank(),
    plot.margin = margin(t = 0, r = -6, b = -6, l = -8, unit =  "pt"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank())+scale_y_continuous(breaks=c(0,1,2))
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
    p_j <- arrangeGrob(p1j,p2j, nrow=2)
    motif_plot_list_p[[count]] <- p_j
    count <- count+1
    }
            p2 <- arrangeGrob(grobs=motif_plot_list_p,
                              nrow=nrow(intersect_matrix_heatmap_i))

        text1 <- textGrob("Co-binding (%)")
        text2 <- textGrob("Motif")


        pdf_i_name <- paste0(id_i,"_cofactor_minireport.pdf")
        
        if( !is.null(pdfName) ){ 
        	pdf_i_name <- pdfName
        	pdftest <- grep(".pdf$",x,perl=TRUE)
        	if( identical(pdftest,integer(0)) ){ 
        		pdf_i_name <- paste0(pdf_i_name,".pdf")
        	}
        }

        blank <- grid.rect(gp=gpar(col="white"))

        def_layout_matrix <- rbind(c(1,1,2),
                                   c(3,3,4), #1
                                   c(3,3,4), #2
                                   c(3,3,4), #3
                                   c(3,3,4), #4
                                   c(3,3,4), #5
                                   c(3,3,4), #6
                                   c(3,3,4), #7
                                   c(3,3,4), #8
                                   c(3,3,4), #9
                                   c(3,3,4), #10
                                   c(3,3,4), #11
                                   c(3,3,4)  #12
                                   )

        if( nrow(intersect_matrix_filter_i)<10 ){
            ixlayout <- ceiling((nrow(intersect_matrix_filter_i)*12)/10)+2
            def_layout_matrix[ixlayout:13,] = 5
          pdf(pdf_i_name,height=height,width=width)
          grid.arrange(text1,text2,p1, p2,blank,
                     layout_matrix = def_layout_matrix)
          dev.off()
          message(paste0("... ... ... miniCofactor report for id '", id_i,"' has been saved as ", pdf_i_name))

        }else{
          pdf(pdf_i_name,height=height,width=width)
          grid.arrange(text1,text2,p1, p2,
                     layout_matrix = def_layout_matrix)
          dev.off()
          message(paste0("... ... ... miniCofactor report for id '", id_i,"' has been saved as ", pdf_i_name))
      }
    
}

# thigs to fix: for GTRD matrix with no meth data.


						  
