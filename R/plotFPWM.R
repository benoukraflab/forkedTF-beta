#' A function to plot the Forked Position Weight Matrix
#'
#' This function plots FPWMs.
#' @param FPWM [FPWM S4 object].
#' @param Methylation [logical] If it set on TRUE, methylation level will be drawn as barplot. Default: TRUE.
#' @param height [numeric] Height (inches) of the pdf output. Default: 12.
#' @param width [numeric] Width (inches) of the pdf output. Default: 7.
#' @param textSize [numeric] Size of the text on the plot. Default: 7.
#' @param scaleBarplot [logical] TRUE to scale the height of the barplot across all the forked motifs. FALSE will scale per-motif. Default: TRUE.
#' @param legend [numeric Vector] Coordinates where the methylation legend will be drawn. Write "none" to prevent the function from ploting the legend. Default: c(0.5,0.95)
#' @param pdfName [character] Name for the output pdf.
#' @examples
#' fpwm <- createFPWM(mainTF ="CEBPB",partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"), cell = "K562", forkPosition = 5)
#' plotFPWM(fpwm, pdfName="fpwm_plot.pdf")
#' @return Plots the FPWM into a PDF file.
#' @export
plotFPWM <- function(FPWM, 
                     Methylation = TRUE, 
                     height = 12, 
                     width = 7,
                     textSize = 7,
                     scaleBarplot = TRUE,
                     legend = c(0.5,0.95),
                     pdfName = "plotFPWM.pdf")
{
 # Making the text size uniform across geom and external theme
  geom_text_size = textSize # geom text size
  text_size = (14/5) * geom_text_size # external theme size


  # getting number of partners 
  nrow <- length(FPWM@id)

if(scaleBarplot == TRUE & Methylation == FALSE){ stop("If scaleBarplot == TRUE, Methylation should be TRUE") }

if(nrow < 2 | nrow > 9){ stop("Function implemented for up to 9 interacting partners.") }

if(scaleBarplot == TRUE & Methylation == TRUE){
  y_ceiling <- lapply(FPWM@betalevel, function(x){
    max(aggregate(as.numeric(x[,2]) , by = list(Category = x[,1]), FUN = sum)[,2])
  })

  y_ceiling_parent <- aggregate(as.numeric(FPWM@parentbeta[,2]), by = list(Category = FPWM@parentbeta[,1]), FUN = sum)[,2]

  y_ceiling <- max(y_ceiling_parent,unlist(y_ceiling)) * 1.5 #y_celing holds the max height for all Cs
}

  
  # Extract the matrix position
  FPWMPO <- FPWM@forked$PO
  from <- min(FPWMPO[duplicated(FPWMPO)])
  to <- max(FPWMPO)
  ix <- cbind(which(FPWMPO %in% from) , which(FPWMPO %in% to))
   
   if (Methylation) {
    #plot beta score
       barplot_color <- c("darkorange1","darkgreen", "dodgerblue1")
       plot_beta_score <- as.data.frame(FPWM@parentbeta)                               
       plot_beta_score$meth <- relevel(plot_beta_score$meth, "beta score 10-90%")
       plot_beta_score$meth <- relevel(plot_beta_score$meth, "beta score>90%")
   
                              
       sum_of_pos <- aggregate(as.numeric(as.character(plot_beta_score$number)),
                              by=list(pos=as.numeric(as.character(plot_beta_score$PO))),
                              FUN=sum)
       sum_of_pos <- sum_of_pos[order(sum_of_pos$pos),]
                              
       colnames(sum_of_pos) <- c("pos", "sum")
       p_ylim <- max(aggregate(as.numeric(FPWM@parentbeta[,2]), by=list(Category=FPWM@parentbeta[,1]), FUN=sum)[,2]) * 1.5

       if(scaleBarplot == TRUE & Methylation == TRUE){
        p_ylim <- y_ceiling
       }

       databox <- plot_beta_score[order(plot_beta_score$meth, decreasing = FALSE),]
       databox$PO <- factor(databox$PO,levels(databox$PO)[order(as.numeric(levels(databox$PO)))])

       legend_meth <- legend

       parentLogo_meth <- ggplot(data = databox,
                   aes(x = PO, y = as.numeric(as.character(number)),
                       fill = meth)) +
        geom_bar(colour="black", stat="identity") +
        scale_fill_manual(values = barplot_color) + ylim(0, p_ylim) +
        theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
              axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "pt"),
              legend.position=legend_meth, 
              text = element_text(size = text_size), legend.title = element_blank(), legend.background = element_blank(),
              legend.box.background = element_rect(colour = "white"),legend.text=element_text(size=22), legend.key.size = unit(1,"line") ) +
        stat_summary(fun = sum, aes(label = stat(sum_of_pos$sum), group = PO), geom = "text",vjust = -0.5,size=geom_text_size)
    } else {
      parentLogo_meth <- ggplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank())
    }
   
    #plot parent logo
    dna_seq = FPWM@forked[1:(from-1), 2:5]
    dna_seq = t(dna_seq)
    logo_method <- "bits"

    parentLogo_motif<- ggplot() + geom_logo(data = dna_seq, method = logo_method, seq_type = "dna") + xlab(FPWM@xid)+
        theme(axis.title.y=element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit =  "pt"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),text = element_text(size = text_size))+scale_y_continuous(limits=c(0,2),breaks=c(0,1,2))
              
    parent_logo <- arrangeGrob(parentLogo_meth,parentLogo_motif, nrow=2)
  
    plotLogo_list <- list()
    for( i in 1:length(FPWM@id) ){
	
        if (Methylation) {
             #plot beta score
             plot_beta_score <- as.data.frame(FPWM@betalevel[[i]])                               
             plot_beta_score$meth <- relevel(plot_beta_score$meth, "beta score 10-90%")
             plot_beta_score$meth <- relevel(plot_beta_score$meth, "beta score>90%")
   
             sum_of_pos <- aggregate(as.numeric(as.character(plot_beta_score$number)),
                              by=list(pos=as.numeric(as.character(plot_beta_score$PO))),
                              FUN=sum)
             sum_of_pos <- sum_of_pos[order(sum_of_pos$pos),]
                              
             colnames(sum_of_pos) <- c("pos", "sum")
             
           p_ylim <- max(aggregate(as.numeric(FPWM@betalevel[[i]][,2]) , by = list(Category = FPWM@betalevel[[i]][,1]), FUN = sum)[,2]) * 1.5
           if(scaleBarplot == TRUE & Methylation == TRUE){
               p_ylim <- y_ceiling
           }
           plot_beta_score_reorder <- plot_beta_score[order(plot_beta_score$meth, decreasing = FALSE),]
           plot_beta_score_reorder$PO <- factor(plot_beta_score_reorder$PO, 
                                             level = sort( unique( as.numeric ( as.character ( plot_beta_score_reorder$PO ) ) ) ) #same as factor reorder for databox
                                       )
             
             p1j <- ggplot(data = plot_beta_score_reorder,
                   aes(x=PO,y=as.numeric(as.character(number)),
                       fill=meth)) +
                   geom_bar(colour="black", stat="identity") +
                   scale_fill_manual(values = barplot_color) + ylim(0, p_ylim) +
                   theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.y=element_blank(),
                   axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
                   plot.margin = margin(t = 0, r = 0, b = 0, l = 5, unit = "pt"),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(),legend.position="none",text = element_text(size = text_size)) +
                   stat_summary(fun = sum, aes(label = stat(sum_of_pos$sum), group = PO), geom = "text",vjust = -0.5,size=geom_text_size)
        } else {
             p1j <- ggplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank())
        }
    
     #plot child logo
    dna_seq = FPWM@forked[ix[i,1]:ix[i,2], 2:5]
    dna_seq = t(dna_seq)
    
    p2j<- ggplot() + geom_logo(data = dna_seq, method = logo_method, seq_type = "dna") + xlab(FPWM@id[[i]]) +
        theme(axis.title.y=element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit =  "pt"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),text = element_text(size = text_size) ) +
              scale_y_continuous(limits=c(0,2),breaks=c(0,1,2)) 
              
    p2j$scales$scales[[1]] <- scale_x_continuous(breaks = 1:length(ix[1,1]:ix[1,2]),
   											 labels=as.character(ix[1,1]:ix[1,2]))
       
              p_j <- arrangeGrob(p1j,p2j,nrow=2)
              
        plotLogo_list[[ FPWM@id[[i]] ]] <- p_j
  }
  
    forkLogo_Plot <- do.call("arrangeGrob", c(plotLogo_list, ncol=1))

    parentLogo_position <-  ceiling( length(forkLogo_Plot) / 2 ) 
    
    parentLogo_list <- plotLogo_list
    
    for(i in 1:length(parentLogo_list)) {
        parentLogo_list[[i]] = ggplot() + theme(axis.title.y=element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit =  "pt"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(),text = element_text(size = text_size))+scale_y_continuous(breaks=c(0,1,2))
    }
    
    parentLogo_list[[parentLogo_position]] <- parent_logo
    
    parentLogo_Plot <- do.call("arrangeGrob", c(parentLogo_list, ncol=1))
    
    
    
    arrow_coord <- data.frame(x = rep(0,length(parentLogo_list)), 
                              y = rep(parentLogo_position,length(parentLogo_list)), 
                              vx = rep(1,length(parentLogo_list)), 
                              vy= 1:length(parentLogo_list) - ( (length(parentLogo_list) / 2) + 0.5 ) )

    hardcoded_y_list <- list(
      c(1,5.4),
      c(.55,3.5,6.4),
      c(0.3,2.5,4.7,6.9),
      c(0.15,1.92,3.69,5.46,7.23),
      c(0.07,1.54,3.01,4.48,5.95,7.42),
      c(0.01,1.27,2.53,3.79,5.05,6.31,7.57),
      c(0.001,1.05,2.15,3.25,4.35,5.45,6.55,7.65), 
      c(0,.89,1.87,2.85,3.83,4.81,5.79,6.77,7.75) )

      test_y <- hardcoded_y_list[[nrow-1]]
      ix_array <- nrow:1
      ix_array <- ix_array[ceiling(nrow/2)]
      arrow_coord$y<-rep(test_y[ix_array],nrow)
      arrow_coord$vy <- test_y-test_y[ix_array]

    percentage_score <- round(unlist(FPWM@score), digits=1)
    percentage_score <- rev(paste(percentage_score,"%",sep=""))

    arrow_plot <- ggplot() + 
        geom_segment(data=arrow_coord, mapping=aes(x=x, y=y, xend=x+vx, yend=y+vy), arrow=arrow(), size=1, color=alpha("darkgrey",1)) + 
        geom_point(data=arrow_coord, mapping=aes(x=x, y=y), size=1, shape=21, fill="white") + ylim(0, 8) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.title.x=element_blank(), axis.text.y=element_blank(),axis.title.y=element_blank(),
        axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit =  "pt"),text = element_text(size = text_size)) + annotate("text", x=rep(.95,length(percentage_score)), 
                                                                                           y = test_y , 
                                                                                           label=percentage_score , size = geom_text_size)
    
    arrow_plot <- do.call("arrangeGrob", c(list(arrow_plot), ncol=1))

    # plotLogo_list
    # parent_logo
    g<-arrangeGrob(parentLogo_Plot,arrow_plot,forkLogo_Plot,
                     layout_matrix = rbind(c(1,2,3))
                     )

    save_plot(file=pdfName, g,nrow=nrow,base_width=11)

}

