#' Produces a miniCofactorReport
#'
#' This function allows you to get a PDF report of top cofactors along with DNA methylation information for a motif of a TF.
#' @param TF [character] Main TF of interest.
#' @param cell [character] Cell of interest.
#' @param filterBy [character] Threshold category to filter cobinding partners.
#' Currently supported are: "mapped.peaks.ratio,effect.size,p.significance,p.value,q.significance,q.value,e.significance,e.value and fraction."
#' @param threshold [numeric] Only the co-factors with co-binding percentages more than this threshold value will be reported. By default the threshold is 0.05.
#' @param Methylation [logical] TRUE to retrieve cytosine methylation information.
#' @param includeMotifOnly [logical] TRUE if you wish to include only peaks that contain the known binding motif.
#' @param height [numeric] Height in inch for the plot.
#' @param width [numeric] Width in inch for the plot.
#' @param pdfName [character] Name of the pdf to be saved.
#' @param server [character] server localtion to be linked, either 'sg' or 'ca'.
#' for 'Singapore' or 'Canada', respectively.
#'
#' @param universe A set of genomic regions that prevent shuffles
#' for occuring outside of it.
#' @param chromSizes A vector containing all the chromosome
#' lengths for the species in consideration.
#' @param shufflesNumber The number of shuffled genomic regions to be created for
#' theorical distribution (higher means more accurate).
#' @param shuffle_seed The random seed to be used for shuffling.
#' @param tail If "lower" then, probabilities are P[X > x],
#' if "higher", P[X <= x], if "both" then higher or lower is selected
#' depending on the number of overlaps vs the theorical mean.
#' @param pAdjust The method that will be used for correcting the p-values.
#' BH, BY and bonferroni are available.
#' @param byChrom Will the shuffles stay in the chromosome they originate (TRUE)
#' or can they be placed everywhere on the genome (FALSE)
#' @param included Represents the fraction of each regions that can
#' be outside of the universe.
#'
#' @return A PDF file with the cofactorReport.
#' @keywords cofactorReport
#' @export
#' @examples
#' miniCofactorReport(TF = "CEBPB",cell = "K562")

miniCofactorReport <- function(TF,
              cell,
              filterBy  = "mapped.peaks.ratio",
              threshold = 0,
              Methylation         = TRUE,
              includeMotifOnly    = TRUE,
              height = 12,
              width = 7,
              shufflesNumber = 100,
              shuffle_seed = 987,
              universe=NULL,
              tail = "lower",
              pAdjust = "BY",
              chromSizes = loadChromSizes("hg38"),
              byChrom = FALSE,
              included = 1,
              pdfName = NULL,
              server = "sg" ) {

  # Setting the max number of binding partners to 10
  NumberofTop = 10
  options(warn = -1)

  # Validate input filterBy Flag
  filterBy_flag <- switch(filterBy,"mapped.peaks.ratio" = 1,
                                   "effect.size"        = 1,
                                   "p.significance"     = 1,
                                   "p.value"            = 1,
                                   "q.significance"     = 1,
                                   "q.value"            = 1,
                                   "e.significance"     = 1,
                                   "e.value"            = 1,
                                   "fraction"           = 1)

  if( is.null(filterBy_flag) ){
    stop("Non-valid ",filterBy," parameter in 'filterBy' option. See available options.")
  }

  if( threshold < 0 ){
    stop("'threshold' should be >= 0.")
  }

  if( missing(cell) ) {
    stop("Unable to proceed, a valid cell-line name should be input to 'cell' parameter.")
  }

  if( !(class(TF) == "character" | class(TF) == "data.frame") ){
        stop("Invalid input for 'TF'. Either write the name of a TF or provide a bedfile in a data.frame format") ; }

  if( class(TF)   == "character" ){
      # Retrieve metadata information for the query TF in a given cell type/tissue from DB
      TF_cell_tissue_name <- dataBrowser(tf = TF, cell_tissue_name = cell, server = server)

      if(is.null(TF_cell_tissue_name)){ stop("Please check the spelling of your TF or cell. This is case-sensitive.") }

      if( dim(TF_cell_tissue_name)[1] != 1 ){ stop("More than one record for the combination of TF and cell tissue")  }

      # Retrieve peak BED summit information from DB for the fetched TF
      TF_cell_peaks <- loadPeaks(id = TF_cell_tissue_name$ID[1], includeMotifOnly = includeMotifOnly, server = server)

      # Convert peak BED to GRanges
      tf_query      <- TF_cell_peaks %>%
        dplyr::mutate( start = start - 99,
                       end   = end   + 100) %>%
        .convertDF_to_GRanges()

      # Get list of enriched co-factors
      enriched_cofactors <- .get_enriched_cofactors(tf_query,
                                                    cell,
                                                    filterBy  = filterBy,
                                                    includeMotifOnly = includeMotifOnly,
                                                    shufflesNumber = shufflesNumber,
                                                    shuffle_seed = shuffle_seed,
                                                    universe = universe,
                                                    tail = tail,
                                                    pAdjust = pAdjust,
                                                    chromSizes = chromSizes,
                                                    byChrom = byChrom,
                                                    included = included,
                                                    server = server)

      message("#################################")
      message("This might take a few minutes!!!")
      message("#################################")


      if( filterBy != "fraction") {
        intersectMatrix_forCofactorReport <-  TFregulomeR::intersectPeakMatrix(user_peak_list_x = list(TF_cell_peaks),
                                                                               #user_peak_x_id   = enriched_cofactors[["TF_cell_tissue_name"]]$ID[1],
                                                                               user_peak_x_id   = TF_cell_tissue_name$ID[1],
                                                                               peak_id_y = enriched_cofactors[["top_cell_TFBS"]]$ID,
                                                                               motif_only_for_id_y = includeMotifOnly,
                                                                               methylation_profile_in_narrow_region = Methylation,
                                                                               server = server) ;
        FPWMcofactorReport_ui(intersectPeakMatrix = intersectMatrix_forCofactorReport,
                              top_tf              = enriched_cofactors[["top_tf"]],
                              top_num             = NumberofTop,
                              filterBy            = filterBy,
                              threshold           = threshold,
                              height              = height,
                              width               = width,
                              pdfName             = pdfName) ;
      } else{
        intersectMatrix_forCofactorReport <-  TFregulomeR::intersectPeakMatrix(user_peak_list_x    = list(TF_cell_peaks),
                                                                               user_peak_x_id   = TF_cell_tissue_name$ID[1],
                                                                               #user_peak_x_id      = enriched_cofactors[["TF_cell_tissue_name"]]$ID[1],
                                                                               peak_id_y           = enriched_cofactors[["cell_TFBS"]]$ID,
                                                                               motif_only_for_id_y = includeMotifOnly,
                                                                               methylation_profile_in_narrow_region = Methylation,
                                                                               server = server) ;
        FPWMcofactorReport_ui(intersectPeakMatrix = intersectMatrix_forCofactorReport,
                              top_tf              = NULL,
                              top_num             = NumberofTop,
                              filterBy            = filterBy,
                              threshold           = threshold,
                              height              = height,
                              width               = width,
                              pdfName             = pdfName) ;
      }
      # Return matrix enrichment
      return( enriched_cofactors[["enrichment_df"]] )
  }

  if( class(TF) == "data.frame" ){
    # Get number of columns
    nb_cols <- TF %>%
      dim() %>%
      head(2)

    # Check number of columns
    if( nb_cols < 3){
      stop("Error in input data.frame coordinates, at least 3 columns should be provided.")
    }
    if( nb_cols > 3 ){
      message("Only the first 3 columns will be used.")
      TF <- TF %>%
        dplyr::select( seq(1,3) )
    }
    # Add column names
    colnames(TF) <- c("chr","start","end")

    # Convert column types
    TF <- TF %>% dplyr::mutate( chr   = chr   %>% as.character )
    TF <- TF %>% dplyr::mutate( start = start %>% as.numeric )
    TF <- TF %>% dplyr::mutate( end   = end   %>% as.numeric )

    # Add 'chr' prefix if missing in the first column.
    chr_records_to_update <- TF %>%
      dplyr::filter( !grepl("^chr",chr) ) %>%
      dim() %>%
      head(1)

    if( chr_records_to_update > 0 ){
      message("Detected",chr_records_to_update,"records with missing 'chr' prefix in chromosome column,",
              "the prefix will be added.")
      TF <- TF %>%
        dplyr::mutate(
          chr = replace(
            chr,
            !grepl("^chr", chr),
            paste0("chr", chr[!grepl("^chr", chr)]  )
          )
        )
    }

    # Test if start is equal to end
    nb_equal_coord <- TF %>%
      dplyr::filter( start == end) %>%
      dim() %>%
      head(1)

    if (nb_equal_coord > 0) {
      stop("Unable to proceed:",nb_equal_coord, "records have the same 'start' and 'end' value.",
           "Please be aware that BED format 'start' is a 0-based index and 'end' a 1-based index.")
    }

    # Test if start is bigger than end
    nb_low_coord <- TF %>%
      dplyr::filter( start > end) %>%
      dim() %>%
      head(1)

    if (nb_low_coord > 0){
      stop("Unable to proceed:",nb_low_coord, "records have a 'start' greater than 'end' value.",
           "Something is weird in your data.")
    }

    # Converto GRanges object
    tf_query <- .convertDF_to_GRanges(TF)

    enriched_cofactors <- .get_enriched_cofactors(tf_query = tf_query,
                                                  cell = cell,
                                                  filterBy  = filterBy,
                                                  includeMotifOnly = includeMotifOnly,
                                                  shufflesNumber = shufflesNumber,
                                                  shuffle_seed = shuffle_seed,
                                                  universe = universe,
                                                  tail = tail,
                                                  pAdjust = pAdjust,
                                                  chromSizes = chromSizes,
                                                  byChrom = byChrom,
                                                  included = included,
                                                  server = server)


    message("#################################")
    message("This might take a few minutes!!!")
    message("#################################")
    if( filterBy != "fraction"){
      intersectMatrix_forCofactorReport <-  TFregulomeR::intersectPeakMatrix(user_peak_list_x = list(TF),
                                             #user_peak_x_id   = enriched_cofactors[["TF_cell_tissue_name"]]$ID[1],
                                             peak_id_y = enriched_cofactors[["top_cell_TFBS"]]$ID,
                                             motif_only_for_id_y = includeMotifOnly,
                                             methylation_profile_in_narrow_region = Methylation,
                                             server = server);
      FPWMcofactorReport_ui(intersectPeakMatrix = intersectMatrix_forCofactorReport,
                            top_tf              = enriched_cofactors[["top_tf"]],
                            top_num             = NumberofTop,
                            filterBy            = filterBy,
                            threshold           = threshold,
                            height              = height,
                            width               = width,
                            pdfName             = pdfName) ;
    }else{
      intersectMatrix_forCofactorReport <-  TFregulomeR::intersectPeakMatrix(user_peak_list_x    = list(TF),
                                             #user_peak_x_id      = enriched_cofactors[["TF_cell_tissue_name"]]$ID[1],
                                             peak_id_y            = enriched_cofactors[["cell_TFBS"]]$ID,
                                             motif_only_for_id_y  = includeMotifOnly,
                                             methylation_profile_in_narrow_region = Methylation,
                                             server = server) ;
      FPWMcofactorReport_ui(intersectPeakMatrix = intersectMatrix_forCofactorReport,
                            top_tf              = NULL,
                            top_num             = NumberofTop,
                            filterBy            = filterBy,
                            threshold           = threshold,
                            height              = height,
                            width               = width,
                            pdfName             = pdfName) ;
    }

    # Return matrix enrichment
    return( enriched_cofactors[["enrichment_df"]] )
  }

}

# helpers -----------------------------------------------------------------

# Convert a data.frame into a GRanges object
# order of column must imperatively be:
# chr-start-end-id-score. I should later
# autodetect regex of columns to span the cases

.convertDF_to_GRanges <- function(bed_df){
  # Get col number
  nb_cols <- bed_df %>%
    dim() %>%
    tail(1)

  # Rename bed columns
  if( nb_cols < 3 ){
    message("Unable to convert to GRanges object, column number < 3.")
    return(NULL)
  } else if(nb_cols == 3) {
    colnames(bed_df) <- c("chr", "start", "end")
  } else if(nb_cols == 5 ) {
    colnames(bed_df) <- c("chr", "start", "end", "id", "score")
  } else {
    message("Unable to convert to GRanges object, not a valid number of columns.")
    return(NULL)
  }

  # Convert to GRanges object
  bed_df <- makeGRangesFromDataFrame(bed_df,
                                     keep.extra.columns      = TRUE,
                                     starts.in.df.are.0based = FALSE)
}

# Get enriched cofactors from catalog
.get_enriched_cofactors <- function(tf_query,
                                cell,
                                NumberofTop = 10,
                                filterBy  = "q.significance",
                                includeMotifOnly    = TRUE,
                                shufflesNumber = 100,
                                shuffle_seed = 987,
                                universe = NULL,
                                tail = "lower",
                                pAdjust = "BY",
                                chromSizes = loadChromSizes("hg38"),
                                byChrom = FALSE,
                                included = 1,
                                server = "sg" ) {
  # Convert 'fraction' to 'nb.overlaps' for compatibility with enrichment
  if(filterBy == "fraction"){
    filterBy <- "nb.overlaps"
  }
  # Retrieve metadata information for the query cell type/tissue (all TFs found) from DB
  cell_TFBS <- dataBrowser(cell_tissue_name = cell, server = server)

  # Retrieve peak BED summit information for the TF collection
  cell_catalog <- list()
  for ( i in seq(1, dim(cell_TFBS)[1]) ) {
    # Retrieve name for the list element
    tf_name <- cell_TFBS %>% dplyr::slice(i) %>% pull("TF")
    # Fetch peaks in DB and coverto to GRanges
    cell_catalog[[ tf_name ]] <- TFregulomeR::loadPeaks(id = dplyr::slice(cell_TFBS,i) %>% dplyr::pull("ID"),
                                                        includeMotifOnly = includeMotifOnly, server = server) %>%
      dplyr::mutate(start = start - 99,
                    end   = end   + 100) %>%
      (.convertDF_to_GRanges)

  }

  # Convert peak list to GRangesList
  cell_catalog  <- GRangesList(cell_catalog)

  ###################################################
  # Compute overlaps

  # Create categories labels and counts
  tf_categories <- names(cell_catalog)
  tf_catCounts  <- sapply(cell_catalog,length)

  # Get overlap counts between TF query and the rest of TFs
  overlaps <- GenomicRanges::countOverlaps(cell_catalog, tf_query, type = "any")

  ###################################################
  # Compute random shuffles for random expectation

  # Init shuffle vector
  shuffleCatCount <- vector()
  shuffleCatCount[tf_categories] <- 0

  message("Starting to compute random shuffles. This will take some minutes...\n")
  set.seed(seed = shuffle_seed)
  shuffles <- replicate(shufflesNumber, shuffle(tf_query, chromSizes, universe, included, byChrom))

  # Compute theorical means from the shuffles overlaps.
  shuffleCatCount <- sapply(shuffles, GenomicRanges::countOverlaps,
                            query   = cell_catalog,
                            type    = "any")

  shuffleCatCount <- rowMeans(shuffleCatCount)
  countsList      <- list(overlaps, shuffleCatCount)

  # Compute all enrichment scores using the (1) theoretical means,
  #                                         (2) the observed counts and
  #                                         (3) the total categories.
  enrichment_df <- extractEnrichment(tf_categories,  tail,
                                     countsList[[1]],countsList[[2]],
                                     tf_catCounts,   pAdjust)

  # Sort the enrichment df by the input threshold filter
  #if( filterBy != "fraction"){
    if( filterBy == "p.value" || filterBy == "e.value" || filterBy == "q.value" ){
      enrichment_df <- enrichment_df %>%
        arrange( !!rlang::sym(filterBy) )
    } else {
      enrichment_df <- enrichment_df %>%
        arrange( desc( !!rlang::sym(filterBy) ) )
    }
    # Get the top co-binding partners
    top_tf            <- enrichment_df  %>% head(NumberofTop)
    top_tf_category   <- top_tf %>% dplyr::pull(category)
    top_cell_TFBS     <- top_tf_category %>% lapply(function(tf){
      cell_TFBS %>% dplyr::filter( TF == tf ) %>% tail(1)
    })
    top_cell_TFBS <- do.call(rbind, top_cell_TFBS)
  #}

  # Objects to return
  enriched_cofactors <- list(
                          #TF_cell_tissue_name = TF_cell_tissue_name,
                          cell_TFBS = cell_TFBS,
                          top_cell_TFBS = top_cell_TFBS,
                          top_tf = top_tf,
                          enrichment_df = enrichment_df
                        )
}

# Create cofactor report and store it as PDF
FPWMcofactorReport_ui <- function(intersectPeakMatrix,
                           top_tf,
                           top_num             = 10,
                           filterBy            = filterBy,
                           threshold           = 0,
        				           height              = height,
	                         width               = width,
        				           pdfName             = pdfName)
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
    # check that top_tf is a valid value
    if( missing(top_tf) ){
      stop("Please provide a 'data.frame' or 'NULL' value to the 'top_tf' argument.")
    }

    # start reporting
    message("Start miniCofactorReport ...")
    message("... The maximum number of cofactors to be reported is ", top_num)
    message("... The minimum selected threshold is '",filterBy,"' >= '",threshold,"'")

    intersectPeakMatrix_i <- intersectPeakMatrix
    id_i                  <- rownames(intersectPeakMatrix_i)                 # mainTF id
    is_from_TFregulomeR   <- intersectPeakMatrix_i[1,1][[1]]@isxTFregulomeID # check if the mm id is correct

		message("... ... Start reporting peak id '",id_i,"' ...")
		# Build intersection fraction matrix as df
  	suppressMessages(intersectPeakMatrix_res_i <- intersectPeakMatrixResult(intersectPeakMatrix_i,
                                                                                return_intersection_matrix = TRUE,
                                                                                angle_of_matrix = "x"))
    intersect_matrix_i           <- intersectPeakMatrix_res_i$intersection_matrix
    rownames(intersect_matrix_i) <- "fraction"
    names_intersect_matrix_i  <- intersectPeakMatrix_res_i$intersection_matrix %>%
                                 names %>%
                                 strsplit("_") %>%
                                 sapply(function(name){name %>% tail(1)}) # NOTE: If the ID nomenclature change, fix it
    intersect_matrix_t_i      <- intersect_matrix_i %>%
                                 t %>%
                                 as.data.frame %>% dplyr::mutate( category = names_intersect_matrix_i,
                                                                  id       = intersect_matrix_i %>% names)
    if( filterBy != "fraction") {
      # Merge with top_tf
      intersect_matrix_t_i <- intersect_matrix_t_i %>%
        dplyr::left_join( top_tf , by = "category")
    } else{
      # Convert theshold
      threshold <- threshold*100
      # Select top
      intersect_matrix_t_i <- intersect_matrix_t_i %>%
        dplyr::arrange(desc(fraction)) %>%
        head(top_num)
    }
    # Filter by threshold
    intersect_matrix_filter_i <- intersect_matrix_t_i %>%
      dplyr::filter( !!rlang::sym(filterBy) >= threshold )

    # Drop all values except the filter value and rename rows
    rownames(intersect_matrix_filter_i) <- intersect_matrix_filter_i %>% dplyr::pull("id")
    intersect_matrix_filter_i           <- intersect_matrix_filter_i %>% dplyr::select( !!rlang::sym(filterBy) )

    if (nrow(intersect_matrix_filter_i) < 1)
    {
        stop("No overlap between your TF of interest and other TFs was found with the specified paramaters.")
    }

    intersect_matrix_order_i   <- intersect_matrix_filter_i[order(intersect_matrix_filter_i[,1]),,drop = FALSE]
    colnames(intersect_matrix_order_i) <- id_i
    # cobinding barplot
    intersect_matrix_heatmap_i <- as.data.frame(matrix(nrow=nrow(intersect_matrix_order_i), ncol = 5))

    colnames(intersect_matrix_heatmap_i) <- c("x","new_x","y","new_y","value")
    intersect_matrix_heatmap_i$y         <- rownames(intersect_matrix_order_i)
    intersect_matrix_heatmap_i$new_y     <- paste(unlist(lapply(rownames(intersect_matrix_order_i),
                                                          function(x) tail(unlist(strsplit(x,split = "_")),1))),
                                                   rev(rownames(intersect_matrix_heatmap_i)), sep = "-")

    intersect_matrix_heatmap_i$x         <- colnames(intersect_matrix_order_i)
    intersect_matrix_heatmap_i$new_x     <- unlist(lapply(colnames(intersect_matrix_order_i),
                                                          function(x) tail(unlist(strsplit(x,split = "_")),1)))

    intersect_matrix_heatmap_i$value     <- intersect_matrix_order_i[,1]
    intersect_matrix_heatmap_i$new_y     <- factor(intersect_matrix_heatmap_i$new_y,
                                                   levels = as.character(intersect_matrix_heatmap_i$new_y))

    cobinding_ylabel              <- as.character(intersect_matrix_heatmap_i$new_y[rev(seq(1,nrow(intersect_matrix_heatmap_i),1))])
    cobinding_ylabel_new          <- paste0(as.character(intersect_matrix_heatmap_i$new_x),
                                            "\n+\n",cobinding_ylabel)

    colors_cobinding              <- colorRampPalette(c("white","#D46A6A", "#801515", "#550000"))(11)
    cobinding_ylabel_new          <- rev(gsub("\\-\\d+$","",cobinding_ylabel_new,perl=TRUE))

    topIX                         <- length(cobinding_ylabel_new)
    cobinding_ylabel_new[ topIX ] <- paste0("All\n",as.character(intersect_matrix_heatmap_i$new_x)[topIX]) # adding the label all

    # Round numbers
    intersect_matrix_heatmap_i$value <- intersect_matrix_heatmap_i$value %>%
      round(digits=3)
    # Get max color limit
    max_color_value <- intersect_matrix_heatmap_i %>%
      dplyr::pull(value) %>% max()
    # Plot barplot
    p1 <- ggplot(intersect_matrix_heatmap_i, aes(x=new_y, y=value, fill = value)) +
      geom_bar(stat="identity",colour="black") +
      scale_fill_gradientn(colours = colors_cobinding,
                           breaks  = c(seq(0, max_color_value, length =  (top_num + 1) )),
                           limits  = c(0, max_color_value)) +
      scale_y_reverse(position = "right") +
      coord_flip() +
      theme(panel.background = element_blank(),
            plot.background = element_blank(),
            axis.title=element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
            axis.text.x = element_text(size=10),axis.text.y = element_text(size=12),
            legend.position = "none") +
      scale_x_discrete(labels=cobinding_ylabel_new)
    # Add '%' to fraction text
    if ( filterBy == "fraction"  ) {
    	p1 <- p1 +
	  geom_text(aes(label=paste0(value,"%")), vjust=-.2, angle =90)
    } else {
        p1 <- p1 + geom_text(aes(label=value), vjust=-.2, angle =90)
    }

    # motifs
    motif_plot_list_p <- list()
    count <- 1
    for (j in order(-seq(1,nrow(intersect_matrix_heatmap_i),1)))
    {
       x_j <- intersect_matrix_heatmap_i$x[j]
       y_j <- intersect_matrix_heatmap_i$y[j]
       # Retrieve the PSSM from the intersectPeakMatrix_i matrix
       motif_j <- t(intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_x@MMmotif@motif_matrix)
       top <- 10
       bottom <- 0
       if (j==1)
       {
              bottom <- 10
              top <- 0
       }

       ##########################################################################################
       #  Edited from plotLogo.R in TFregulomeR-dev 2019-09-18 :

       if (is_from_TFregulomeR)  { MM_object = intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_x }
       if (!is_from_TFregulomeR) { MM_object = intersectPeakMatrix_i[x_j,y_j][[1]]@MethMotif_y }

       y_max        <- 2
       logo_method  <- "bits"

       # get beta score matrix, motif matrix and nb sites
       MMBetaScore  <- MM_object@MMBetaScore
       MMmotif      <- MM_object@MMmotif
       motif_matrix <- t(MMmotif@motif_matrix)
       nb_sites     <- MMmotif@nsites %>% as.character %>% paste0(" sites")
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
      #labs(subtitle = nb_sites) +
      theme(axis.title.y=element_blank(), plot.subtitle = element_text(hjust = 0.5),
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
      p2 <- arrangeGrob(grobs = motif_plot_list_p,
                        nrow  = nrow(intersect_matrix_heatmap_i))

      text1 <- textGrob(paste0("Co-binding (",filterBy,")"))
      text2 <- textGrob("Motif")


      pdf_i_name <- paste0(id_i,"_cofactor_minireport.pdf")

      if( !is.null(pdfName) ){
        	pdf_i_name <- pdfName
        	pdftest <- grep(".pdf$",pdf_i_name,perl=TRUE)
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
          message("... ... ... miniCofactor report for id '", id_i,"' has been saved as ", pdf_i_name)

        }else{
          pdf(pdf_i_name,height=height,width=width)
          grid.arrange(text1,text2,p1, p2,
                     layout_matrix = def_layout_matrix)
          dev.off()
          message("... ... ... miniCofactor report for id '", id_i,"' has been saved as ", pdf_i_name)
      }

}

# thigs to fix: for GTRD matrix with no meth data.


