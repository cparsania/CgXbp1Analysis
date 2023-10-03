
##########################################################################################################################################
#-------------------------------------------------------------plot_sd_distr_plot ------------------------------------------------------###
##########################################################################################################################################


## funct. plot sd distr
plot_sd_distr_plot <- function(x , col_names = NULL) {
        ## x is a dataFrame or tibble having rows are genes and columns are variables.
        #x <- polII_FPKM_tbl
        # col_names <- c("polII_0_hr","polII_2_hr","polII_4_hr","polII_6_hr","polII_8_hr")
        
        stopifnot(!is.character(col_names))
        x <- as_tibble(x) %>% select_if(isIntOrDoub)
        colsForSd <- colnames(x)
        
        if (!is.null(col_names)) {
                colsForSd = col_names
                x <- as_tibble(x) %>% select(col_names)
        }
        
        ## only double variable will be selected.
        stopifnot(ncol(x) > 0 && nrow(x) > 0)
        
        sdCol <-
                x %>% mutate(stdev = rowSds(as.matrix(.[colsForSd]))) %>% arrange(desc(stdev)) %>% select(stdev)
        ggFrmtStDev <- melt(sdCol)
        gp <-
                ggplot(ggFrmtStDev) + geom_histogram(
                        aes(x = value) ,
                        fill = "cyan3" ,
                        col = "black",
                        binwidth = 2.5
                ) + xlab("stdev across all timepoints")
        gp <- myGGPlotAnnotations(gp)
        return(gp)
}


isIntOrDoub <- function(x) {
        return(is.double(x) | is.integer(x))
}

##########################################################################################################################################
#-------------------------------------------------------------get_tidy_fpkm_mat ------------------------------------------------------###
##########################################################################################################################################


#' Title : get_tidy_fpkm_mat
#' @description : given the FPKM mat of several different conditions where rownames are genes and columns are condition wise FPKM values, convert it in to tidy format
#' @param x : A numeric matrix or data.frame have rownames indicating genenames and columns are fpkm of individual samples
#' @param organism : A character vector of length 1 or equal to number of columns in the x indicating organism name for each column. (cg, an etc.) , default "NA"
#' @param strain  : A character vector of length 1 or equal to number of columns in the x indicating strain information for each column. (e.g wt, deletion etc.),  default "NA"
#' @param condition : A character vector of length 1 or equal to number of columns in the x indicating experimental condition for each column. (e.g rpmi, thp1 etc.) , default "NA"
#' @param timepoint  : A character vector of length 1 or equal to number of columns in the x indicating timepoint for each column. (e.g. 0hr, 2hr, etc.), default "NA"
#' @param moleculeType : A character vector of length 1 or equal to number of columns in the x indicating moleculeType for each column. (e.g. rnaseq , polII, h3k4me3 etc.), default "NA"
#' @param replicate : A character vector of length 1 or equal to number of columns in the x indicating replicate informations for each column. (e.g rep1,rep2, etc.), default "NA"
#' @param format : Either of the 'wide' or 'long' indicating return format of the data. If 'wide', sample wise fpkm will be displyed in each column and gene names will appear only once, while long will combine fpkm of all samples under one column and gene names will be repeated , default 'long'.
#' @return : A tibble having
#' @export
#'
#' @examples
get_tidy_fpkm_mat <-
        function(x ,
                 organism = "NA",
                 strain = "NA",
                 condition = "NA",
                 timepoint = "NA",
                 moleculeType  = "NA" ,
                 replicate = "NA" ,
                 format = "long") {
                # test variable
                # x = cg_xbp1_polII_FPKM_mat
                # organism = "NA" #"cg"
                # strain = "NA" #"xbp1"
                # condition = "NA" #"thp1"
                # timepoint =  paste(rep(seq(0,8,by = 2) , times = 2) ,"hr" , sep="")
                # replicate = rep(c("rep1" , "rep2") , each = 5)
                
                ## x must be data.frame or matrix have all columns numeric(int or double). Given rownames will be used as each row identifier
                stopifnot(is.numeric(as.matrix(x)))
                
                ## convert to tibble taking rownames as genenames
                x <-
                        as.data.frame(x) %>% rownames_to_column("geneName") %>% as_tibble()
                
                ## check function arguments. length of organism, strain, condition, timepoint and replicate variables must be either 1 or equal to number of columns - 1. First column is geneNames will not be considered
                stopifnot(lengths(
                        list(
                                organism,
                                strain,
                                condition,
                                timepoint,
                                replicate,
                                moleculeType
                        )
                ) == 1 |
                        lengths(
                                list(
                                        organism,
                                        strain,
                                        condition,
                                        timepoint,
                                        replicate,
                                        moleculeType
                                )
                        ) == ncol(x) - 1)
                
                ## check 'format' argument. It cannot be other than wide or long
                stopifnot(format %in% c("wide" , "long"))
                
                ## create new column names, generated from function arguments.
                new_col_names <-
                        paste(organism ,
                              strain,
                              condition ,
                              timepoint ,
                              moleculeType ,
                              replicate ,
                              sep = "_")
                
                ## number of unique column names must be equal to number of columns in x
                if (length(unique(new_col_names)) != ncol(x) - 1) {
                        stop(
                                "Cannot generate unique column names from the combinations of organism, strain, condition, timepoint, moleculeType and replicate arguments equal to number of columns in x"
                        )
                }
                
                ## Print the message for user
                message(
                        "\nNew column names are\n------------------------\n",
                        paste("-", new_col_names, collapse = "\n"),
                        "\n-----------------------\n",
                        sep = ""
                )
                message(
                        "Make sure that x have same column order as displayed above. If not, change the arguments accordingly and rerun the function\n"
                )
                
                Sys.sleep(time = 1)
                
                x <-
                        x %>% `colnames<-` (c("geneName", new_col_names))
                
                if (format == "wide") {
                        return(as_tibble(x))
                }
                
                ## make tidy format
                x <-
                        x %>% gather("Cond" , "FPKM" , new_col_names) %>%
                        separate(
                                col = "Cond",
                                sep = "_" ,
                                into = c(
                                        "organism",
                                        "strain",
                                        "condition",
                                        "timepoint",
                                        "moleculeType",
                                        "replicate"
                                )
                        )
                
                return(x)
                
        }




###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------get_one_to_one_corr ---------------------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#' Title get_one_to_one_corr
#'
#' @param dat : tibble or data.frame having columns given in set1 and set2
#' @param set1 : A character vector, must  have identical length to set2.
#' @param set2 : A character vector, must  have identical length to set1.
#' @description : Given the same length of two vectors refering to column names of given dataset (dat), function will calculate the correlation between respective elements from set1 and set2. Note : Rows containing na will be removed from the given data (dat)
#' (E.g. Correlation between first element , between second element and so on for two vectors.)
#' @return : a data.frame
#' @export
#'
#' @examples
get_one_to_one_corr <- function(dat, set1 , set2) {
        # dat = polII_FPKM_tbl
        # set1 = rep_set1
        # set2 = rep_set2
        
        if (!is.tibble(dat) & !is.data.frame(dat)) {
                stop("is.tibble(dat) or is.data.frame(data) one of these must be TRUE")
        }
        
        ## remove na rows if any
        dat <- na.omit(dat)
        
        if (length(set1) != length(set2)) {
                stop("length(set1) == length(set2) is not TRUE")
        }
        corrs <- sapply(seq_along(set1), function(i) {
                #i = 1
                set1_data <- dat[[set1[i]]]
                set2_data <- dat[[set2[i]]]
                return(cor(set1_data, set2_data))
        })
        
        return(cbind.data.frame(set1, set2 , cor = round(corrs, 4)))
        
}


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------counts_by_sd_cutoffs --------------------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


## check number of genes at each interval of SD cutoff

#' Title : counts_by_sd_cutoffs
#' @description : given an tibble or data frame (x) along with column names (y) to calculate SD by row, function returns plot of counts at various sd cutoffs (sds_cutoffs)
#' @param x  : An object of class tibble or data.frame.
#' @param y : A character vector for indicating column names from x to calculate standard deviation.
#' @param sds_cutoffs : A numeric/double vector indicating values for sd cutoffs
#' @param plot_title : A string for plot title
#' @return : An object of ggplot
#' @export
#'
#' @examples
counts_by_sd_cutoffs <-
        function(x ,
                 y ,
                 sds_cutoffs = c(1:10),
                 plot_title = "SD cutoffs vs Counts") {
                ## test variables
                # x = polII_FPKM_tbl
                # y = syms(c("xbp1_thp1_set1_0h","xbp1_thp1_set1_2h","xbp1_thp1_set1_4h","xbp1_thp1_set1_6h","xbp1_thp1_set1_8h"))
                # plot_title = "SD cutoffs vs Counts"
                # sds_cutoffs <- seq(from = 1, to = 10, by = 0.5)
                
                
                ## conver var y to symbols
                y = syms(y)
                
                ## logic goes here
                counts <- sapply(sds_cutoffs , function(elem) {
                        #elem <- 1.5
                        x %>%
                                dplyr::select_if(is.double) %>%
                                na.omit() %>%
                                dplyr::rowwise() %>%
                                dplyr::mutate(stdev = sd(c(!!!y))) %>%
                                dplyr::filter(stdev > elem) %>%
                                dplyr::ungroup() %>%
                                dplyr::summarise(count = dplyr::n())
                }, USE.NAMES = T)
                
                sds_cutoff_vs_counts <-
                        cbind.data.frame(sds_cutoffs , counts = unlist(counts, use.names = F))
                gp <- sds_cutoff_vs_counts %>% ggplot() +
                        geom_bar(aes(x =  as.factor(sds_cutoffs) , y = counts) ,
                                 stat = "identity" ,
                                 fill = "cyan3") +
                        xlab("SD Cutoffs") + ylab("Nuber of genes") + ggtitle(plot_title)
                
                #gp <- myGGPlotAnnotations(gp)
                
                return(gp)
                
        }



###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------get_column_averaged_lineplot -----------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

### get_column_averaged_lineplot
#' Title : get_column_averaged_lineplot
#' @description : given a data.frame or tibble, function will calculate column wise means with in the group. Groups must be defined by one of the column of x indicated by group_by_column_name. Returns the average line plot across the columns.
#'
#' @param x : an object of class data.frame or tibble.
#' @param group_by_column_name : column names indicating groups
#' @param columns_to_plot : A character vector of column names for which average line plot required.
#'
#' @return : list of ggplot. Each plot denotes individual group.
#' @export
#'
#' @examples
get_column_averaged_lineplot <-
        function(x , columns_to_plot , group_by_column_name) {
                ## test variables
                # x <-  for_averaged_line_plot
                # group_by_column_name <- "heatMapClust"
                # columns_to_plot  <- c(paste(seq(0,8,2),"h.x",sep = ""))
                
                ## validate arguments
                stopifnot(is.tibble(x) | is.data.frame(x))
                
                ## convert column name to symbol
                columns_to_plot <- rlang::enquos(columns_to_plot)
                group_by_column_name <-
                        rlang::enquo(group_by_column_name)
                
                ## subset required columns
                x <-
                        x %>% 
                        tibble::as_tibble() %>% 
                        dplyr::select(!!!columns_to_plot ,!!!group_by_column_name)
                
                molten_data <-
                        x  %>% dplyr::group_by(!!!group_by_column_name) %>%
                        dplyr::summarise_if(is.numeric, funs(mean)) %>% ## take all numeric columns and mean it
                        melt()
                
                ## create plot list
                spltd <-
                        split(molten_data, molten_data[[group_by_column_name]])
                plotList <- list()
                for (i in names(spltd)) {
                        #i <- names(spltd)[1]
                        pp <-
                                ggplot(
                                        spltd[[i]] ,
                                        aes_string(
                                                x = "variable" ,
                                                y = "value" ,
                                                group = group_by_column_name
                                        )
                                ) + geom_line(linetype = 1, size = 2) + geom_point(size = 6, color = "black")
                        pp <- pp + ggtitle(i)
                        plotList[[i]] <- pp
                }
                
                return(plotList)
        }


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------get_enrichment_heatmap_list -------------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#' Title get_enrichment_heatmap_list
#' @description Given the list of normalizedMatrix class objects, function will return list of EnrichmentHeamtmap objects. Return objects can be directly passed to Heatmap::draw method to plot several heatmaps.
#'
#' @param x is a list where each element is an object belonging to class normalizedMatrix. Objects of normalizedMatrix class can be derived from function EnrichedHeatmap::normalizeToMatrix().
#' @param names a character vector of length x, where each element represent the heatmap legend name for corresponding element of x.
#' @param titles a character vector of length x, where each element represent the heatmap title for corresponding element of x.
#' @param ... pass to EnrichedHeatmap
#'
#' @return an object of class EnrichedHeatmapList, which can directly pass to Heatmap::draw method
#' @export
#'
#' @examples
get_enrichment_heatmap_list <- function(x, names, titles, ...) {
        # x <- xx$norm_matrix
        
        ll <- length(x)
        
        ## first heatmap
        ehml <-
                EnrichedHeatmap(
                        mat = x[[1]],
                        name = names[[1]],
                        column_title = titles[[1]],
                        show_heatmap_legend = T,
                        ...
                )
        
        ## several other heatmaps if length of x > 1.
        if (ll > 1) {
                for (i in 2:ll) {
                        ehml <- ehml +
                                EnrichedHeatmap(
                                        mat = x[[i]],
                                        name = ifelse(length(names) >= i, names[i], "NA"),
                                        column_title = ifelse(length(titles) >= i, titles[i], "NA"),
                                        show_heatmap_legend = ifelse(length(names) >= i, TRUE, FALSE),
                                        ...
                                ) ## legend will be shown only if the name is given for a heatmap.
                }
        }
        
        return(ehml)
}


###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------get_tf_summit_heatmap -------------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#' Title get_tf_summit_heatmap
#' @description  given the set of genes, peak annotations and signal file, the function will return the output of EnrichedHeatmap::EnrichedHeatmap().
#' @param query_genes : is a cheracter vector of gene names for which heatmap will be plotted. If the gene name is not found in the column gene of peak_anno_file, it will be eliminated from the plot. 
#' @param peak_anno_file : full path of the file, which must have 4 tab deliminated columns, named :  chr, summit_location, peak_name, gene. query_genes will be searched in the column gene. If found, they will be displayed in final plot. 
#' @param signal_file : full path to .bw file. It will be used to plot the chip signal.
#' @param flank_bp : an integer indicating number of basepairs extension on either side of summit. 
#' @param legend_name :  a string denoting name of the legend. 
#' @param hm_axis_name : cheracter vector of length 3 denoting names of start summit and end position in the heatmap. 
#' @param hm_axis_name_font_size : an object of gpar() , default gpar(fontsize = 12, fontface = "bold")
#' @param axis_name_angle : an integer refering angle of hm_axis_name
#' @param ... : will be passed to EnrichedHeatmap::EnrichedHeatmap()
#'
#' @return : an object of class EnrichedHeatmap::EnrichedHeatmap()
#' @export
#'
#' @examples
get_tf_summit_heatmap <- function(query_genes , 
                                  peak_anno_file ,
                                  signal_file ,
                                  flank_bp = 1500, 
                                  legend_name = "binding signal", 
                                  hm_axis_name = NULL, 
                                  hm_axis_name_font_size  =gpar(fontsize = 12, fontface = "bold"), 
                                  axis_name_angle = 90 , ...){
        
        ###### test param 
        
        # query_genes <- c("CAGL0A00737g", "CAGL0A00737g", "CAGL0A01782g", "CAGL0A01826g", "CAGL0A02211g", "CAGL0A02233g", "CAGL0A02321g", "CAGL0A03212g", "CAGL0B01012g", "CAGL0C05489g", "CAGL0D04708g", "CAGL0E01353g", "CAGL0E03355g", "CAGL0E03674g", "CAGL0F00583g", "CAGL0F01419g", "CAGL0F02717g", "CAGL0F04499g", "CAGL0G00242g", "CAGL0G02453g", "CAGL0G03487g", "CAGL0G08624g", "CAGL0H08393g", "CAGL0I02178g", "CAGL0I04862g", "CAGL0I09086g", "CAGL0I09724g", "CAGL0J01661g", "CAGL0J02948g", "CAGL0J08162g", "CAGL0J08184g", "CAGL0K05753g", "CAGL0K12210g", "CAGL0L03696g", "CAGL0L04884g", "CAGL0L07766g", "CAGL0M03465g", "CAGL0M09020g", "CAGL0M10417g")
        # 
        # peak_anno_file <- "../../../../Cg_XBP1_Chip/XBPMyc2_Day4_peak_with_genes.txt"
        # signal_file <- "../../../../Cg_XBP1_Chip/XBPMyc2_Day4_CL1019Mix_ACGTAGCTC//XBPMyc2_Day4_CL1019Mix_ACGTAGCTC_normalized.bw"
        # flank_bp <- 1500
        # 
        # legend_name <-  "binding signal"
        # hm_axis_name <- c(paste("-",flank_bp/1000 ," kb",sep = ""), "summit", paste("+",flank_bp/1000," kb", sep = ""))
        # hm_axis_name_font_size = gpar(fontsize = 12, fontface = "bold")
        # axis_name_angle = 90 
        
        ###### funct param
        query_genes = query_genes
        peak_anno_file = peak_anno_file
        signal_file = signal_file
        flank_bp = flank_bp
        legend_name = legend_name
        hm_axis_name =  hm_axis_name
        
        if(is.null(hm_axis_name)){
                hm_axis_name = c(paste("-",flank_bp/1000 ," kb",sep = ""), "summit", paste("+",flank_bp/1000," kb", sep = "")) 
        }
        
        hm_axis_name_font_size  = hm_axis_name_font_size
        axis_name_angle = axis_name_angle
        
        ### process peak anno file 
        peak_anno <- read_delim(peak_anno_file , delim = "\t") %>% 
                mutate(gene =  strsplit(gene , ","))
        
        ## if more than one gene, convert them to list 
        peak_anno_sub <- peak_anno  %>% 
                filter(map_lgl(gene, function(.){any(. %in% query_genes)} )) %>% 
                unnest()
        
        ## prepare gr by summit 
        peak_summit_gr <- GenomicRanges::GRanges(seqnames = peak_anno_sub$chr , 
                                                 ranges = IRanges(start =  peak_anno_sub$summit_location , end = peak_anno_sub$summit_location) , 
                                                 peak_name  = peak_anno_sub$peak_name , 
                                                 genes_asso = peak_anno_sub$gene)
        
        ## import signals
        sig <- rtracklayer::import(signal_file)
        nn <- EnrichedHeatmap::normalizeToMatrix(signal = sig, target = peak_summit_gr, value_column = "score", smooth = T, extend = flank_bp )
        
        hm <- EnrichedHeatmap(mat = nn , 
                              name =  legend_name, 
                              axis_name = hm_axis_name, 
                              axis_name_gp = hm_axis_name_font_size, 
                              axis_name_rot = axis_name_angle, 
                              top_annotation = HeatmapAnnotation(lines = anno_enriched(yaxis_facing = "right" , yaxis_side = "left")) , ...) 
        
        return(hm)
}



###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------get_cg_polII_heatmap -------------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#' Title get_cg_polII_heatmap
#'
#' @param x a characeter vector of the genes
#' @param type either of the "wt" or "xbp1"
#' @param ... will be passed to ComplexHeatmap::Heatmap()
#'
#' @return an boject of class ComplexHeatmap
#' @export
#'
#' @examples
get_cg_polII_heatmap <- function(x , type  = "wt" , ...) {
        #x <- c("CAGL0I07293g", "CAGL0C00968g", "CAGL0H00110g", "CAGL0E01661g", "CAGL0C01133g", "CAGL0M03773g", "CAGL0C03575g" , "AED1")
        #type = "wt"
        
        wt_fpkm_mat_file  <-
                "data/rds/cg_polII_wt_thp1_sample_wise_fpkm_mat.rds"
        cat("\nLoading object:" , wt_fpkm_mat_file, "\n")
        load(wt_fpkm_mat_file)
        wt_mat <- sample_wise_fpkm_mat_final
        rm(sample_wise_fpkm_mat_final)
        
        xbp1_fpkm_mat_file <-
                "data/rds/cg_polII_xbp1_thp1_sample_wise_fpkm_mat.rds"
        cat("\nLoading object:" , xbp1_fpkm_mat_file, "\n")
        load(xbp1_fpkm_mat_file)
        xbp1_mat <- sample_wise_fpkm_mat_final
        rm(sample_wise_fpkm_mat_final)
        
        ## wt data heatmap format
        wt_hm_format <- wt_mat %>%
                tidyr::gather("key", "value", -geneName) %>%
                separate(
                        col = "key",
                        into = c(
                                "org",
                                "strain",
                                "cond",
                                "time_point",
                                "molecule",
                                "replicate"
                        )
                ) %>%
                dplyr::select(geneName, time_point, value) %>%
                tidyr::spread(key = "time_point", value = "value") %>% dplyr::rename("0.5h" = "0h")
        
        ## xbp1 data heatmap format
        ## note : use replicate 1
        
        xbp1_hm_format <- xbp1_mat %>%
                tidyr::gather("key", "value", -geneName) %>%
                separate(
                        col = "key",
                        into = c(
                                "org",
                                "strain",
                                "cond",
                                "time_point",
                                "molecule",
                                "replicate"
                        )
                ) %>%
                dplyr::filter(replicate == "set1") %>%
                dplyr::select(geneName, time_point, value) %>%
                tidyr::spread(key = "time_point", value = "value") %>% dplyr::rename("0.5h" = "0h")
        
        ## select data for HM
        if (type == "wt") {
                hm_data <- wt_hm_format
        } else if (type == "xbp1") {
                hm_data <- xbp1_hm_format
        } else {
                stop("parameter 'type' must be either of the \"wt\" or \"xbp1\"")
        }
        
        ## subset genes for heatmap
        hm_data <-
                tibble(geneName = x) %>% left_join(hm_data) %>% tidyr::drop_na()
        
        genes_not_found  <- x[!x %in% hm_data$geneName]
        if (length(genes_not_found) >= 1) {
                warning("following genes not found in the data. " ,
                        genes_not_found)
        }
        
        hm <-
                hm_data %>% as.data.frame() %>% 
                tibble::column_to_rownames("geneName") %>% ComplexHeatmap::Heatmap(...)
        return(hm)
        
}

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------tibble_to_rowwise_sd_zs -------------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


tibble_to_rowwise_sd_zs <-
        function(x ,
                 cols_to_use ,
                 std_dev_colname =  "sd" ,
                 zscore_colname_suffix = "_zscore") {
                # x = xx
                # cols_to_use = my_cols
                # std_dev_colname = "sd"
                # zscore_colname_suffix = "_zs"
                #
                cols <- dplyr::enquos(cols_to_use)
                std_dev_colname <- enquo(std_dev_colname)
                
                #cols <- c("SRR3384897","SRR4454155" ,"SRR4454585","SRR4456864","SRR5882561")
                #cols <- parse_exprs(cols)
                
                ## cols must be numeric
                if (!(map_lgl(x %>% dplyr::select(!!!cols) , is.numeric) %>% all)) {
                        stop(paste0(
                                as.character(cols) ,
                                "all must be of type numeric." ,
                                collapse = " "
                        )) ## other than column 1, all must be of type numeric
                }
                
                ## given column names cannot be `n` or `sd`
                if (any(cols %in% c("n", std_dev_colname))) {
                        stop("column names cannot be `n` or " , std_dev_colname)
                }
                
                xx <-
                        x %>% rownames_to_column() %>% dplyr::select(rowname,!!!cols) %>% dplyr::ungroup()
                
                ## add standard deviation
                dd <- xx  %>%
                        tidyr::nest(!!!cols) %>%  ## nest sample variable in the one column called `data`
                        dplyr::mutate(zscore = map(data , function(elem) {
                                ## iterate over each elem of data to calculate rowwise zscore
                                elem_v <- as.numeric(elem)
                                zs <-
                                        scale(elem_v)  %>% as.numeric()
                                names(zs) <-
                                        paste(!!!cols , zscore_colname_suffix , sep = "")
                                return(as.tibble(as.list(zs)))
                        })) %>%
                        dplyr::mutate(!!std_dev_colname := purrr::map_dbl(data , ~
                                                                                  (sd(.x)))) %>% ## add standard deviation column of named `std_dev_colname`
                        tidyr::unnest() %>%
                        dplyr::select(-!!std_dev_colname , dplyr::everything())
                
                ## join with orig data by roworder
                dd <-
                        x  %>% rownames_to_column() %>% dplyr::left_join(dd)  %>% dplyr::select(-rowname)
                
                return(dd)
        }


## cg_common id to gene symbol mapping

get_cg_gene_names_mapping <- function() {
        id_mapping_file <-
                "../../../1_genome/CBS138_s02-m07-r06/annotation/cagl_to_common_name_mapping"
        id_map <-
                read_delim(id_mapping_file ,
                           delim = "\t" ,
                           col_names = T)
        return(id_map)
        
}

###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------heatmap_with_column_boxplot -------------------------------------------------
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#' Title given the matrix get heatmap with boxplot on the top
#'
#' @param xx
#' @param col_title
#' @param legend_name
#'
#' @return
#' @export
#'
#' @examples
heatmap_with_column_boxplot <-
        function(xx ,
                 col_title =  "" ,
                 legend_name = "" ,
                 colors = c('#fdd49e',
                            '#fdbb84',
                            '#fc8d59',
                            '#ef6548',
                            '#d7301f',
                            '#990000')) {
                col_min = min(xx)
                col_max = max(xx)
                hm_colrs <-
                        circlize::colorRamp2(
                                breaks = seq(col_min, col_max , length.out = length(colors)) ,
                                colors = colors
                        )
                ha <-
                        HeatmapAnnotation(box = anno_boxplot(
                                xx %>% as.matrix() ,
                                axis = T,
                                axis_gp = gpar(size = 10)
                        ),
                        which = "column")
                ComplexHeatmap::Heatmap(
                        xx %>%  as.matrix() ,
                        top_annotation = ha,
                        cluster_rows = T ,
                        cluster_columns = F ,
                        top_annotation_height =  unit(5 , "cm") ,
                        column_title = col_title ,
                        col = hm_colrs ,
                        name = legend_name
                )
        }



