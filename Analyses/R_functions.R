


##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

#### Barplot
plotGeneTPM<-function(object, gene, logMode=FALSE){
    if (grepl("MSTRG", gene)){
    geneName<-subset(geneFeatures, mstrg_ID == gene)$gene_ID
  } else {geneName<-subset(geneFeatures, gene_ID == gene)$mstrg_id}
    annot <- subset(ref_to_gene_to_annots, gene_id == gene | gene_name == gene)$description
    p <- ggplot(subset(object, gene_id == gene), aes(sample, TPM, fill = Tissue)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) + facet_grid(~Source, scales="free_x", space = "free_x") + labs(title = paste(gene," (", geneName,")\n",annot, sep = "")) + theme(axis.text.x=element_text(angle = 45, hjust = 1))
    if (logMode)
    {
        p <- p + scale_y_log10()
    }
    if (logMode)
    {
        p <- p + ylab("log10 TPM")
    } else {
        p <- p + ylab("TPM")
    }
    return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

####### Boxplot 
plotBoxTPM<-function(object, gene){
    if (grepl("MSTRG", gene)){
        geneName<-subset(geneFeatures, mstrg_ID == gene)$gene_ID
    } else {geneName<-subset(geneFeatures, gene_ID == gene)$mstrg_id}
    l5.annot <- subset(geneFeatures, mstrg_ID == gene | gene_ID == gene)$GenBank_description
    p <- ggplot(subset(object, gene_id == gene), aes(sample, TPM, colour = Tissue)) + 
        facet_grid(.~Source, scales = "free") +
        geom_boxplot(width = 0.2, fill = "black", position = "dodge", outlier.size = NULL) + 
        geom_jitter(width = 0.3) + theme_bw() +
        labs(title = paste(gene," (", geneName,"). ", "\nAaegL5 description: ", l5.annot ,sep = "")) + 
        theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
    p <- p + ylab("TPM")
return(p)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## extract gene IDs based on GO term and enriched gene list:
extract_GO_genes = function(go_term, gene_set){
    rownames(subset(GO_info, row.names(GO_info) %in% gene_set & grepl(go_term, GO_info$V2)))
    
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## A function to calculate the tissue specificity index (based on CummerBund's S function)
calcSpecificity<-function(matrix,logMode=T,pseudocount=1,relative=FALSE){
    tpms<-matrix
    if(logMode){
        tpms<-log10(tpms+pseudocount)
    }
    tpms<-t(makeprobs(t(tpms)))
    d<-diag(ncol(tpms))
    res<-apply(d,MARGIN=1,function(q){
        JSdistFromP(tpms,q)
    })
    colnames(res)<-paste(colnames(tpms))
    
    if(relative){
        res<-res/max(res)
    }
    1-res
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## estimate number of cluster for K-means
findK<-function(object, k.range=c(2:20), logMode=T, pseudocount=1,...){
    require(cluster)
    m<-as.data.frame(object)
    m<-m[rowSums(m)>0,]
    if(logMode){
        m<-log10(m+pseudocount)
    }
    n<-JSdist(makeprobs(t(m)))
    myWidths<-c()
    for (k in k.range){
        #print(k)
        myWidths<-c(myWidths,pam(n,k,...)$silinfo$avg.width)
    }
    plot(k.range,myWidths)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Modifications of functions to compare groups of lists 
## (from http://stackoverflow.com/questions/23559371/how-to-get-the-list-of-items-in-venn-diagram-in-r)
Intersect <- function (x) {  
    # Multiple set version of intersect
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        intersect(x[[1]], x[[2]])
    } else if (length(x) > 2){
        intersect(x[[1]], Intersect(x[-1]))
    }
}
#
Union <- function (x) {  
    # Multiple set version of union
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        union(x[[1]], Union(x[-1]))
    }
}
#
Setdiff <- function (x, y) {
    # Remove the union of the y's from the common x's. 
    # x and y are lists of characters.
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Ouput the color IDs used by ggplot
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Miscellaneous operators
'%!in%' <- function(x,y)!('%in%'(x,y))

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

# TPM function
tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## function for plotting comprehensive gene expression data from the aegypti atlas

compGeneBarplot<-function(gene, logMode=FALSE){
    
    ## Now the expression plots
    gene_sub_tpm.c = subset(compRNAseq, ref_gene_id == gene)
    
    # gene_sub_tpm.c = summarySE(gene_sub_tpm, measurevar = "TPM", groupvars = c("gene_id", "ref_gene_id", "description", "gene_name", "Sample_Name", "alt_name", "dev_stage", "sex", "tissue", "strain", "sprot_Top_BLASTX_hit_description", "Source", "Source_description", "coordinates", "sprot_Top_BLASTX_hit_description", "description", "gene_ontology_pfam"))
    
    gene_sub_tpm.c$alt_name = factor(gene_sub_tpm.c$alt_name, levels = c("0-2hrs", "2-4hrs", "4-8hrs", "8-12hrs", "12-16hrs", "16-20hrs", "20-24hr", "24-28hrs", "28-32hrs", "32-36hrs", "36-40hr", "40-44hrs", "44-48hrs", "48-52hrs", "52-56hr", "56-60hr", "60-64hrs", "64-68hrs", "68-72hrs", "72-76hrs", "1st instar", "2nd instar", "3rd instar", "4th instar", "early", "middle", "late", "virgin", "0hpm", "6hpm", "24hpm", "not BF", "12hr post-BF", "24hr post-BF", "36hr post-BF", "48hr post-BF", "48hrs post-BF", "60hr post-BF", "72hr post-BF", "96hrs post-BF", "larvae (m & f)", "adult female", "female", "male", "males and females", "ovaries", "testes", "testes and AG", "mated", "Acc_glnds", "carcass"))
    
    gene_sub_tpm.c$tissue = factor(gene_sub_tpm.c$tissue, levels = c("embryos", "larvae", "pupae", "proboscises", "malpighian tubule", "SGs", "Pr","MPs", "fat body", "antennae", "rostrums", "brains", "carcass", "GD male", "forelegs", "midlegs", "hindlegs", "gonads", "ovaries", "AG", "male accessory gland", "accessory glands", "abdominal tips", "testes", "lower female reproductive tract", "sperm"))
    
    Akbari_adults = ggplot(subset(gene_sub_tpm.c, Source == "Cal Tech 2013" & dev_stage == "adult"), aes(alt_name, TPM, fill = sex)) +
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        facet_grid(.~tissue, scales = "free", space = "free_x") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        scale_fill_manual(values = wes_palette("Royal1")) +
        labs(title = paste(subset(gene_sub_tpm.c, Source == "Cal Tech 2013")$Source_description, ": adult tissues", sep = ""))
    
    Akbari_embryos_larvae_pupae = ggplot(subset(gene_sub_tpm.c, Source == "Cal Tech 2013" & tissue == "embryos"| Source == "Cal Tech 2013" & tissue == "larvae" | Source == "Cal Tech 2013" & tissue == "pupae"), aes(alt_name, TPM, fill = tissue)) +
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        facet_grid(.~tissue, scales = "free", space = "free_x") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), legend.position="none", axis.title.x = element_blank()) +
        scale_fill_manual(values = wes_palette("Royal2")) +
        labs(title = paste(gene," (",gene_sub_tpm.c$gene_name,")\n",
                           "Genomic location: ", gene_sub_tpm.c$coordinates, "\n",
                           "Description:", "\n-",
                           gene_sub_tpm.c$sprot_Top_BLASTX_hit_description, " (SwissProt) \n-",
                           str_wrap(gene_sub_tpm.c$description, width = 110)," (AaegL5) \n",
                           "Gene Ontology:", str_wrap(gene_sub_tpm.c$gene_ontology_pfam, width = 110), "\n\n",
                           subset(gene_sub_tpm.c, Source == "Cal Tech 2013")$Source_description, ": embryos, larvae, and pupae", 
                           sep = ""))
    
    Matthews = ggplot(subset(gene_sub_tpm.c, Source == "Rockefeller 2016"), aes(alt_name, TPM, fill = sex)) +
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        facet_grid(.~tissue, scales = "free", space = "free_x") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        scale_fill_manual(values = wes_palette("Royal1")) +
        labs(title = subset(gene_sub_tpm.c, Source == "Rockefeller 2016")$Source_description)
    
    Rockefeller = ggplot(subset(gene_sub_tpm.c, Source == "Rockefeller 2017"), aes(alt_name, TPM, fill = sex)) +
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        facet_grid(.~tissue, scales = "free", space = "free_x") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        scale_fill_manual(values = wes_palette("Royal1")) +
        labs(title = subset(gene_sub_tpm.c, Source == "Rockefeller 2017")$Source_description)
    
    Alfonso_Parra = ggplot(subset(gene_sub_tpm.c, Source == "Cornell 2016"), aes(alt_name, TPM)) +
        geom_bar(position=position_dodge(), stat="identity", fill = "#0298c5") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        #         facet_grid(.~tissue, scales = "free") +
        theme_bw() +
        #         scale_fill_manual(values = wes_palette("Royal1")) +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        labs(title = subset(gene_sub_tpm.c, Source == "Cornell 2016")$Source_description)
    
    Degner = ggplot(subset(gene_sub_tpm.c, Source == "Cornell 2018"), aes(alt_name, TPM)) +
        geom_bar(position=position_dodge(), stat="identity", fill = "#008831") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        facet_grid(.~tissue, scales = "free", space = "free_x") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        labs(title = subset(gene_sub_tpm.c, Source == "Cornell 2018")$Source_description)
    
    Sutton = ggplot(subset(gene_sub_tpm.c, Source == "Oxford 2015"), aes(alt_name, TPM)) +
        geom_bar(position=position_dodge(), stat="identity", fill = "#cd0027") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        facet_grid(.~tissue, scales = "free", space = "free_x") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        labs(title = "Sperm")
    
    NMSU = ggplot(subset(gene_sub_tpm.c, Source == "NMSU" | Source == "NIAID"), aes(alt_name, TPM)) +
        geom_bar(position=position_dodge(), stat="identity", fill = "#d36700") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        facet_grid(.~tissue, scales = "free", space = "free_x") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        labs(title = "FB, MT, and Sal. Glnds")
    
    vtech = ggplot(subset(gene_sub_tpm.c, Source == "Virginia tech"), aes(alt_name, TPM)) +
        geom_bar(position=position_dodge(), stat="identity", fill = "#9187c6") +
        geom_errorbar(aes(ymin=TPM-se, ymax=TPM+se), width=.2, position=position_dodge(.9)) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust = 1), axis.title.x = element_blank()) +
        labs(title = "Early embryos")
    
    plots <- align_plots(Akbari_embryos_larvae_pupae, Akbari_adults, Matthews, Rockefeller, Alfonso_Parra, align = 'v', axis = 'l')
    # first_row <- plot_grid(plots[[1]])
    second_row <- plot_grid(plots[[1]], align = 'h')
    third_row <- plot_grid(plots[[2]], NMSU, align = 'h', rel_widths = c(2.5, 1))
    fourth_row <- plot_grid(plots[[3]])
    fifth_row <- plot_grid(plots[[4]])
    sixth_row <- plot_grid(plots[[5]], Degner, Sutton, vtech, align = 'h', rel_widths = c(1, 1, 1, 1), nrow = 1)
    p <- plot_grid(second_row, third_row, fourth_row, fifth_row, sixth_row, ncol = 1, rel_heights = c(3.75, 2.5, 2.5, 2, 2))
    return(p)
}

#### HEATMAP.3

## pulled from here, and then tweaked slightly: http://www.biostars.org/p/18211/
 
# CODE
 
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.1,
                      #cexRow = 0.2 + 1/log10(max(nr,2)),
                      #cexCol = 0.2 + 1/log10(max(nc,2)),
        cexRow = 0.2,
        cexCol = 0.2,                  

        scaleRangeMin,
        scaleRangeMax,


    cex.main = 1,
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
 
    invalid <- function (x) {
      if (missing(x) || is.null(x) || length(x) == 0)
          return(TRUE)
      if (is.list(x))
          return(all(sapply(x, invalid)))
      else if (is.vector(x))
          return(all(is.na(x)))
      else return(FALSE)
    }



    x <- as.matrix(x)
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }

    retval <- list()


    scale <- if (symm && missing(scale))
        "none"
    else match.arg(scale)

    dendrogram <- match.arg(dendrogram)

    trace <- match.arg(trace)

    density.info <- match.arg(density.info)

    if (length(col) == 1 && is.character(col))
        col <- get(col, mode = "function")

    if (!missing(breaks) && (scale != "none"))
        warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")

    if (is.null(Rowv) || is.na(Rowv))
        Rowv <- FALSE

    if (is.null(Colv) || is.na(Colv))
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv))
        Colv <- FALSE

    if (length(di <- dim(x)) != 2 || !is.numeric(x))
        stop("`x' must be a numeric matrix")

    nr <- di[1]
    nc <- di[2]

    if (nr <= 1 || nc <= 1)
        stop("`x' must have at least 2 rows and 2 columns")
    #print(paste("nr:", nr, "nc:", nc, "cexCol:", cexCol, "cexRow:", cexRow))
    #stop("debug")



    if (!is.numeric(margins) || length(margins) != 2)
        stop("`margins' must be a numeric vector of length 2")

    if (missing(cellnote))
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))

    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% c("both", "row"))) {
            if (is.logical(Colv) && (Colv))
                dendrogram <- "column"
            else dedrogram <- "none"

            warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting row dendogram.")
        }
    }

    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv))
                dendrogram <- "row"
            else dendrogram <- "none"

            warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                dendrogram, "'. Omitting column dendogram.")
        }
    }
 
   if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
 
   if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc)
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm)
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }

    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()

    x <- x[rowInd, colInd]  # rearrange matrix according to dendrograms
    x.unscaled <- x

    cellnote <- cellnote[rowInd, colInd]  # also rearrange the cellnotes

    # get labels 
    if (is.null(labRow))
        labRow <- if (is.null(rownames(x)))
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol))
        labCol <- if (is.null(colnames(x)))
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]


    ## do scaling of matrix according to Z-scores
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }

    # number of breaks
    if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
        if (missing(col) || is.function(col))
            breaks <- 16
        else breaks <- length(col) + 1
    }

    # set breakpoints
    if (length(breaks) == 1) {
        if (missing(scaleRangeMin))
            scaleRangeMin = min(x, na.rm=na.rm)

        if (missing(scaleRangeMax))
            scaleRangeMax = max(x, na.rm=na.rm)


        if (!symbreaks) {
            breaks <- seq(scaleRangeMin, scaleRangeMax, length=breaks);
        } else {
            #extreme <- max(abs(x), na.rm = TRUE)
            extreme = max(abs(c(scaleRangeMin,scaleRangeMax)), na.rm=na.rm)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }

    nbr <- length(breaks)
    ncol <- length(breaks) - 1

    if (class(col) == "function")
        col <- col(ncol)

    min.breaks <- min(breaks)
    max.breaks <- max(breaks)

    # adjust for out-of-range given break settings
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks

    # layout height
    if (missing(lhei) || is.null(lhei))
        lhei <- c(keysize, 4)

    # layout width
    if (missing(lwid) || is.null(lwid))
        lwid <- c(keysize, 4)

    # define the layout
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
 
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || ncol(ColSideColors) != nc)
                stop("'ColSideColors' must be a matrix of ncol(x) ", nc, " columns")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
            #lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
            side_height = min(side.height.fraction*nrow(ColSideColors), 1);
            lhei=c(lhei[1], side_height, lhei[2])
        }
 
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || nrow(RowSideColors) != nr)
                stop("'RowSideColors' must be a matrix of nrow(x) ", nr, " rows.  It currently has ", nrow(RowSideColors), " rows.")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
            #lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
            side_width = min(side.height.fraction*ncol(RowSideColors), 1);
            lwid <- c(lwid[1], side_width, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
 
    if (length(lhei) != nrow(lmat))
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat))
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    

    op <- par(no.readonly = TRUE)
    on.exit(par(op))
 
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
 
    ###########################################
    ## Draw the colorbars for the annotations:
    ########################################### 

    if (!missing(RowSideColors)) {
        if (!is.matrix(RowSideColors)){
                par(mar = c(margins[1], 0, 0, 0.5))
                image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        } else {
            par(mar = c(margins[1], 0, 0, 0.5))
            rsc = t(RowSideColors[rowInd, , drop=F])
            rsc.colors = matrix()
            rsc.names = names(table(rsc))
            rsc.i = 1
            for (rsc.name in rsc.names) {
                rsc.colors[rsc.i] = rsc.name
                rsc[rsc == rsc.name] = rsc.i
                rsc.i = rsc.i + 1
            }
            # print(rsc)
            rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
            #print("RSC: ", rsc)
            #print(rsc.colors)    
            image(1:nrow(rsc), 1:ncol(rsc), rsc, col = as.vector(rsc.colors), axes = FALSE, xlab="", ylab="")
        
            # add labels
            if (length(colnames(RowSideColors)) > 0) {  
                #axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
                #axis(1, 0:(nrow(rsc)-1), colnames(RowSideColors), las = 2, tick = T) # ncol because transposed
                axis(1, 1:ncol(RowSideColors), labels=colnames(RowSideColors), las=2, cex.axis=0.5, tick=F, xlab="", ylab="")

            }
        }
    }
    


    if (!missing(ColSideColors)) {
 
        if (!is.matrix(ColSideColors)){
            par(mar = c(0.5, 0, 0, margins[2]))
            image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        } else {
            par(mar = c(0.5, 0, 0, margins[2]))
            csc = ColSideColors[, colInd, drop=F]
            csc.colors = matrix()
            csc.names = names(table(csc))
            csc.i = 1
            for (csc.name in csc.names) {
                csc.colors[csc.i] = csc.name
                csc[csc == csc.name] = csc.i
                csc.i = csc.i + 1
            }
            csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
            #print(csc)
            image(1:nrow(t(csc)), 1:ncol(t(csc)), t(csc), col = as.vector(csc.colors), axes = FALSE, xlab="", ylab="")

            # add labels
            if (length(rownames(ColSideColors)) > 0) {
                #axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
                axis(2, 1:(nrow(ColSideColors)), labels=rownames(ColSideColors), las = 2, tick = FALSE, cex.axis=0.5)
            }
        }
    }
 


    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr"))
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    
    # draw the central heatmap
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
    
    # store the matrix drawn
    retval$carpet <- x
    
    # store the dendrograms
    if (exists("ddr"))
        retval$rowDendrogram <- ddr
    if (exists("ddc"))
        retval$colDendrogram <- ddc
    
    # store the breaks
    retval$breaks <- breaks
    
    # store the colormap used
    retval$col <- col
    
    # specially color in the na values  
    if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", col = na.color, add = TRUE)
    }

    # X-axis column labels
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol)

    # X-axis title
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)

    # Y-axis row labeling
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
        cex.axis = cexRow)

    # Y-axis title
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
 
    if (!missing(add.expr))
        eval(substitute(add.expr))
    if (!missing(colsep))
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep))
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    

    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    
    # column trace
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }

    # row trace
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }

    # add cell labels
    if (!missing(cellnote))
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), col = notecol, cex = notecex)

    ###########################
    ## Plot the row dendrogram
    ###########################

    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()

    #############################
    ## Plot the column dendrogram
    #############################

    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()

    if (!is.null(main))
        title(main, cex.main=cex.main) #cex.main = 1.5 * op[["cex.main"]])


    ############################
    ## Add the Color Chart
    ############################

    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(c(x,breaks), na.rm = TRUE)
            max.raw <- max(c(x,breaks), na.rm = TRUE)
        }
 
        message('for plotting:: min.raw: ', min.raw, ' max.raw: ', max.raw);
        
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row")
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column")
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, KeyValueName, line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()

    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], high = retval$breaks[-1], color = retval$col)

    invisible(retval)
}



# EXAMPLE USAGE
 
# example of colsidecolors rowsidecolors (single column, single row)
#mat <- matrix(1:100, byrow=T, nrow=10)
#column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
#column_annotation <- as.matrix(column_annotation)
#colnames(column_annotation) <- c("Variable X")
 
#row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
#row_annotation <- as.matrix(t(row_annotation))
#rownames(row_annotation) <- c("Variable Y")
 
#heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
 
# multiple column and row
#mat <- matrix(1:100, byrow=T, nrow=10)
#column_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), ncol=2)
#colnames(column_annotation) <- c("Variable X1", "Variable X2")
 
#row_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), nrow=2)
#rownames(row_annotation) <- c("Variable Y1", "Variable Y2")
 
#heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
 


#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------

# Borrowing the following from gplots  (gplots isn't compatible with R 3.0 (yet), and so bypassing it for now).

colorpanel = function (n, low, mid, high) 
{
    if (missing(mid) || missing(high)) {
        low <- col2rgb(low)
        if (missing(high)) 
            high <- col2rgb(mid)
        else high <- col2rgb(high)
        red <- seq(low[1, 1], high[1, 1], length = n)/255
        green <- seq(low[3, 1], high[3, 1], length = n)/255
        blue <- seq(low[2, 1], high[2, 1], length = n)/255
    }
    else {
        isodd <- odd(n)
        if (isodd) {
            n <- n + 1
        }
        low <- col2rgb(low)
        mid <- col2rgb(mid)
        high <- col2rgb(high)
        lower <- floor(n/2)
        upper <- n - lower
        red <- c(seq(low[1, 1], mid[1, 1], length = lower), seq(mid[1, 
            1], high[1, 1], length = upper))/255
        green <- c(seq(low[3, 1], mid[3, 1], length = lower), 
            seq(mid[3, 1], high[3, 1], length = upper))/255
        blue <- c(seq(low[2, 1], mid[2, 1], length = lower), 
            seq(mid[2, 1], high[2, 1], length = upper))/255
        if (isodd) {
            red <- red[-(lower + 1)]
            green <- green[-(lower + 1)]
            blue <- blue[-(lower + 1)]
        }
    }
    rgb(red, blue, green)
}


greenred = function (n)  {
    colorpanel(n, "green", "black", "red")
}

odd = function (x) {
    x%%2 == 1
}

even = function (x) {
    x%%2 == 0
}


plot_counts_matrix_log2_dist = function(matrix_file) {

    
    data = read.table(file=matrix_file, com='', row.names=1, header=T)

    conditions = colnames(data)
    colors = rainbow(length(conditions))


    plot(density(log2(data[,1])), col=colors[1], main=matrix_file, xlab='log2(frag_counts)', ylab='density')

    for (i in 2:length(data[1,])) {

        points(density(log2(data[,i])), type='l', col=colors[i])

    }

    legend('topright', conditions, col=colors, pch=15)

}


matrix_to_color_assignments = function(matrix_m, col=NULL, by=c("matrix", "row", "col")) {

    if (! is.matrix(matrix_m))
        stop("Error, matrix_to_color_assignments() requires a matrix as parameter.")
    num_colors = 0
    
    if (is.null(col)) {
        num_colors = min(nrow(matrix_m), ncol(matrix_m))
        col = rainbow(num_colors)
    }
    else {
        num_colors = length(col)
    }
    
    by = match.arg(by)
    
    if (by == "matrix") {

        min_val = min(matrix_m)
        matrix_m = matrix_m - min_val
        max_val = max(matrix_m)
        matrix_m = matrix_m / max_val * num_colors
        #print(matrix_m)
        matrix_m = apply(matrix_m, 1:2, function(x) ifelse (x<1, as.character(col[1]), as.character(col[x])));
        
        matrix_m = matrix(as.character(matrix_m), nrow=dim(matrix_m)[1])
    }
    else {

        row_or_col_only_color_selector_func = function(x) { 
                a = min(x); 
                b = max(x); 
                c = (x-a)/(b-a) * num_colors;
                c = round(c);
                c = ifelse (c<1, 1, c); 
                #print(paste(c("color selection: (a)", a, " (b)", b, " (c)", paste(c, sep=',')))); 
                colors = as.character(col[c]);
                return(colors);
        }
    
        if (by == "row") {
            matrix_m = t(apply(matrix_m, 1, row_or_col_only_color_selector_func));
        }
        else {
            # by column
            matrix_m = apply(matrix_m, 2, row_or_col_only_color_selector_func);
        }
    }
    
    #print(matrix_m)
    return(matrix_m)
}

sample_matrix_to_color_assignments = function(sampleAnnotationsMatrix, colors) {

    if (missing(colors))
        colors = rainbow(nrow(sampleAnnotationsMatrix))

    nsamples = nrow(sampleAnnotationsMatrix);

    if (length(colors) < nrow(sampleAnnotationsMatrix))
        stop("Error, only ", length(colors), " colors specified, but have ", nsamples, " samples");

    for (i in 1:nrow(sampleAnnotationsMatrix)) {
        c = colors[i]
        sampleAnnotationsMatrix[i,] = sapply(sampleAnnotationsMatrix[i,], function(x) ifelse( x, as.character(c), 'white'))
    }

    return(sampleAnnotationsMatrix);

}
