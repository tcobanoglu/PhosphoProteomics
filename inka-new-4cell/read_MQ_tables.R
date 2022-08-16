#!/usr/bin/env Rscript

library(data.table)
library(stringr)

# This path needs to be set correctly in order
# for the script to find all necessary resources

install.location <- "./"

# ---------------------------------------------------------------------------------------
# I want oto get rid of this one !!!
mapping.file <- paste0(install.location,"/valid_mapping.txt")

# ---------------------------------------------------------------------------------------

manning.kinase.file <- paste0(install.location,"/Kinases_HGNC2017.01.20_v2.Rdata")

complete.hgnc.path <- paste0(install.location,"/","complete-HGNC-21Apr2016.csv")
#
# ----------
#

write_normalization_table <- function(mq.design = "./experimentalDesignTemplate.txt",
                                   mq.evidence = "./evidence.txt",
                                   normf.median = 1e7,
                                   mq.version = "1.5",
                                   out.normf = "normfactors.txt"

                                   ) {

#   ---- experimentalDesignTemplate.txt file is read in
    tab.design <- read.delim(mq.design, header = TRUE)
    experiments <- unique(tab.design[,"Experiment"])
    no.raw <- nrow(tab.design)

    ## calculating normalization factor

    cat("Calculating normalization factors...")

    evidence <- read.delim(mq.evidence, stringsAsFactors=FALSE)

    ## check for missing Experiment column
    exp.col <- "Experiment" %in% colnames(evidence)

    if (!exp.col) {
        stop("\n", "The evidence file is missing an 'Experiment' column. Please add.", "\n")
    } else {
        if (mq.version == "1.4") {
            evidence.ok <- subset(evidence, Reverse == "" & 
                              Contaminant == "" & Intensity > 0 & MS.MS.Count > 0)
        } else {
            names(evidence) <- str_replace(names(evidence),"MS.MS..ount", "MS.MS.Count")
            evidence.ok <- subset(evidence, Reverse == "" & Potential.contaminant == "" & 
                                    Intensity > 0 & MS.MS.Count > 0)
        }

        m <- sapply(as.character(experiments), function(x) median(evidence.ok$Intensity[evidence.ok$Experiment == x]))
        normf <- if(normf.median > 0) normf.median/m else median(m)/m
        write.table(data.frame(experiment=make.names(experiments), normf), 
                file=out.normf, sep='\t', row.names=FALSE, quote=FALSE, qmethod="double")

        cat(" Done.\n")
    }

}

#
# ------------
#

extract_mq_version <- function ( pp.file ){

    first.line <- readLines(pp.file, n=1) 
    c.names <- unlist(strsplit(first.line,split="\t"))
    found <- grep("ontaminant",c.names)
    if (length(found)==0){
        version <- "1.5"
    }
    else {
        if (c.names[found] == "Potential contaminant" | c.names[found] == "Potential.contaminant"){
            version <- "1.5"
        }
        else {
            version <- "1.4"
        }
    }
    cat("Max Quant version detected as :",version,".\n")
    return(version)
}

#
# --------
#

to_kinase <- function(g.names_HGNC, kinases, sorted = TRUE) {

    tmp <- g.names_HGNC

    for (i in 1:length(tmp)) {
        if (tmp[i] != "") {
            gs <- unlist(strsplit(tmp[i], ";"))

            ind <- match(gs, kinases);

            a  <- NULL
            for (g in gs) {
                if (g %in% kinases) {
                    a <- c(a, g)
                }
            }

            if (is.null(a)) {
                tmp[i]  <-  ""
            }
            else {
                if (sorted) {
                    tmp[i]  <-  paste(sort(a), collapse = ";")
                } else {
                    tmp[i]  <-  paste(a, collapse = ";")
                }
            }
        }
    }

    return(tmp)
}

#
# --------
#

make_gene_list <- function (uids,res.type,mapping=UniProt_data){

    genes <- match(uids,mapping$UniProt.Entry)
    matched <- !is.na(genes)
    if (sum(matched)>0){
        u.genes <- unique(mapping[[ res.type ]][genes[matched]])
        retval <- paste0(u.genes,collapse=";")
    } else {
        retval <- ""
    }

    return(retval)
    
}

#
# --------
#

mapping_gn <- function(uids) {

    g.names <- lapply(strsplit(uids,split=";"),make_gene_list,res.type="gene")

    return(unlist(g.names))
}

#
# --------
#

mapping_gn_HGNC <- function(uids) {

    g.names <- lapply(strsplit(uids,split=";"),make_gene_list,res.type="HGNC")

    return(unlist(g.names))
}

#
# -----
#

loadf <- function(f) {

        env <- new.env()
        obj <- load(f, env)
        env[[obj]]
}

#
# -----
#

write_pp_report <- function(file.in = "./modificationSpecificPeptides.txt",
                                     mq.design = "./experimentalDesignTemplate.txt",
                                     mq.evidence = "./evidence.txt",
                                     out.normf = "./normfactors.txt",
                                     mq.version = "1.5",
                                     evquant.command = evquant.path,
                                     kinase.file = "../../Manning.Kinases_12_2015.txt",
                                     mapping = NA) {

    cat("Preparing pp-Report file.\n")

    tab.pg <- read.delim(file.in,stringsAsFactors=FALSE)

    if (mq.version == "1.4") {
        tab.tmp <- subset(tab.pg, Reverse == "" & Contaminant == "" & Intensity > 0 & Phospho..STY. > 0)
    } else {
        tab.tmp <- subset(tab.pg, Reverse == "" & Potential.contaminant == "" & Intensity > 0 & Phospho..STY. > 0)
    }

    cat("Mapping genes...\n")

    kinases <- as.character(loadf(kinase.file)[,1])

    tab.pg.selected <- tab.tmp[order(-tab.tmp[,"Intensity"]),]

    cat(" * Step 1\n")
    g.names <- paste0("=\"", mapping_gn(tab.pg.selected$Proteins),  "\"")

    cat(" * Step 2\n")
    g.names_HGNC <- mapping_gn_HGNC(tab.pg.selected$Proteins)

    cat(" * Step 3\n")
    g.names_kinase <- to_kinase(g.names_HGNC, kinases)


    tab.design <- read.delim(mq.design, header = TRUE)

    experiments <- gsub("[\\+-/ ()]", ".", unique(tab.design[,"Experiment"]))

    col.int <- paste0("Intensity.", experiments)

    # normalize
    normf <- read.delim(out.normf, header = TRUE)
    rownames(normf) <- paste0("Intensity.", normf[,"experiment"])

    expr <- tab.pg.selected[, col.int]
    for (ex in colnames(expr)) {
        expr[, ex] <- expr[, ex] * normf[ex, "normf"]
    }
    colnames(expr) <- paste0("Norm.", colnames(expr))

    Protein.Groups <- paste0("=\"", as.character(tab.pg.selected[,"Protein.Groups"]), "\"")

    Gene.names <- paste0("=\"", as.character(tab.pg.selected[,"Gene.Names"]), "\"")

    Proteins <- tab.pg.selected[, c("Proteins")]

    cat(" Ready.\n")

# ----------
canonize_evidence_names <- function(ev){

    old.colnames <- colnames(ev)

    idx <- grep("MS\\.MS\\.Count",colnames(ev),ignore.case=T)
    old.colnames[ idx ] <- "MS.MS.Count"

    return(old.colnames)
}
# ----------

    cat("Creating phospho peptide count for peptide table ...\n")

    ev <- fread(mq.evidence,check.names=TRUE)

    setnames(ev,canonize_evidence_names(ev))

    cat("Processing evidence file..")
    evpp <- ev[Phospho..STY. > 0, .(Mod..peptide.ID, Experiment,  MS.MS.Count)]
    cat("Ready.\n")

    pp.count_long <- evpp[, sum(MS.MS.Count), by =  .(Mod..peptide.ID, Experiment)]

    setnames(pp.count_long, c("ppModPeptide.ID", "Sample", "MS.MS.Count"))
    pp.count_wide <- dcast(pp.count_long, ppModPeptide.ID ~ Sample, value.var = "MS.MS.Count")
    pp.count_wide[is.na(pp.count_wide)] <- 0
    names(pp.count_wide)[-1] <- paste0("OPL.ppCount.", names(pp.count_wide)[-1])


    cat("Done.\n")

    tab.out <- cbind(tab.pg.selected[, c("id", "Sequence", "Modifications",
                                         "Mass","Proteins")],
                     list(Gene.names = g.names,
                          HGNCplus.gene.names = paste0("=\"", g.names_HGNC, "\""),
                          Sorted.kinase.names = paste0("=\"", g.names_kinase,  "\""),
                          MQ.gene.names = Gene.names),
                     tab.pg.selected[, c("Protein.Names",
                                         "Phospho..STY.",
                                         "MS.MS.Count",
                                         "Score",
                                         "Phospho..STY..site.IDs",
                                         col.int)],
                     expr)

    colnames(tab.out) <- str_replace(colnames(tab.out),"id","ppModPeptide.ID")

    tab.out <- merge(tab.out, pp.count_wide, by = "ppModPeptide.ID", all.x = TRUE)

    write.table(tab.out,
                file = "pp-Peptide-report.txt",
                sep = "\t",
                row.names = FALSE,
                quote = FALSE,
                qmethod = "double")

}

#
# ---------------
#

cleanup_files <- function (flist) {

    result <- file.remove(flist)
    if ( sum(result) != length(flist)){
        cat("Failed removing ",length(flist) - sum(result)," files.\n")
    }

}

#
# ---------------
#
cat("Preparing gene mapping structures...")

load(paste0(install.location,"/","UniProt_data_08Jun2016.Rdata"))
hgnc <- read.delim(complete.hgnc.path,header=TRUE,sep=";",stringsAsFactors=FALSE)
colnames(UniProt_data)[1] <- "accession"
colnames(UniProt_data)[3] <- "gene"
UniProt_data$HGNC <- rep("",nrow(UniProt_data))
hg.match.indx <- match(UniProt_data$UniProt.Entry,hgnc$uniprot_ids)
hg.matched <- !is.na(hg.match.indx)
UniProt_data$HGNC[hg.matched] <- hgnc$symbol[ hg.match.indx[hg.matched] ]
still.unmapped <-  UniProt_data$HGNC == ""
hgnc.index <- match( UniProt_data$gene, hgnc$symbol )
hgnc.ok <- !is.na(hgnc.index)
shift.id <- still.unmapped & hgnc.ok
UniProt_data$HGNC[shift.id] <- UniProt_data$gene[shift.id]

cat(" Ready.\n")

mqver <-  extract_mq_version("./modificationSpecificPeptides.txt")

write_normalization_table(mq.version=mqver)

write_pp_report(kinase.file = manning.kinase.file, mapping = UniProt_data ,mq.version=mqver)

cat("Cleaning up temporary files.\n")
temp.files <- "normfactors.txt" 
cleanup_files(temp.files)
cat("Ready.\n")
