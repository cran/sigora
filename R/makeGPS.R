require(slam)


#' Create your own Signature Object.
#'
#' Given a repository of gene-pathway associations either in a tab delimited
#' file with three columns (pathwayID,pathway Description,Gene) or a
#' corresponding dataframe, this function identifies all Gene Pair Signatures
#' (pairs of genes that are as a combination unique to a single pathway) and
#' Pathway Unique Genes (genes that are uniquely associated with a single
#' pathway) and stores them in a format that is usable by \code{sigora}.
#' Please also see the "details" and "note" sections below.
#'
#' The primary purpose of \code{makeGPS} is to convert a user-supplied
#' gene-pathway association table to a repository of weighted Gene Pair
#' Signatures (GPS) that are unique features of pathways. Such GPS can than be
#' used for signature (gene-pair) based analyses using \code{sigora}.
#' Additionally, the resulting object also retains the original "single
#' gene"-"pathway" associations for the purpose of followup analyses, such as
#' comparison of sigora-results to traditional methods. \code{ora} is an
#' implementation of the traditional (individual gene) Overrepresentation
#' Analysis.
#' @export
#' @param pathwayTable A data frame describing gene-pathway associations in
#' following format: pathwayID,pathwayName,Gene. Either pathwayTable or fn
#' should be provided.
#' @param fn Where to find the repository.Either pathwayTable or fn should be
#' provided.
#' @param maxLevels For hierarchical repositories, the number of levels to
#' consider.
#' @param saveFile Where to save the object as an rda file.
#' @param repoName Repository name.
#' @param maxFunperGene A cutoff threshold, genes with more than this number of
#' associated pathways are excluded to speed up the GPS identification process.
#' @param maxGenesperPathway A cutoff threshold, pathways with more than this
#' number of associated genes are excluded to speed up the GPS identification
#' process.
#' @param minGenesperPathway A cutoff threshold, pathways with less than this
#' number of associated genes are excluded to speed up the GPS identification
#' process.
#' @return A GPS repository, to be used by \code{sigora} and \code{ora.}
#' @note This function relies on package \code{slam}, which should be installed
#' from CRAN.  It is fairly memory intensive, and it is recommended to be run
#' on a machine with at least 6GB of RAM. Also, make sure to save and reuse the
#' resulting GPS repository in future analyses!
#' @seealso \code{\link{sigora}}, \code{\link{sigora-package}}
#' @references Foroushani AB, Brinkman FS and Lynn DJ
#' (2013).\dQuote{Pathway-GPS and SIGORA: identifying relevant pathways based
#' on the over-representation of their gene-pair signatures.}\emph{PeerJ},
#' \bold{1}
#' @keywords functions
#' @examples
#'
#' data(nciTable); data(idmap)
#' ## what the input looks like:
#' head(nciTable)
#' ## create a SigObject. use the saveFile parameter for reuse.
#' \donttest{
#' nciH<-makeGPS(pathwayTable=load_data('nciTable'))
#' ils<-grep("^IL",idmap[,"Symbol"],value=TRUE)
#' ilnci<-sigora(queryList=ils,GPSrepo=nciH,level=3)
#' }
makeGPS <-
  function(pathwayTable = NULL,
           fn = NULL,
           maxLevels = 5,
           saveFile = NULL,
           repoName = 'userrepo',
           maxFunperGene = 100,
           maxGenesperPathway = 500,
           minGenesperPathway = 10) {
    ##`pathwayTable: a dataframe with three columns: pathwayId, ,pathwayName,gene
    ##`fn : tab delimited file with three columns in the following order tab delimited pathway ids, pathway names, genes.
    ##` saveFile: where to store the results. (as rda file)
    ##` repoName: the repository name. Eg 'KEGG2016'
    ##` maxFunperGene: a cutoff threshold, genes with more than this number of associated pathways are excluded to speed up the GPS identification process.
    ##` maxGenesperPathway: a cutoff threshold, pathways with more than this number of associated genes  are excluded to speed up the GPS identification process.
    ##` minGenesperPathway: a cutoff threshold, pathways with less than this number of associated genes  are excluded to speed up the GPS identification process.
    if (is.null(pathwayTable)) {
      fG <-
        utils::read.table(fn,
                   header = TRUE,
                   sep = '\t',
                   quote = '@')
    } else{
      fG <- pathwayTable
    }
    ##,fileEncoding     = "utf8")
    colnames(fG) <- c('pwys', 'nms', 'gns')
    fG <- unique(fG)
    valGenes <-
      names(table(fG$gns))[which(table(fG$gns) < (1 + maxFunperGene))]
    valPathways <-
      names(table(fG$pwys))[which(table(fG$pwys) < (1 + maxGenesperPathway) &
                                    table(fG$pwys) > (minGenesperPathway -
                                                        1))]

    fG <- fG[which(fG$gns %in% valGenes & fG$pwys %in% valPathways) , ]
    ##Encoding(levels(fG$pwys)) <- "latin1"
    ##levels(fG$pwys) <- iconv(levels(fG$pwys),"latin1","UTF-8")
    ## fGM<-fG[grep('MUS',fG$gns),]
    ## fG<-fG[-grep('MUS',fG$gns),]
    L1 <- NULL
    L2 <- NULL
    L3 <- NULL
    L4 <- NULL
    L5 <- NULL
    t1 <- Sys.time()
    makedictionary <- function(y1) {
      upw <- unique(y1$pwys)
      ugn <- unique(y1$gns)
      dy1 <- cbind(match(y1$pwys, upw), match(y1$gns, ugn))
      colnames(dy1) <- c('pwys', 'gns')
      res <- list(upw, ugn, dy1)
      invisible(res)
    }
    makesigs <- function(f1) {
      gs <- as.character(unique(f1$gns))
      ps <- as.character(unique(f1$pwys))
      si <- match(f1$gns, gs)
      sj <- match(f1$pwys, ps)
      sv <- rep(1, nrow(f1))
      s <-
        slam::simple_triplet_matrix(
          i = si,
          j = sj,
          v = sv,
          dimnames = list('rownames' = gs, 'colnames' = ps)
        )

      M <- slam::tcrossprod_simple_triplet_matrix(s)

      ##M<-(m%*%t(m))
      rownames(M) <- gs
      colnames(M) <- gs
      ##rm(s)
      PU <- which(diag(M) == 1)
      #PUG<-cbind(si[PU],sj[PU])
      m <- as.matrix(s)
      #rownames(m)<-gs
      #colnames(m)<-ps
      #browser()
      xx <- which(m[PU, , drop = FALSE] == 1, arr.ind = TRUE)[, 2]
      PUG <- cbind(PU, xx[names(PU)])
      #	PUG <-
      # cbind(PU,which(m[PU,,drop=FALSE]==1,arr.ind = TRUE)[,2] )
      GP <- which(M == 1, arr.ind = TRUE)
      GP <- GP[GP[, 1] < GP[, 2], ]
      S <- vector(length = nrow(GP))
      S1 <- PUG[match(GP[, 1], PUG[, 1]), 2]
      S2 <- PUG[match(GP[, 2], PUG[, 1]), 2]
      S <- ifelse(is.na(S1), S2, S1)
      rm(M)
      gc(TRUE)
      #    m<-as.matrix(s)
      for (h1 in which(is.na(S))) {
        S[h1] <- which(m[GP[h1, 1], ] + m[GP[h1, 2], ] == 2)
      }
      GPS <- cbind(GP, S)
      print(Sys.time() - t1)
      rm(m)
      degs <- table(as.character(f1$gns))
      pwyszs <- table(as.character(f1$pwys))
      rownames(GPS) <- NULL
      rownames(PUG) <- NULL
      invisible(list(
        'GPS' = GPS,
        'PUG' = PUG,
        'gs' = gs,
        'ps' = ps,
        'degs' = degs,
        'pwyszs' = pwyszs
      ))

    }##end internal function

    L1 <- makesigs(fG)
    fG2 <- fG[-which(fG$pwys %in% (L1$ps[unique(L1$GPS[, 'S'])])), ]
    if (nrow(fG2) > 1) {
      L2 <- makesigs(fG2)
      fG3 <- fG2[-which(fG2$pwys %in% (L2$ps[unique(L2$GPS[, 'S'])])), ]
      if (nrow(fG3) > 1) {
        L3 <- makesigs(fG3)
        fG4 <-
          fG3[-which(fG3$pwys %in% (L3$ps[unique(L3$GPS[, 'S'])])), ]
        if (nrow(fG4) > 1) {
          L4 <- makesigs(fG4)
          fG5 <-
            fG4[-which(fG4$pwys %in% (L4$ps[unique(L4$GPS[, 'S'])])), ]
          if (nrow(fG5) > 1) {
            L5 <- makesigs(fG5)
          }
        }
      }
    }
    res <- list()
    res[['origRepo']] <- makedictionary(fG)
    res[['L1']] <- L1
    res[['L2']] <- L2
    res[['L3']] <- L3
    res[['L4']] <- L4
    res[['L5']] <- L5
    res[['repoName']] <- repoName
    res[['pathwaydescriptions']] <- unique(fG[, seq_len(2)])
    res[["call"]] <- as.character(match.call())
    x1 <- as.character(repoName)
    if (!is.null(saveFile)) {
      cmd2 <- paste("save(", x1 , ", file='", saveFile, "')", sep = "")
      assign(x1, res)
      eval(parse(text = cmd2))
    }
    ##save(rp,file=saveFile)
    invisible(res)
  }
