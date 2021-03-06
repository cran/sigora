\name{genesFromRandomPathways}
\alias{genesFromRandomPathways}

\title{
Function to randomly select genes associated with randomly
pathways. 
}
\description{
This function first randomly selects a number (np) of pathways, then
randomly selects a number (ng) of genes that are associated with at
least one of the selected pathways.  The function can be used to compare
Sigora's performance to traditional overrepresentation tests.
 }
\usage{
genesFromRandomPathways(seed = 1234, GPSrepo, np, ng)
}
\arguments{
  \item{seed}{
A random seed.
}
  \item{GPSrepo}{
A signature repository (created by ..) or one of the precompiled options.
}
  \item{np}{
How many pathways to select.
}
  \item{ng}{
Number of genes to be selected.
}
}
\value{
\item{selectedPathways}{A vector containing the "np" originally selected pathways.}
\item{genes}{ A vector containing the "ng" selected genes from selectedPathways.}
}
\references{
Foroushani AB, Brinkman FS and Lynn DJ (2013).\dQuote{Pathway-GPS and SIGORA: identifying relevant pathways based on the over-representation of their gene-pair signatures.}\emph{PeerJ}, \bold{1}
}

\seealso{
\code{\link{sigora-package}}
}
\examples{
## select 50 genes from 3 human KEGG pathways
a1<-genesFromRandomPathways(seed=12345,kegH,3,50)
## originally selected pathways:
a1[["selectedPathways"]]
## what are the genes
a1[["genes"]]
## sigora's results
sigoraRes <- sigora(GPSrepo =kegH, queryList = a1[["genes"]],
        level = 4)
## compare to traditional methods results:
oraRes <- ora(a1[["genes"]],kegH)
dim(oraRes)
oraRes
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}
