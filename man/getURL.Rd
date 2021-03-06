\name{getURL}
\alias{getURL}

\title{
Highlight the relevant genes for a specific pathway in its pathway diagram
}
\description{
This function highlights the genes involved in the present GPS for a pathway of interest in its diagram. Please note that this functionality is only implemented for results of Reactome or KEGG based analyses.
}
\usage{
getURL(yy, i)
}
\arguments{
  \item{yy}{
A sigora analysis result object (created by \code{sigora}).
}
  \item{i}{
The rank position of the pathway of interest in summary_results.
}
}
\value{
The URL of the pathway diagram, where the relevant genes from your original query list are highlighted.
}
\seealso{
\code{\link{sigora}}
}
\examples{
a1<-genesFromRandomPathways(seed=12345,kegH,3,50)
## originally selected pathways:\cr
a1[["selectedPathways"]]
## what are the genes
a1[["genes"]]
## sigora's results with this input:\cr
sigoraRes <- sigora(GPSrepo =kegH, queryList = a1[["genes"]],level = 2)
## Diagram for the most significant result, where the genes from our list are highlighted in red:
getURL(sigoraRes,1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{functions}

