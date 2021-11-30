## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library("sigora")

## -----------------------------------------------------------------------------
data(kegH)
set.seed(seed = 12345)
a1 <- genesFromRandomPathways( kegH, 3, 50)
## originally selected pathways:
a1[["selectedPathways"]] ## what are the genes
a1[["genes"]]
## Traditional ora identifies dozens of statistically significant pathways!

head(ora(a1[["genes"]], kegH))
## Now let us try sigora with the same input:

sigoraRes <-
  sigora(GPSrepo = kegH,
         queryList = a1[["genes"]],
         level = 4)
## Again, the three originally selected pathways were:
a1[["selectedPathways"]]


## -----------------------------------------------------------------------------
data(nciTable) ## what does the input look like?
head(nciTable) ## create a SigObject. use the saveFile parameter for future reuse.
data(idmap)
nciH <- makeGPS(pathwayTable = nciTable, saveFile = NULL)
ils <- grep("^IL", idmap[, "Symbol"], value = TRUE)
ilnci <- sigora(queryList = ils,
                GPSrepo = nciH,
                level = 3)


## -----------------------------------------------------------------------------
sigRes <-
  sigora(kegH,
         queryList = a1$genes,
         level = 2,
         saveFile = NULL)


## -----------------------------------------------------------------------------
data('kegH')
data('idmap')
barplot(
  table(kegH$L1$degs),
  col = "red",
  main = "distribution of number of functions per gene in KEGG human pathways.",
  ylab = "frequency",
  xlab = "number of functions per gene"
)
## creating your own GPS repository
nciH <- makeGPS(pathwayTable = load_data('nciTable'))
ils <- grep("^IL", idmap[, "Symbol"], value = TRUE)
## signature overrepresentation analysis:
sigRes.ilnci <- sigora(queryList = ils,
                       GPSrepo = nciH,
                       level = 3)


## -----------------------------------------------------------------------------
sessionInfo()

