
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("trackViewer")

setwd("C:/Users/SETU/Documents")
library(trackViewer)

SNP <- c(19,	58,	80,	103,	104,	133,	134,	142,	211,	212,	215,	417,	440,	446,	452,	477,	478,	484,	493,	496,	498,	501,	505,	547,	570,	613,	614,	642,	655,	679,	681,	701,	764,	796,	856,	950,	954,	969,	981,	1263,	1264)

Name <- c("T19R",	"F58V",	"D80A",	"G103D",	"W104R",	"F133L",	"Q134X",	"G142D",	"N211K",	"L212V",	"D215G",	"K417N",	"N440K",	"G446S",	"L452R",	"S477N",	"T478K",	"E484K/A",	"Q493R",	"G496S",	"Q498R",	"N501Y",	"Y505H",	"T547K",	"A570D",	"Q613H",	"D614G",	"V642F",	"H655Y",	"N679K",	"P681H/R",	"A701V",	"N764K",	"D796Y",	"N856K",	"D950N",	"Q954H",	"N969K",	"L981F",	"S1263L",	"V1264L")

Name2 <- c("S1", "Furine cleavage site", "S2", "RBD")

sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0(Name)))
features <- GRanges("chr1", IRanges(c(1, 685, 687, 238), 
                                    width=c(684, 2, 587, 293),
                                    names=paste0(Name2)))

features$fill <- c("blue", "#51C6E6", "#DFA32D", "#800000")

features$height <- c(0.02, 0.04, 0.02, 0.04)

sample.gr$color <- rep(list(c("red")), length(SNP), replace=TRUE)



frq <- c(1,	1,	1,	1,	1,	1,	1,	1,	1,	2,	1,	4,	2,	2,	1,	3,	4,	4,	3,	3,	3,	4,	3,	4,	1,	1,	4,	1,	4,	3,	5,	1,	2,	4,	2,	2,	3,	4,	4,	1,	1)

sample.gr$score<- frq
sample.gr$score
pdffile <- "rong3.pdf";
pdf(pdffile, 6, 6);
par(mar=c(2, 2, 2, 2));

lolliplot(sample.gr, features, cex=.6)


dev.off()
