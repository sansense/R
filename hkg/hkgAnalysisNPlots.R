library(ReadqPCR)
library(NormqPCR)
library(reshape)
library(ggplot2)
library(xlsx)

### Tomorrow DO IT ####
# 1. ordered exprs are useless. Remove them!!!
# 2. Make different methods for different part of the current analysis.
##################


# Runs the Vandesompele hkg Analysis, returns the hkg-rankings
# ourHKGs is vector containg all the potential HKGs
# Creates plots for hkgs from which could be chosen the best hkg
# writes the arith and geom corrected exprs from an excel file.

# pcrMatrix => 1st row is genes/mirs, samples in column 2-onwards
hkgAnalysisNPlots = function(pcrMatrix, ourHKGs) {
  # format the pcrMatrix. The row-names should be gene_name
  rownames(pcrMatrix) = pcrMatrix[,1]
  label = colnames(pcrMatrix)[1] # gene/mir/X/...
  pcrMatrix = pcrMatrix[,-1]
  
  # Transform it in the format of  pcr type files
  # pcr-type => Sampe, Detector, Cq cols only
  # values will be written to file, then that file is read by read.qPCR!
  
  for (i in 1:ncol(pcrMatrix)) {
    xx=cbind(names(pcrMatrix)[i],row.names(pcrMatrix),pcrMatrix[i])
    names(xx)=c("Sample","Detector","Cq")
    if(i == 1) {
      qpcr=xx
    } else {
      qpcr=rbind(qpcr,xx)
    }
  }
  
  # note that you dont want the row-names which are genenames
  write.table(qpcr,"qpcr.txt",sep="\t",row.names=FALSE)
  
  # read the PCR
  pcr=read.qPCR("qpcr.txt")
  
  # get "only HKGs" for the HKG analysis
  pcr4hkg=pcr[ourHKGs, ]
  
  #Ranking of hkgs
  hkgsRanking=selectHKs(pcr4hkg,method="geNorm",Symbol=featureNames(pcr4hkg),log=TRUE)
  #hkgsRanking
  
  # This is the msg to print to user
  msg = "================================================
  # [Vandesompele et. al.] recommend a cut-off value of 0.15 for the pairwise variation. Below this bound the inclusion of an additional housekeeping gene is not required. 
  # This is an arbitrary cutt-off. So we would plot all the gene-combinations (2:max_len) to see for which length, the HKGs are stable graphically.
  # [Vandesompele] http://genomebiology.com/2002/3/7/research/0034
  # DOI: doi:10.1186/gb-2002-3-7-research0034  
================================================"
  cat(msg)
  
  hkgPlots.arith = list() # AM-corrected plots
  hkgPlots.geom = list() # GM-corrected plots
  
  totHKGs = length(hkgsRanking$ranking)
  for (i in 1:(totHKGs-1)) {
    hkgs = hkgsRanking$ranking[1:(i+1)]  # minimum 2 KHGs is reqd for this method
    # cat(i, hkgs, "\n")
    
    # Aritmetic mean of ge-norm
    pcr.arith = deltaCq(pcr, hkgs=hkgs, calc="arith", combine=TRUE)
    pcr.arith.exprs = ordered_exprs(exprs(pcr.arith), pcrMatrix) # get the ordered-expression value    
    write.xlsx(pcr.arith.exprs, file="pcr.arith.exprs.xlsx", sheetName=paste(toString(i+1), "HKG"), append=(i>1))
    hkgPlots.arith[[i]] = plotHKGs(pcr.arith.exprs[ourHKGs,], hkgs)
    
    # Gometric mean of ge-norm
    pcr.geom=deltaCq(pcr,hkgs=hkgs,calc="geom",combine=TRUE)
    pcr.geom.exprs = ordered_exprs(exprs(pcr.geom), pcrMatrix) # get the ordered-expression value
    write.xlsx(pcr.geom.exprs, file="pcr.geom.exprs.xlsx", sheetName=paste(toString(i+1), "HKG"), append=(i>1))
    hkgPlots.geom[[i]] = plotHKGs(pcr.geom.exprs[ourHKGs,], hkgs)      
  }
  
  # Arrange the arith and geom on the same grid
  require("gridExtra")
  len = length(hkgPlots.arith)
  
  # make the variations as titles of plot
  variation = rev(hkgsRanking$variation) # The variations are in the reverese order of HKGs
  main = paste("variation", names(variation), "=", variation)
  
  pdf("hkgPlots.pdf")
  # First plot uncorrected one
  uncorrectedPlot = plotHKGs(exprs(pcr4hkg), "Uncorrected/Original")
  print(uncorrectedPlot)         
  
  for(i in 1:len) {
    my.legend = g_legend(hkgPlots.arith[[i]] + theme(legend.position="bottom")) # take any one of "bottom positioned" legends.
    p1 = hkgPlots.arith[[i]] 
    p2 = hkgPlots.geom[[i]]
    grid.arrange(my.legend, arrangeGrob(p1 + theme(legend.position="none"),
                                       p2 + theme(legend.position="none"),
                                       nrow=1),  nrow=2,heights=c(1, 10), main=main[i])
    # arrange the AM and GM corrected plots side-by-side 
    #grid.arrange(hkgPlots.arith[[i]], hkgPlots.geom[[i]], ncol=2, main=main[i])
  }
  dev.off()
  
  return(hkgsRanking)
}

# preserve the right order using the col-order of original matrix
ordered_exprs = function (pcrUnsorted, orig) {
  pcrSorted=pcrUnsorted[order(match(rownames(pcrUnsorted), rownames(orig))), ] # row sorting  
  pcrSorted=pcrSorted[ ,order(match(colnames(pcrSorted), colnames(orig)))] # col sorting
  return(pcrSorted)
}

# plots HKGs
# hkgMatrix = The hkgMatrix to plot
# hkg = the names of hkgs. This is used to create the title of plot.
plotHKGs = function(hkgMatrix,hkg) {  
  hm=melt(hkgMatrix, id.vars=1:ncol(hkgMatrix), varnames=c("HKG", "sample"))
  
  p1 = ggplot(hm, aes(sample, value, color=HKG)) 
  p2 = p1 + geom_line(aes(group=HKG)) # line  
  p3 = p2 + geom_point() # also print points
  #p4 = p3 + scale_y_discrete() # y-ticks at equal interval, easy to visualize
  p4 = p3 + scale_y_continuous(breaks=-40:40) # y-ticks at equal interval, easy to visualize
  
  p5 = p4 + theme(axis.text.x = element_text(face="bold", angle=90, vjust=0.5))#, legend.position=c(0.5,.8)) # x-axis legend is bold, 90-degree rotated, and the legend is in the middle of tick (vjust=0.5)
  
  #p5 = p4 + theme(axis.text.x = element_text(face="bold", angle=90, vjust=0.5), legend.position=c(0.5,.8), legend.background = element_rect(fill = NA, colour = NA)) # x-axis legendd is bold, 90-degree rotated, and the legend is in the middle of tick (vjust=0.5)
  
  
  
  # creating the main title of plot
  argname = substitute(hkgMatrix) # to create the "main title", take the "name" of matrix passed 
  main = paste(c(argname, hkg), collapse=" \n") # "main title is the matrix-name + hkg-names" 
  p6 = p5 + ggtitle(main) 
  return(p6)
}

# extracts the legend of a ggplot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

