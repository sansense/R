library(ReadqPCR)
library(NormqPCR)
library(reshape)
library(ggplot2)
library(xlsx)


###################################################################################
# This convertst the data-matrix in pcrBatch object type.
# pcrMatrix => genes/mirs are in rows; samples are in columns
get.pcrBatch = function(pcrMatrix) {  
  # Transform the pcrMatrix in the format of pcr type files
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
  
  return(pcr)
}


###################################################################################
# Runs the Vandesompele hkg Analysis, returns the hkg-rankings
# pcrBatch = batch-PCR object; allHKGs is vector containg all the potential HKGs
hkgAnalysis = function (pcrBatch, allHKGs) {
  
  pcr=pcrBatch
  
  # get "only HKGs" for the HKG analysis
  pcr4hkg=pcr[allHKGs, ]
  
  #Ranking of hkgs
  hkgsRanking=selectHKs(pcr4hkg,method="geNorm",Symbol=featureNames(pcr4hkg),log=TRUE)  
  
  # This is the msg to print to user
  msg = "================================================
  # run plotAllHKGs = function(pcrBatch, hkgsRanking) on data to see the plotting
  [Vandesompele et. al.] recommend a cut-off value of 0.15 for the pairwise variation. Below this bound the inclusion of an additional housekeeping gene is not required. 
  # This is an arbitrary cutt-off. So we would plot all the gene-combinations (2:max_len) to see for which length, the HKGs are stable graphically.
  # [Vandesompele] http://genomebiology.com/2002/3/7/research/0034
  # DOI: doi:10.1186/gb-2002-3-7-research0034  
================================================"
  cat(msg)
  
  return(hkgsRanking)
}  


###################################################################################
# plots all the HKG combinations for visualization
# writes AM and GM corrected files (pcr.arith.exprs.xlsx etc.) to disk 
plotAllHKGs = function(pcrBatch, hkgsRanking) {  
  
  pcr = pcrBatch
  allHKGs = hkgsRanking$ranking # list of all HKGs
  
  hkgPlots.arith = list() # AM-corrected plots
  hkgPlots.geom = list() # GM-corrected plots
  
  #totHKGs = length(hkgsRanking$ranking)
  for (i in 1:(length(allHKGs) - 1)) {
    hkgs = hkgsRanking$ranking[1:(i+1)]  # minimum 2 KHGs is reqd for this method
    
    # Aritmetic mean of ge-norm
    pcr.arith = deltaCq(pcr, hkgs=hkgs, calc="arith", combine=TRUE)
    pcr.arith.exprs = exprs(pcr.arith) # AM-corrected exprs  
    write.xlsx(pcr.arith.exprs, file="pcr.arith.exprs.xlsx", sheetName=paste(toString(i+1), "HKG"), append=(i>1))
    hkgPlots.arith[[i]] = plotHKG(pcr.arith.exprs[allHKGs,], hkgs)
    
    # Gometric mean of ge-norm
    pcr.geom=deltaCq(pcr,hkgs=hkgs,calc="geom",combine=TRUE)
    pcr.geom.exprs = exprs(pcr.geom) # GM-corrected exprs 
    write.xlsx(pcr.geom.exprs, file="pcr.geom.exprs.xlsx", sheetName=paste(toString(i+1), "HKG"), append=(i>1))
    hkgPlots.geom[[i]] = plotHKG(pcr.geom.exprs[allHKGs,], hkgs)      
  }
  
  # Arrange the arith and geom on the same grid
  require("gridExtra")
  len = length(hkgPlots.arith)
  
  # make the variations as titles of plot
  variation = rev(hkgsRanking$variation) # The variations are in the reverese order of HKGs!!!
  main = paste("variation", names(variation), "=", variation)
  
  pdf("hkgPlots.pdf")
  
  # First plot the uncorrected one
  uncorrectedPlot = plotHKG(exprs(pcr)[allHKGs,], "Original / Uncorrected")
  print(uncorrectedPlot)         
  
  for(i in 1:len) {
    p1 = hkgPlots.arith[[i]] 
    p2 = hkgPlots.geom[[i]]
    
    # take any one of "bottom positioned" legends.
    my.legend = g_legend(p1 + theme(legend.position="bottom")) 
    
    # arrange AM-corrected and GM-corrected plots side-by-side
    arrangePlots = arrangeGrob(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), nrow=1)     
  
    # arrange on the grid. Legend goes on the top row (with the proportion of rows 10:1)
    grid.arrange(my.legend, arrangePlots, nrow=2,heights=c(1, 10), main=main[i])
  }
  
  dev.off()
} 

###################################################################################
# plots a single HKG
# hkgMatrix = The hkgMatrix to plot
# hkg = the names of hkgs. This is used to create the title of plot.
plotHKG = function(hkgMatrix,hkg) {  
  hm=melt(hkgMatrix, id.vars=1:ncol(hkgMatrix), varnames=c("HKG", "sample"))
  
  p1 = ggplot(hm, aes(sample, value, color=HKG)) 
  p2 = p1 + geom_line(aes(group=HKG)) # line  
  p3 = p2 + geom_point() # also print the points
  p4 = p3 + scale_y_continuous(breaks=-40:40) # y-ticks at equal interval, easy to visualize
  
  # x-axis legend is bold, 90-degree rotated, and the legend is in the middle of tick (vjust=0.5)
  p5 = p4 + theme(axis.text.x = element_text(face="bold", angle=90, vjust=0.5))
  
  # unimplemented: can change the main legened background and positiong (for example, right on the graph)
  # p5 = p4 + theme(axis.text.x = element_text(face="bold", angle=90, vjust=0.5), legend.position=c(0.5,.8), legend.background = element_rect(fill = NA, colour = NA)) 
  
  # creating the main title of plot
  argname = substitute(hkgMatrix) # to create the "main title", take the "name" of matrix passed 
  main = paste(c(argname, hkg), collapse=" \n") # "main title is the matrix-name + hkg-names" 
  p6 = p5 + ggtitle(main) 
  return(p6)
}


###################################################################################
# extracts the legend of a ggplot
# http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

