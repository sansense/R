library(reshape)
library(ggplot2)
library(gridExtra)

# returns plots of selected genes from result of ibm-analysis
# additionally, saves the pdf-plots in "filname".
plotSelectedGenes = function (res, to.select.genes = NULL, filename=NULL) {  
  s.genes = resSelected(res, to.select.genes)
  
  genePlots = list()
  for(i in 1:nrow(s.genes)) {
    genePlots[[i]] = GenePlotter(s.genes[i,])
  }  
  
  if(!is.null(filename)) {
    filename = paste(filename, "pdf", sep=".")
    #graphics.off()
    pdf(filename)    
    # print(genePlots) ### this prints the blank list-number on stdout (the output is correctly printed to pdf file). This is annoying if there are too many plots in the list as we see a lot of list-nums.
    for (i in 1:length(genePlots))  print(genePlots[[i]])
    dev.off()
  }
  
  return(genePlots)
}

GenePlotter = function (x) {   
  len = length(x)
  pts.data = x[, 11:(len-4)]
  annot = x[, 1:10]
  stat = x[((len-3):len)]
  
  #pts.data.m = melt(pts.data, variable_name="Patients") 
  pts.data.m = t(pts.data) # this is a matrix
  pts.data.m = data.frame(Patients = rownames(pts.data.m), value = pts.data.m[,1])
  
  pts.data.m$pt.side = substr(pts.data.m$Patients, 6, 6)  # C or T side?
  
  p1 = ggplot(pts.data.m, aes(x=Patients, y=value)) + geom_bar(stat="identity", aes(fill=pt.side))   
  p2 = p1 + geom_line(stat = "hline", yintercept = "mean", aes(colour = pt.side, group = pt.side ), size=I(2), linetype="dashed")
  
  # wraps the title with "\n"
  gene = annot$geneassignment #genename
  wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")
  p3 = p2 + ggtitle(paste(strwrap(gene,65), collapse="\n"))
  p4 = p3 + annotation_custom(grob = tableGrob(t(stat)), xmin = 0, ymin = 0)
  
  return(p4)
}

# this is the same resSelected as in ibmGEX.R
resSelected = function(res, to.select.genes) {
  to.select.genes = paste(" ", to.select.genes, sep="") # only take those where gname is the start of wor
  s.index = sapply(to.select.genes, grep, res$geneassignment, ignore.case=T)  # sapply gives the gene-name as the names of list-element. Still have to figure out how sapply works as I was expecting a matrix-sort of op. In any case, sapply is better than lapply as the names of genes are included in the op (and this seems to be the only diff between sapply and lapply!)
  
  # which genes were not found?
  s.index.length = sapply(s.index, length) # how many matches for each gene? ( == length of vector)
  missing.genes = names(s.index.length[!(s.index.length)])
  cat("missing genes:", missing.genes, "\n")
  
  s.genes = lapply(s.index, function (x) res[x,]) # here sapply gives crazy results!
  s.genes = as.data.frame(do.call(rbind, s.genes)) # convert that to a data.frame
  
  return(s.genes)
}