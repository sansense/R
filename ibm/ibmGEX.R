# Auto and Allo pts
auto="P05|P06|P08"
allo="P02|P03|P04|P07"

# time points in months
tPoint=c(1, 3, 10, 12)


# write the result for individual months for a particular kind of pt (auto| allo| etc.)
# also write the combined summary (resCombined) results in another file (only summary + annot; no pt data)
writeResults = function (pt.type, mo, select.genes = NULL) {
  resCombined = ibm.exprs[, NULL] # get only the rownames to intialize the results
  for (i in mo) {
    label = paste(substitute(pt.type), i, "m", sep="")  
    cat("prcoessing", label, "\n")
    
    pts = my.pts(ibm.exprs, pt.type, as.numeric(i))
    res.ttest = my.ttest(pts, label)
    if(is.null(res.ttest)) next;
    res.ttest = apply(res.ttest, 2, unlist) #unlist the list item so that it can be printed later
    resCombined = cbind(resCombined, res.ttest)  # combined result to print at end
    res = cbind(ibm.annot, pts, res.ttest)
    #res1 = data.frame(lapply(res, function(x) (unlist(x)))) # remove the list element so as to write it properly
    write.table(res, file=paste(label, ".txt", sep=""), sep="\t", col.names=NA)
    
    # if genes has to be selected
    if(!is.null(select.genes)) {
      res.selected = resSelected(res, select.genes)
      write.table(res.selected, file=paste(label, ".selected.txt", sep=""), sep="\t", col.names=NA)
    }
  }
  
  fname = paste(substitute(pt.type), paste(mo, collapse="_"), "m.combined.txt", sep="")
  write.table(resCombined, file=fname, sep="\t", col.names=NA)
  
}   


# returns only the selected pts according to "auto| allo| any_other_set" for given months
# x = df containing the pts info
# pts = which pts are to be returned (default = * =>  ALL)
# mo = a selection of ("pre", "all" or a vector of months)
my.pts = function (x, pts="*", mo) {
  pt.ids = grep(pts, colnames(x), value=T)   
  pt.sel = NULL
  for(i in 1:length(mo)) {    
    if(mo[i] == "all") {
      pt.sel = pt.ids
      break
    }  
    if(mo[i] == "pre") {
      m_ = mo[i];  
    } else if (mo[i]>=10) {
      m_ = paste(mo[i], "[CT]", sep="")
    } else { # 0 < mo <10
      m_=paste("0", mo[i], "[CT]", sep="") #padding because initial digit is not present
    }    
    if(is.null("pt.sel")) { 
      pt.sel = grep(m_, pt.ids, value=T)
    } else{
      pt.sel = c(pt.sel, grep(m_, pt.ids, value=T))
    }
  }
  
  return (x[, pt.sel, drop=F]) 
}


# paired ttest for the data; returns statistic, p-val, adj.p.val and estimate (as list)
# dat = data; name=name you want to append at the end of the various stats
my.ttest = function (dat, name="") {
  ct=dat[,grep("C$", colnames(dat)), drop=F] #drop=F allows even one col matrix to remain matrix
  tx=dat[,grep("T$", colnames(dat)), drop=F]
  cat("using CT: ", colnames(ct), "\n")
  cat("using TX: ", colnames(tx), "\n")
  if(length(colnames(ct)) < 2) { # ttest needs at least 2 observations.
    #myt=ct[,-(1:ncol(ct))] # return the matrix with same rownams but zero cols
    cat("Skipping: not possible to do ttest for <2 data point\n")
    return(NULL)
  }
  myt=t(sapply(1:nrow(ct), function(i) t.test(tx[i,], ct[i,], paired=T)))
  rownames(myt) = rownames(dat)
  myt=myt[,c("statistic", "p.value", "estimate")]
  adj.p.value=p.adjust(myt[, "p.value"], method="BH") #BH adjustment
  myt=cbind(myt, adj.p.value)
  colnames(myt) = paste(colnames(myt), toupper(name), sep=".")
  return(myt)
}

# returns selected genes only from the whole (or annotation)
# res = result (w/ or only the annot)
# to.select.genes = gene names (vector) to be grep'ed. In our case, a space will be added before the gn, so as to look for a word begining with that gn. Normal grep can match even in the middle of word. Next iteration will be a better way to search it using regex.

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

