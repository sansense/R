# formats the pcr-data and hkgs with proper names
# data.path and hkg.path *must* be  provided
# By default, median.corrects the duplicated data.

format.PCR.data = function (data.path="", hkg.path="", median.correct.data="Y") {
  dat = read.table(data.path, sep="\t", header=T)
  hkg = scan(file = hkg.path, what="", sep="\n")

  # make proper names
  colnames(dat) = make.names(colnames(dat))
  dat[,1] = make.names(dat[,1])  # title, whether it is gene/mir/X/...
  hkg = make.names(hkg)

  if(median.correct.data == "Y") { # requires median correcting of the duplicates
    cat("median correcting the duplicate-data")
    dat = medianCorrectDuplicates(dat) 
  }
  
  return(list(pcr.data=dat, hkgs=hkg))  
}

# Since the 1stCol can have duplicated data, we need to median correct all the dups before proceeding for any analysis
medianCorrectDuplicates = function (dat) {
  genenames = dat[,1] # first col is gene/mir name
  label = colnames(dat)[1] # what is the label of first col? gene/mir/X/...
  
  dups = unique(genenames[which(duplicated(genenames))])
  dups = as.vector(dups)
  
  for (i in dups) {
    dups.data = dat[(dat[,1] == i), ] 
    # median of dups
    dups.data.median = sapply(dups.data[,-1], median)
    
    names(i)  = label; # set the correct label, which is needed for rbinding later with original format of df
    # dat.dups.median = as.data.frame(t(c(i, dups.data.median))) # this doesnwt work as c coereces everything to char
    dat.dups.median = data.frame(t(i), t(dups.data.median)) # # vectros are like row-data with only 1Col. Make them as col-data.
    
    dat.rest = dat[!(dat[,1] == i), ] # rest of dat data
    dat = rbind(dat.rest, dat.dups.median)
  }
  
  return(dat)
}

