###################################################################################
# extracts the legend of a ggplot
# http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
g_legend<-function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
