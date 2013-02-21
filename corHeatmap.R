# creates a heatMap of the correlation matrix
# Adapted from http://www.r-bloggers.com/simplest-possible-heatmap-with-ggplot2/
# inputp: heat = a df containig the data. The correlation between its cols will be calculated. 
# df = col of data. The rowname of df will be ignored for heatmap.
corHeatmap = function (heat) {

require("ggplot2")
require("reshape2")
require("RColorBrewer")
  
  
#calculate the correlation matrix
corData=cor(heat)
#save the colnames for later restoring the correct order
colOrig=colnames(corData) 

# Define palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# note that the colnames needs a letter in the begining otherwise the axes becomes continuous!
colnames(corData) = paste("A", (1:ncol(corData)))
rownames(corData) = paste("A", (1:ncol(corData)))

zp1 <- ggplot(melt(corData), aes(x = Var1, y = Var2, fill = value))
zp1 <- zp1 + geom_tile()
zp1 <- zp1 + scale_fill_gradientn(colours = myPalette(100))
zp1 <- zp1 + coord_equal()
zp1 <- zp1 + scale_x_discrete(label=colOrig, name="")
zp1 <- zp1 + scale_y_discrete(label=colOrig, name="")
#zp1 <- zp1 + opts(axis.text.x=theme_text(angle=-90, hjust=.5))
zp1 <- zp1 + opts(axis.text.x=theme_text(angle=-90))
zp1 <- zp1 + opts(title="Correlation Heatmap")                  
return(zp1)
}
