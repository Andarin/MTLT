# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEEstLiterature
# -----------------------------------------------------------------------------
# Description:   Creates 2 plots which compare estimates of chromosomal re-
#                arrangements between estimates in the thesis (Stat. and 
#                 Num. Est.) and estimates from different papers
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        2 Plots comparing recent estimates with the literature
# -----------------------------------------------------------------------------
# Keywords:      data mining, preprocessing, descriptive-statistics, tree, 
#                discrete
# -----------------------------------------------------------------------------
# See also:      -
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEEstLiterature.RData 
# -----------------------------------------------------------------------------

# Load package
require(ggplot2)

# Load data
load('COAPMODEEstLiterature.RData')

ggplot(estimates.literature, aes(species, inv)) + 
  geom_bar(aes(fill=paper), position="dodge", stat="identity") +
  xlab('') +
  ylab('Number of inversions') +
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(colour = "white"), 
        legend.title=element_blank(), 
        axis.text=element_text(size = 14), 
        axis.title = element_text(size = 14), 
        axis.title.y = element_text(vjust = 1.5),
        axis.title.x = element_text(vjust = 0.2),
        plot.title = element_text(size = 14, vjust=2))

ggplot(estimates.literature, aes(species, transl)) + 
  geom_bar(aes(fill=paper), position="dodge", stat="identity") +
  xlab('') +
  ylab('Number of translocations') +
  theme(axis.line = element_line(colour = "black"), 
        panel.background = element_rect(fill = "white"), 
        panel.grid = element_line(colour = "white"), 
        legend.title=element_blank(), 
        axis.text=element_text(size = 14), 
        axis.title = element_text(size = 14), 
        axis.title.y = element_text(vjust = 1.5),
        axis.title.x = element_text(vjust = 0.2),
        plot.title = element_text(size = 14, vjust=2))