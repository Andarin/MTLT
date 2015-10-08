# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEEntropyScore
# -----------------------------------------------------------------------------
# Description:   For over 200 genome-genome comparisons with different phylo-
#                genetic distances, the Chi-square and G-Statistic log-p-
#                values are calculated to measure the entropy due to 
#                gene events and chromosomal rearrangements.
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        Plot showing the development of log-p-value statistics (y) for
#                Chi-square and G-Statistic for genome pairings with different
#                phylogenetic distances (x)
# -----------------------------------------------------------------------------
# Keywords:      scatterplot, discrete, distance, chi-square, Kullback-Leibler 
# -----------------------------------------------------------------------------
# See also:      -
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEEntropyScore.RData
# -----------------------------------------------------------------------------

# Load package
require(ggplot2)

# Load data
# Dataset contains information of for different 2 genome pairings
# concerning the phylogenetic distance, 
# for the genome-genome contingency table in genes the log-p-value 
# of Chi-square test and G-Test (Kullback-Leibler),
# Furthermore the p-value based on contigency table with pseudo-count
load('COAPMODEEntropyScore.RData')

# Define functions
fancy_scientific <- function(l) {
    # Function copied from Brian Diggs,
    # https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
} 

plot.label <- list(bquote(chi^2-Score),bquote(G-Score))

ggplot() + 
    geom_point(data = df.entropy.score, 
               aes(x     = phylo_distance, 
                   y     = pseudocount, 
                   group = variable, 
                   fill  = variable, 
                   color = variable, 
                   shape = variable)
               ) +
    ggtitle("Testing different scores for 189 species comparisons") +
    xlab("Phylogenetic distance in mya") +
    scale_y_continuous(name = "Log p-value", 
                       labels = fancy_scientific)+
    scale_fill_manual(name   = "Score", 
                      values = c(2, 3), 
                      breaks = c("chi2_nodup","G_nodup"), 
                      labels = plot.label) +
    scale_color_manual(name   = "Score", 
                       values = c(2, 3), 
                       breaks = c("chi2_nodup","G_nodup"), 
                       labels = plot.label) +
    scale_shape_manual(name   = "Score", 
                       values = c(21, 24), 
                       breaks = c("chi2_nodup","G_nodup"), 
                       labels = plot.label) +
    theme(axis.line        = element_line(colour = "black"), 
          panel.background = element_rect(fill = "white"), 
          panel.grid       = element_line(colour = "white"), 
          legend.title     = element_blank(), 
          axis.title       = element_text(size = 16), 
          title            = element_text(size = 12),
          axis.title       = element_text(size = 16, vjust=3), 
          title            = element_text(size = 12),
          plot.title       = element_text(vjust=2, size = 16))