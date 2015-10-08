# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEAmniotaChr
# -----------------------------------------------------------------------------
# Description:   Probability for different chromosome numbers for
#                the ancestral Amniota genome, inferred with ChromEvol 2
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        Plot shows the probability (y) for different chromosome 
#                numbers (x) for the ancestral Amniota genome
# -----------------------------------------------------------------------------
# Keywords:      plot, probability, Markov, transition probability, 
#                distribution
# -----------------------------------------------------------------------------
# See also:      COAPMODEChromEvolAIC
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEAmniotaChr.RData
# -----------------------------------------------------------------------------

# Load package
require(ggplot2)

# Load data
# Probability of different Amniota chromosome numbers according to the optimal
# parametrisation of the continous Markov-chain model with transition pro-
# bability at 0.525
load('COAPMODEAmniotaChr.RData')

ggplot() + 
    geom_segment(aes(x    =  df.nChrAmn.p$nChr[which.max(df.nChrAmn.p$P)], 
                     y    = 0, 
                     xend = df.nChrAmn.p$nChr[which.max(df.nChrAmn.p$P)], 
                     yend = max(df.nChrAmn.p$P)),
                 size     = 1.3, 
                 colour   = '#CCCCCC', 
                 linetype = 'dashed') +
    geom_line(data    = df.nChrAmn.p, 
              mapping = aes(x=nChr, y=P), 
              size    = 1.3) + 
    coord_cartesian(xlim = c(1,45)) +
    ggtitle('Posterior probability of Amniota chromsome number\nfor fission/fusion rate = 0.525') +
    xlab('Number of chromosomes for Amniota') +
    ylab('Probability') +
    theme(axis.line        = element_line(colour = "black"), 
          panel.background = element_rect(fill = "white"), 
          panel.grid       = element_line(colour = "white"), 
          legend.title     = element_blank(), 
          axis.text        = element_text(size = 14), 
          axis.title       = element_text(size = 14), 
          axis.title.y     = element_text(vjust = 1.5),
          axis.title.x     = element_text(vjust = 0.2),
          plot.title       = element_text(size = 16, vjust=2))