# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEAICChromEvol2ForFusFisRatesAmniota
# -----------------------------------------------------------------------------
# Description:   Using ChromEvol 2, a continuous Markov chain estimator,
#                the optimal transition probability for chromosome numbers
#                is inferred by searching the minimum AIC.
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        Plot showing the development of AIC (y) for continuous Markov-
#                chain models with different transition probabilities (x)
# -----------------------------------------------------------------------------
# Keywords:      plot, AIC, Markov, transition probability, optimization
# -----------------------------------------------------------------------------
# See also:      COAPMODEProbAmniotaChromosomeNumber
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEAICChromEvol2ForFusFisRatesAmniota.RData
# -----------------------------------------------------------------------------

# Load package
require(ggplot2)

# Load data
# Contains two datasets: 
# df.aic: the AIC for different values of the scale parameter
#     in the gamma distribution describing the inversion size distribution in 
#     genes
# df.rect: contains information for the range of optimal Amniota
#     starting chromosome numbers
load('COAPMODEAICChromEvol2ForFusFisRatesAmniota.RData')

plot.aic = ggplot() + 
    geom_rect(data    = df.rect, 
              mapping = aes(xmin = rect.min, 
                            xmax = rect.max, 
                            ymin = 0, 
                            ymax = 10000, 
                            fill = AmniotaChr), 
              alpha    = 0.8, 
              linetype = 'blank') +
    geom_line(data    = df.aic, 
              mapping = aes(x=GainRate, y=AIC), 
              size    = 1.3) +
    scale_x_log10(breaks       = c(0.25, 0.5, 1, 2), 
                  minor_breaks = c(0.25, 0.5, 1, 2)
                  ) +
    scale_fill_gradient(limits = c(1, 22), 
                        low    = "#FFFF77", 
                        high   = "#FF3333") +
    coord_cartesian(xlim = c(0.15, 2), 
                    ylim = c(3250,3700)
                    ) +
    ggtitle('AIC of fission/fusion rate as found by ChromEvol 2
            \n(Zones indicate most likely chromosome number for Amniota)') +
    xlab('Fission / fusion rate per million years') +
    theme(axis.line        = element_line(colour = "black"), 
          panel.background = element_rect(fill = "white"), 
          panel.grid       = element_line(colour = "white"), 
          legend.title     = element_blank(), 
          axis.text        = element_text(size = 14), 
          axis.title       = element_text(size = 14), 
          plot.title       = element_text(size = 16, vjust=2),
          legend.position  = "None")

# Add annotation to graphic: hat colour represents which starting 
# chromosome number in Amniota
pos.annotation = data.frame(c(23, 22, 21, 20, 1), 
                            c(0.18, 0.33, 0.518, 0.589, 1.2))
for (recRow in seq_len(nrow(pos.annotation))) {
    plot.aic = plot.aic + 
                   annotate("text", 
                            x         = pos.annotation[recRow, 2], 
                            y         = 3650, 
                            label     = pos.annotation[recRow, 1], 
                            colour    = "black", 
                            text.size = 8)
}
plot.aic