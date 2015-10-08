# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEInversionSizeDistribution
# -----------------------------------------------------------------------------
# Description:   For Gamma distributions with different scales (x), 2 different
#                ways to calculate the KS-statistic are shown (based on
#                synteny blocks or genes). Furthermore, G-statistic development
#                is shown over different scales. Finally, the optimal
#                gamma distribution is simulated and plotted.
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        5 Plots. First 4 show simulation results for Gamma distri-
#                butions with different scales (x). Scales are drawn from a 
#                uniform distributions; y shows different ways to calculate
#                KS-statistic or alternative measure (G-statistic).
#                Last plot simulates and shows optimal gamma distribution.
# -----------------------------------------------------------------------------
# Keywords:      plot, simulation, uniform, gamma, 
#                Kolmogorov-Smirnov test
# -----------------------------------------------------------------------------
# See also:      -
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEGammaParameters_block.csv,
#                COAPMODEGammaParameters_block_allData.csv,
#                COAPMODEGammaParameters.csv,
#                COAPMODEGammaParameters_allData.csv,
#                COAPMODEGammaParameters_allData.csv
# -----------------------------------------------------------------------------

# Load package
require(ggplot2)

# Load data
# 1) KS based on blocks
df.gammaBlock    = read.table('COAPMODEGammaParameters_block.csv', 
                              sep='\t', 
                              header = T, 
                              stringsAsFactors = F)
df.gammaBlockAll = read.table('COAPMODEGammaParameters_block_allData.csv', 
                              sep='\t', 
                              header = T, 
                              stringsAsFactors = F)
# Change from log to normal scale
df.gammaBlock[,c(4,6,7,8)] = exp(df.gammaBlock[,c(4,6,7,8)])
df.gammaBlockAll[,c(4)]    = exp(df.gammaBlockAll[,c(4)])

# 2) Compare to KS based on genes
df.gamma    = read.table('COAPMODEGammaParameters.csv', 
                         sep='\t', 
                         header = T, 
                         stringsAsFactors = F)
df.gammaAll = read.table('COAPMODEGammaParameters_allData.csv', 
                         sep='\t', 
                         header = T, 
                         stringsAsFactors = F)
# Change from log to normal scale
df.gamma[,c(4,6,7,8)] = exp(df.gamma[,c(4,6,7,8)])
df.gammaAll[,c(4)]    = exp(df.gammaAll[,c(4)])

# 3) Look at g-Measure
df.gamma_gmAll = read.table('COAPMODEGammaParameters_GMeasure_allData.csv', 
                            sep='\t', 
                            header = T, 
                            stringsAsFactors = F)
# Change score from log to normal scale
df.gamma_gmAll[,c(4)] = exp(df.gamma_gmAll[,c(4)])

# Define plotting function

plot.gamma.parameters = function(df, x, xLab, y, yLab, 
                                 plotTitle,
                                 yMin = NA, yMax = NA,
                                 pointsDf = NULL) {
  p = ggplot() 
  if (!is.null(pointsDf)) {
      p = p + 
              geom_point(data    = pointsDf, 
                         mapping = aes_string(x='gammaScale', y='cumAbs'), 
                         size    = 2.0, 
                         alpha   = 0.4, 
                         colour = '#CCCCCC'
                         )
  }
  p = p + 
      geom_line(data    = df, 
                mapping = aes_string(x=x, y=y), 
                size    = 1.3) +
      ggtitle(plotTitle) +
      xlab(xLab) +
      ylab(yLab) +
      theme(axis.line = element_line(colour = "black"), 
            panel.background = element_rect(fill = "white"), 
            panel.grid = element_line(colour = "white"), 
            legend.title=element_blank(), 
            axis.text=element_text(size = 14), 
            axis.title = element_text(size = 14), 
            axis.title.y = element_text(vjust = 1.5),
            axis.title.x = element_text(vjust = 0.2),
            plot.title = element_text(size = 14, vjust=2))
  if (!is.na(yMin)) {
      p = p + 
          geom_line(data    = df, 
                    mapping = aes_string(x='gammaScale', y=yMin), 
                    size    = 1.0, 
                    colour   = '#CCCCCC'
                    )
  }
  if (!is.na(yMax)) {
      p = p + 
        geom_line(data    = df, 
                  mapping = aes_string(x    = 'gammaScale', 
                                       y    = yMax), 
                                       size = 1.0, 
                                       colour = '#CCCCCC')
  }
  return(p)
}

# Kolmogorov-Smnirnoff based on blocks
plot.gamma.parameters(df.gammaBlock, 'gammaScale', 'Gamma(shape = 1, scale = x)', 
                      'cumAbsMean', 'Average KS-statistic', 
                      paste('Average KS-statistic for all MagSimus species ',
                            'comb.\n(gray dots indicate individual ',
                            'simulations)', sep=''), 
                      pointsDf = df.gammaBlockAll)
plot.gamma.parameters(df.gammaBlock, 'gammaScale', 'Gamma(shape = 1, scale = x)', 
                      'HS.MM_mean', 'KS-statistic', 
                      paste('Average KS-statistic for HS-MM\n(gray lines min ',
                            'and max value in simulations)', sep=''), 
                      'HS.MM_min', 'HS.MM_max')

# Kolmogorov-Smnirnoff based on genes
plot.gamma.parameters(df.gamma, 'gammaScale', 'Gamma(shape = 1, scale = x)', 
                      'cumAbsMean', 'Average KS-statistic', 
                      paste('Average KS-statistic for all MagSimus species ',
                            'comb.\nwith each block weighted by the number of',
                            'genes\n(gray dots indicate individual ',
                            'simulations)', sep=''), 
                      pointsDf = df.gammaAll)


# Entropy score based on gMeasure
df.gamma_gm = matrix(rep(NA, 2*length(unique(df.gamma_gmAll$gammaScale))), 
                     ncol=2)
i = 1
for (scale in unique(df.gamma_gmAll$gammaScale)) {
  rows.sel        = df.gamma_gmAll$gammaScale==scale
  cumScale        = prod(df.gamma_gmAll[rows.sel,'cumAbs'])^(1/sum(rows.sel))
  df.gamma_gm[i,] = c(scale, cumScale)
  i = i +1
}
df.gamma_gm = as.data.frame(df.gamma_gm)
colnames(df.gamma_gm) = c('gammaScale', 'cumAbsMean')
plot.gamma.parameters(df.gamma_gm, 'gammaScale', 'Gamma(shape = 1, scale = x)', 
                      'cumAbsMean', 'Geom. mean G-Measure', 
                      paste('Geom. mean of G-Measure for all MagSimus species',
                            ' comb.\n(gray dots indicate individual ',
                            'simulations)', sep=''), 
                      pointsDf = df.gamma_gmAll)

# Plot optimal gamma function
df.gamma.theoretic = data.frame(x = seq(0.01, 100, 0.01), 
                                y = pgamma(seq(0.01, 100, 0.01), 
                                           shape = 1, 
                                           scale = df.gammaBlock$gammaScale[
                                             which.min(
                                               df.gammaBlock$cumAbsMean
                                               )
                                             ]
                                           )
                                )
plot.title = paste0('CDF of selected inversion size distribution:\nGamma(1, ', 
                    round(df.gammaBlock$gammaScale[
                      which.min(df.gammaBlock$cumAbsMean)],
                          2), ')'
                    )
p                  = plot.gamma.parameters(df.gamma.theoretic, 
                                           'x', 
                                           'Inversion size in genes', 
                                           'y', 
                                           '',
                                           plot.title
                                           )
p + ylim(0,1) + xlim(0,100)
