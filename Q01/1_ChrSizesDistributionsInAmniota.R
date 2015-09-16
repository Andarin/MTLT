# -----------------------------------------------------------------------------
# Quantlet:      ChrSizesDistributionsInAmniota
# -----------------------------------------------------------------------------
# Description:   Based on the chromosome sizes expressed in relative
#                gene content, the discrete CDFs of 21 real Amniota species
#                (green) and 5 simulated (orange) are displayed.
#                An estimated Amniota start genome is shown, which is the mean
#                of the real Amniota species CDFs.
#                The black dashed line indicates chromosomes of uniform sizes.
#                The blue line shows the limit distribution of
#                proportional translocation sampling
#                The purple line shows the limit distribution of
#                uniform translocation sampling
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        Plot of Discrete Cumulative Density Functions for
#                different real, simulated and estimated Amniota species
# -----------------------------------------------------------------------------
# Keywords:      plot, cdf, sampling, discrete, linear interpolation,
#                simulation
# -----------------------------------------------------------------------------
# See also:      ChrSizesMinMaxAmniota, ChrSizesInDifferentScalings
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      1_df_chromosome_sizes.RData
# -----------------------------------------------------------------------------

require(ggplot2)

# 4 data frames with the following content:
# - real Amniota chromosome sizes
# - mean simulated chromosome sizes for 5 selected Amniota species
# - Limit distributions of uniform and proportional sampling for translocations
# - estimated Amniota chromosome sizes

load('1_df_chromosome_sizes.RData')

ggplot() + 
    geom_line(data = df.dist, 
              aes(x = chrom_cum, 
                  y = gene_cum, 
                  group = species, 
                  colour = species), 
              size=1.3) + 
    geom_point(data = df.amniota.spec, 
               aes(x      = chrom_cum, 
                   y      = gene_cum, 
                   group  = species, 
                   colour = species), 
               show_guide = FALSE) +
    geom_point(data = df.amniota.est, 
               aes(x = chrom_cum, 
                   y = gene_cum), 
               fill   = 'red', 
               colour ='black', 
               shape  = 21, 
               size   = 5) + 
    geom_abline(slope    = 1, 
                size     = 1.5, 
                linetype = "dashed") +
    xlab("Chromosome density") + 
    ylab("Cumulative density of gene content") + 
    ggtitle(paste("Gene content distribution over chromosomes plotted for ",
                  "Amniota species,\nsimulated Amniota genome and limit ",
                  "distributions for different translocation mechanisms", 
                  sep='')
            ) + 
    theme(axis.line        = element_line(colour = "black"), 
          panel.background = element_rect(fill = "white"), 
          panel.grid       = element_line(colour = "white"), 
          legend.title     = element_blank(), 
          legend.position  = "none",
          axis.text        = element_text(size = 14), 
          axis.title       = element_text(size = 16), 
          title            = element_text(size = 12),
          axis.title       = element_text(size = 16,),
          axis.title.y     = element_text(vjust=2),
          title            = element_text(size = 12),
          plot.title       = element_text(vjust=2, size = 16)) + 
    scale_x_continuous(expand = c(0.003, 0)) + 
    scale_y_continuous(expand = c(0.027, 0)) +
    scale_colour_discrete()