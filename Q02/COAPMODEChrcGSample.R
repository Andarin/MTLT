# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEChrcGSample
# -----------------------------------------------------------------------------
# Description:   The minimal and maximal chromosome sizes in different
#                Amniota genomes are shown, in both number of genes as well
#                as number of nucleotides.
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        2 Plots of minimal and maximal chromosome sizes with different
#                y-scale
# -----------------------------------------------------------------------------
# Keywords:      plot, extreme-value, descriptive-statistics, 
#                distribution, bandwidth
# -----------------------------------------------------------------------------
# See also:      COAPMODEChromosomeSizes, ChrSizesInDifferentScalings
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEChrcGSample.RData
# -----------------------------------------------------------------------------


# Load package
require(ggplot2)

# Load data
# Data contains summary statistics for every chromosome of 21 selected 
# Amniota genomes
load("COAPMODEChrcGSample.RData")
magsimus.species = c("homo_sapiens", "mus_musculus", "canis_familiaris", 
                     "monodelphis_domestica", "gallus_gallus")

# Define functions
format.species.labels = function(sentence, sep = "_") {
  # Function to reformat labels
  vec.sentence = gsub(sep, " ", sentence, fixed = T)
  out = paste(toupper(substring(vec.sentence, 1, 1)), 
              substring(vec.sentence, 2), sep = "", 
              collapse = "")
  return(out)
}

create.cap.df = function(df.input, criteria, magsimus.species.underscore) {
  # Function to create min/max chromosome size data frame and graphic
  my.min          = aggregate(chr.cG.comp[,criteria], 
                              list(chr.cG.comp$species), 
                              function(x) return(min(x, na.rm = T)))
  my.min$category = 'min'
  my.max          = aggregate(chr.cG.comp[,criteria], 
                              list(chr.cG.comp$species), 
                              function(x) return(max(x, na.rm = T)))
  my.max$category = 'max'
  out.df                   = rbind(my.min, my.max)
  names(out.df)[1]         = 'species'
  out.df$MagSimus          = 2
  ind.sel                  = out.df$species %in% magsimus.species.underscore
  out.df$MagSimus[ind.sel] = 3
  out.df$MagSimus          = out.df$species %in% magsimus.species.underscore
  
  # Format species names to look nicer
  levels(out.df$species)   = unlist(
    lapply(
      levels(out.df$species), 
      FUN= function(x) format.species.labels(x)
    )
  )
  if (criteria == 'physical_length') {
    out.df$x = out.df$x / 1000000
    lab.y = 'Length in Mb'
  } else {
    lab.y = 'Length in genes'
  }
  p = ggplot(out.df, 
             aes(x=species, y=x, group=category, colour=category, size=MagSimus)
  ) + 
    geom_point() + 
    xlab('') +
    ylab(lab.y) +
    scale_size_manual(name = '', 
                      values = c(3,7), 
                      breaks = c(T,F),
                      labels = c('Selected for\nnum. optimization',
                                 'Not selected for\nnum. optimization')
    )+
    scale_color_manual(name = '', 
                       values = c('tomato1', 'deepskyblue2'), 
                       breaks = c('max','min'),
                       labels = c('biggest chromosome',
                                  'smallest chromosome')
    ) +
    ggtitle('Length of shortest and longest
                    chromosome in different species') + 
    theme(axis.text.x = element_text(angle = 325, 
                                     vjust = 0.92, 
                                     hjust=0.0), 
          axis.line   = element_line(colour = "black"), 
          axis.text   = element_text(size = 11), 
          axis.title  = element_text(size = 16, vjust = 2), 
          title       = element_text(size = 12),
          plot.title  = element_text(vjust=2),
          panel.background   = element_rect(fill = "white"), 
          panel.grid.major.x = element_line(color="lightgray")
    )
  return(list(out.df, p))
}
cap.gene = create.cap.df(chr.cG.comp, "gene_length", magsimus.species)
cap.gene[[2]]
cap.phys = create.cap.df(chr.cG.comp, "physical_length", magsimus.species)
cap.phys[[2]]