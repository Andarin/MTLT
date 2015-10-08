# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEChrSizesInDifferentScalings
# -----------------------------------------------------------------------------
# Description:   For 21 (or 5) Amniota species, all chromosomes are plotted
#                in two different scalings (x and y-axis). The linear
#                regression model is indicated per species and for the over-all
#                dataset (black dashed). For the latter, the linear model
#                equation with R^2 fit is drawn on the plot.
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        4 Plots showing relationships between x = gene length of a 
#                chromosome and y = chromosome length / non-coding length of
#                the chromosome
# -----------------------------------------------------------------------------
# Keywords:      plot, linear regression, R-squared, estimation, approximation
# -----------------------------------------------------------------------------
# See also:      COAPMODEChrSizesMinMaxAmniota, 
#                COAPMODEChrSizesDistributionsInAmniota
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
#                Some function adapted from Fernando G Taboada
#                https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEChrcGSample.RData
# -----------------------------------------------------------------------------

# Load package
require(ggplot2)

# Load data
# Data contains summary statistics for every chromosome of 21 selected 
# Amniota genomes
load('COAPMODEChrcGSample.RData')

# Define functions
format.species.labels = function(sentence, sep = "_") {
  # Function to reformat labels
  vec.sentence = gsub(sep, " ", sentence, fixed = T)
  out = paste(toupper(substring(vec.sentence, 1, 1)), 
              substring(vec.sentence, 2), sep = "", 
              collapse = "")
  return(out)
}

get.equation = function(model) {
    # Function adapted from https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
    eq.parts = list(a  = format(coef(model)[1], 
                                digits = 2, 
                                big.mark = ',', 
                                scientific=99),
                    b  = format(abs(coef(model)[2]), 
                                digits = 2, 
                                big.mark = ',', 
                                scientific=99),
                    r2 = format(summary(model)$r.squared, 
                                digits = 3)
                    )
  
    if (coef(model)[2] >= 0)  {
        eq = substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
                        eq.parts)
    } else {
        eq = substitute(italic(y) == a - b %.% italic(x)*","~~italic(R)^2~"="~r2,
                        eq.parts)    
    }
  
    as.character(as.expression(eq));                 
}

plot.relations = function(df, x, y, z, 
                          title, annotate.size = 5, 
                          xfac = 0.4, yfac = 10/11) {
    # Function to plot multiple LM and plot equation of regression
    # Function partly adapted from https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
    my.lm = lm(as.formula(paste0(y,'~',x)), df)
    p = ggplot(df, 
               aes_string(x=x, y=y, color=z)
               ) + 
        ggtitle(title) + 
        stat_smooth(method = "lm", se = F) +
        geom_point() +
        geom_abline(intercept = my.lm$coefficients[1], 
                    slope = my.lm$coefficients[2], 
                    linetype="dashed") +
        annotate("text", 
                 x = max(df[x], na.rm = T)*xfac, 
                 y = max(df[y], na.rm = T)*yfac, 
                 label = get.equation(my.lm), 
                 colour="black", 
                 size = annotate.size, 
                 parse=TRUE) +
        theme(axis.line = element_line(colour = "black"), 
              panel.background = element_rect(fill = "white"), 
              panel.grid = element_line(colour = "white"), 
              legend.title=element_blank(), 
              legend.position = "none",
              axis.text=element_text(size = 14), 
              axis.title = element_text(size = 16), 
              title = element_text(size = 12),
              axis.title = element_text(size = 16, vjust=3), 
              title = element_text(size = 12),
              plot.title = element_text(vjust=2, size = 16)) +
        ylab(format.species.labels(y)) +
        xlab(format.species.labels(x))
    return(p)
}

# Plot physical size to gene number data
plot.relations(chr.cG.comp, 'gene_length', 'physical_length', 'species', 
               'Relationship between physical length\nand gene number')
# Now only for the 5 MagSimus species
magsimus.species.underscore = c("homo_sapiens", "mus_musculus", 
                                "canis_familiaris", "monodelphis_domestica", 
                                "gallus_gallus")
data.MS = chr.cG.comp[chr.cG.comp$species %in% magsimus.species.underscore, ]
plot.relations(data.MS, 'gene_length', 'physical_length', 'species', 
    paste('Relationship between physical length and gene number\n(only ',
          'MagSimus species)', sep=''))

# Plot non-coding area to gene number data
plot.relations(chr.cG.comp, 'gene_length', 'noncoding_length', 'species', 
               'Relationship between non-coding area and gene number')
# Now only for the 5 MagSimus species
plot.relations(data.MS, 'gene_length', 'noncoding_length', 'species', 
               'Relationship between non-coding area and gene number')
