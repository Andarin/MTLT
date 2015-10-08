# -----------------------------------------------------------------------------
# Quantlet:      COAPMODECreateStatisticsDatasetForGenomes
# -----------------------------------------------------------------------------
# Description:   Creates dataset for 21 Amniota species, which includes
#                various statistics for their chromosomes, partially downloaded
#                from the Ensembl Project. The gene name data here is only a
#                proof of concept: small sample of the actually gene name data 
#                used which could not be uploaded due to data restrictions.
#                The descriptive data is visualized with the quantlets
#                indicated in "See also".
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        1 Plot showing the phylogenetic tree of 21 Amniota
#                A dataset with multiple chromosome statistics 
# -----------------------------------------------------------------------------
# Keywords:      data mining, preprocessing, descriptive-statistics, tree, 
#                discrete
# -----------------------------------------------------------------------------
# See also:      COAPMODEMultipleLinearRegression, COAPMODEChrcGSample
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODECleanGenomesSample.RData, COAPMODEChrPhysicalSize.csv, 
#                COAPMODEAmniotaGood.nwk
# -----------------------------------------------------------------------------


# Load packages
require(biomaRt)
require(reshape2)
require(ape)

# Load data
load('COAPMODECleanGenomesSample.RData')
chr.physical_size = read.csv('COAPMODEChrPhysicalSize.csv')
tree.amniota.good = read.tree('COAPMODEAmniotaGood.nwk')
magsimus.species.underscore = c('homo_sapiens', 'mus_musculus', 
                                'canis_familiaris', 'monodelphis_domestica', 
                                'gallus_gallus')


# Download Ensembl data
# listMarts()
# ensembl=useMart("ensembl")
# listDatasets(ensembl)

ensembl = useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
species = as.character(listDatasets(ensembl)$dataset)
all     = data.frame('species'=NA, 'chromosome_name'=NA, 
                     'ensembl_gene_id'=NA, 'start_position'=NA, 
                     'end_position'=NA)
names(all) = c('species','chromosome_name', 'ensembl_gene_id', 
               'start_position', 'end_position')
for (i in seq_along(species)) {
  ensembl_rec = useDataset(species[i],mart=ensembl)
  rec         = getBM(attributes = c('chromosome_name', 'ensembl_gene_id', 
                                     'start_position', 'end_position'), 
                      mart=ensembl_rec)
  rec$species = species[i]
  all         = rbind(all, rec)
}

all = all[-1,]
row.names(chr.physical_size) = chr.physical_size$X
all$species                  = substr(all$species, 1, nchar(all$species) - 13)
# felis_catus has some doubling in the data base, hence delete them
all        = unique(all)
names.all  = names(table(all$species))
names.real = colnames(chr.physical_size)
names.species = unlist(lapply(strsplit(names.real, '_'), 
                              function(x) return (x[length(x)]))
                       )
for (names.all.ind in 1:length(names.all)) {
  ind.rec       = which(substr(names.all[names.all.ind], 
                               2, 
                               nchar(names.all[names.all.ind]))==names.species)
  real.name.rec = names.real[ind.rec]
  all$species[all$species == names.all[names.all.ind]] = real.name.rec
}
species.good       = tree.amniota.good$tip.label
species.names.corr = tolower(species.good)
species.names.corr[species.names.corr=="gorilla.gorilla.gorilla"] = 
    "gorilla.gorilla"
species.names.corr[species.names.corr=="canis.lupus.familiaris"] = 
    "canis.familiaris"

all=all[gsub('_', '.', all$species) %in% species.names.corr,]

plot(tree.amniota.good, show.node.label = F, 
     edge.width = 0.2, cex = 0.6, label.offset = 3)
nodelabels(text = tree.amniota.good$node.label, 
           bg = 'white', font = 1, cex = 0.5, adj = 0.5)

all.cG.good        = all.cleanGenomes.sample[names(all.cleanGenomes.sample) %in% 
                                          species.good]
all.cG.good        = all.cG.good[order(names(all.cG.good))]
names(all.cG.good) = levels(as.factor(all$species))
for (id.genome in seq_along(all.cG.good)) {
    rec.name   = names(all.cG.good)[id.genome]
    rec.genome = cbind(rep(rec.name, 
                           nrow(all.cG.good[[id.genome]])),
                       all.cG.good[[id.genome]])
    if (id.genome == 1) {
        all.cG.good.long = rec.genome
    } else {
        all.cG.good.long = rbind(all.cG.good.long, rec.genome)
    }
}

all.cG.good.long = unique(all.cG.good.long)
all.cG.good.long = all.cG.good.long[,-ncol(all.cG.good.long)]
names(all.cG.good.long) = c(names(all)[c(1,2,4,5)],
                            'orientation',
                            'ensembl_gene_id')

all.cG.good.long$species = as.character(all.cG.good.long$species)
all.cG                   = all.cG.good.long

# all.cG has shortest transcription gene positions
names(all.cG)[grep('position', names(all.cG))] = 
    paste(names(all.cG)[grep('position', names(all.cG))], '_shortest', sep='')

all.cG$start_position_shortest = all.cG$start_position
all.cG$end_position_shortest   = all.cG$end_position
all.allG = all
all.cG   = merge(all.cG, 
               all.allG[, c('ensembl_gene_id', 
                            'start_position', 
                            'end_position')],
               by = 'ensembl_gene_id', all.x = T)

create.chr.data = function(all, chr.physical_size) {
    chr.physical_size_sel = chr.physical_size[colnames(chr.physical_size) %in% 
                                                  names(table(all$species))]
    # filter NA-rows
    chr.physical_size_sel = 
        chr.physical_size_sel[!apply(is.na(chr.physical_size_sel), 
                                     1, 
                                     function(x) all(x)),]
    for (rname in row.names(chr.physical_size_sel)) {
        dupl = (toupper(row.names(chr.physical_size_sel)) == toupper(rname))
        if (sum(dupl) < 2) {
            next
        }
        for (i in seq_len(ncol(chr.physical_size_sel))) {
            if (!(is.na(chr.physical_size_sel[tolower(rname),i])) 
                & !(is.na(chr.physical_size_sel[toupper(rname),i]))) {
                stop('Cannot resolve merging conflict')
            }
            if (!(is.na(chr.physical_size_sel[tolower(rname),i]))) {
                rec.name = tolower(rname)
            } else {
                rec.name = toupper(rname)
            }
            print(paste(i, rname, rec.name))
            chr.physical_size_sel[toupper(rname),i] = 
                chr.physical_size_sel[rec.name,i]
      }
      chr.physical_size_sel = 
          chr.physical_size_sel[row.names(chr.physical_size_sel) 
                                != tolower(rname),]
    }
    
    all['species_f'] = as.factor(all$species)
    all['gene_length'] = all$end_position - all$start_position
    chr.gene_size = matrix(rep(NA, length(levels(all$species_f))), 
                           nrow = 1, 
                           dimnames = list(c('blubb'),levels(all$species_f)))
    for (specie in levels(all$species_f)) {
        # chromosome size by gene number
        t.rec = table(all[all$species_f==specie,'chromosome_name'])
        t.rec = t.rec[sapply(names(t.rec), nchar)<=4]
        for (rec_chr in names(t.rec)) {
            rec_chr_upper = toupper(rec_chr)
            if (!(rec_chr_upper %in% row.names(chr.gene_size))) {
                print(paste(specie, rec_chr_upper))
                chr.gene_size = rbind(chr.gene_size, NA)
                row.names(chr.gene_size)[nrow(chr.gene_size)] = rec_chr_upper
            }
            chr.gene_size[rec_chr_upper, specie] = t.rec[rec_chr]
        }
    }
    chr.gene_size = chr.gene_size[-1,]
    chr.gene_size = chr.gene_size[row.names(chr.gene_size)!='MT',]
    chr.gene_size = data.frame(chr.gene_size)
    # chromosome size by coding/non-coding base pairs
    all.selection = all[nchar(all$chromosome_name)<=4,]
    chr.l.coding_size = aggregate(gene_length ~ species + chromosome_name, 
                                  data = all.selection, FUN = 'sum')
    names(chr.l.coding_size)[3] = 'coding_length'
    chr.l.coding_size$chromosome_name = 
        toupper(chr.l.coding_size$chromosome_name)
    row.names(chr.physical_size_sel) = 
        toupper(row.names(chr.physical_size_sel))
    row.names(chr.gene_size) = toupper(row.names(chr.gene_size))
    chr.l.gene_size = melt(as.matrix(chr.gene_size))
    names(chr.l.gene_size) = c(names(chr.l.coding_size)[c(2,1)], 'gene_length')
    chr.l.physical_size = melt(as.matrix(chr.physical_size_sel))
    names(chr.l.physical_size) = c(names(chr.l.coding_size)[c(2,1)], 
                                   'physical_length')
    
    chr.data = merge(chr.l.gene_size, 
                     chr.l.physical_size, 
                     by=c('species', 'chromosome_name'), all.x=T, all.y=F)
    chr.data = merge(chr.data, chr.l.coding_size, 
                     by=c('species', 'chromosome_name'), all.x=T, all.y=F)
    chr.data$coding.to.total = chr.data$coding_length/chr.data$physical_length
    chr.data = chr.data[!apply(is.na(chr.data), 1, function(x) all(x)),]
    first.nine.chr = nchar(as.integer(levels(chr.data$chromosome_name)))==1
    levels(chr.data$chromosome_name)[first.nine.chr] = 
        paste0('0', levels(chr.data$chromosome_name)[first.nine.chr])
    chr.data$chromosome_name = factor(chr.data$chromosome_name,
                                      levels(chr.data$chromosome_name)[order(
                                          levels(chr.data$chromosome_name))], 
                                      ordered = T)
    chr.data$noncoding_length = chr.data$physical_length-chr.data$coding_length
    chr.data$chr_gene_number = chr.data$chr_noncoding_number = 
        chr.data$chr_physical_number = NA
    chr.data$gene_density = chr.data$noncoding_density = 
        chr.data$physical_density = NA
    for (specie in levels(chr.data$species)) {
        rec_rows = which(chr.data$species==specie)
        # chr gene length number
        rec_gene_number = order(chr.data[rec_rows,'gene_length'], 
                                decreasing = T, na.last = NA)
        chr.data[rec_rows,'chr_gene_number'][rec_gene_number] = 
            1:length(rec_gene_number)
        chr.data$gene_density[rec_rows] = 
            chr.data[rec_rows,'gene_length']/sum(chr.data[rec_rows,
                                                          'gene_length'], 
                                                 na.rm=T)
        # chr physical length number
        rec_physical_number = order(chr.data[rec_rows,'physical_length'], 
                                    decreasing = T, na.last = NA)
        chr.data[rec_rows,'chr_physical_number'][rec_physical_number] = 
            1:length(rec_physical_number)
        chr.data$physical_density[rec_rows] = 
            chr.data[rec_rows,'physical_length']/sum(chr.data[
                rec_rows,'physical_length'], na.rm=T)
        # chr noncoding length number
        rec_noncoding_number = order(chr.data[rec_rows,'noncoding_length'], 
                                     decreasing = T, na.last = NA)
        chr.data[rec_rows,'chr_noncoding_number'][rec_noncoding_number] = 
            1:length(rec_noncoding_number)
        chr.data$noncoding_density[rec_rows] = chr.data[
            rec_rows,'noncoding_length']/sum(chr.data[rec_rows,
                                                      'noncoding_length'], 
                                             na.rm=T)
    }
    
    cumulate = function (df, order.var1, order.var2, cum.var) {
        df=df[order(df[,order.var1], df[,order.var2]),]
        df$ycum = NA
        for (line in seq_len(nrow(df))) {
            if (line == 1) {
              df[line,'ycum'] = df[line,cum.var]
            } else if (df[line-1,order.var1] != df[line,order.var1] ) {
              df[line,'ycum'] = df[line,cum.var]
            } else {
              df[line,'ycum'] = df[line,cum.var] + df[line-1,'ycum']
            }
        }
        df = df[,names(df)!=cum.var]
        return(df)
    }
    
    small.df = function(df, vars) {
        df.short = df[,vars]
        cdsa3 = aggregate(df.short[,vars[3]], 
                          by=list(df.short[,vars[1]]), 
                          function(x) return(sum(x, na.rm = T)))
        cdsa2 = aggregate(df.short[,vars[2]], 
                          by=list(df.short[,vars[1]]), 
                          function(x) return(max(x, na.rm = T)))
        for (species in levels(df.short$species)) {
          df.short[df.short$species==species,vars[2]] = 
              df.short[df.short$species==species,vars[2]]/
              cdsa2[cdsa2[,1]==species,2]
          df.short[df.short$species==species,vars[3]] = 
              df.short[df.short$species==species,vars[3]]/
              cdsa3[cdsa3[,1]==species,2]
        }
        names(df.short) = c("species","x","y")
        df.short = cumulate(df.short, 'species', 'x', 'y')
        return(df.short)
    }
    chr_gene_number.rel = c()
    for (species in levels(chr.data$species)) {
        temp = chr.data[chr.data$species==species,'chr_gene_number']
        chr_gene_number.rel = c(chr_gene_number.rel, temp/max(temp, na.rm = T))
    }
    chr.data$chr_gene_number.rel = chr_gene_number.rel
    
    chr.data=chr.data[order(chr.data$species, chr.data$chr_gene_number.rel),]
    chr.data.short = small.df(chr.data, 
                              c("species","chr_gene_number","gene_length"))
  #   chr.data.short.50 = small.df(chr.data[which(chr.data$gene_length>=50),], 
  #                                c("species","chr_gene_number","gene_length"))
    gene_length.on.non_coding = data.frame()
    for (specie in levels(chr.data$species)) {
        tmp = chr.data[chr.data$species==specie,]
        rec.lm = lm(formula = noncoding_length ~ gene_length, data = tmp)
        rec.sum = summary(rec.lm)    
        if (length(gene_length.on.non_coding) == 0) {
            gene_length.on.non_coding = 
                rbind(c(specie,(rec.sum$coefficients)['gene_length',
                                                    c('Estimate','Pr(>|t|)')], 
                      rec.sum$adj.r.squared))
        } else {
          gene_length.on.non_coding = 
                rbind(gene_length.on.non_coding, rbind(c(specie,
                                        (rec.sum$coefficients)['gene_length', 
                                        c('Estimate','Pr(>|t|)')], 
                                                      rec.sum$adj.r.squared)))
          }
    }
    gene_length.on.non_coding = as.data.frame(gene_length.on.non_coding)
    names(gene_length.on.non_coding) = c('species', 'slope', 'pvalue',
                                         'adj_R_squared')
    
    
    return(list(chr.data=chr.data, 
                chr.data.short=chr.data.short, 
  #               chr.data.short.50=chr.data.short.50, 
                gene_length.on.non_coding=gene_length.on.non_coding))
}

# If still included, exclude Y and W chromosome
all.cG.woYW            = all.cG[all.cG$chromosome_name!='Y' & 
                                all.cG$chromosome_name!='W',]
chr.physical_size.woYW = chr.physical_size[row.names(chr.physical_size)!='Y' & 
                                           row.names(chr.physical_size)!='W',]

# Only coding genes (~20,000 for human)
chr.data.cG = create.chr.data(all.cG.woYW, chr.physical_size.woYW)

chr.cG.comp = 
    chr.data.cG$chr.data[!(is.na(chr.data.cG$chr.data$gene_length) | 
                           is.na(chr.data.cG$chr.data$physical_length) | 
                           is.na(chr.data.cG$chr.data$coding_length)), ]
save(chr.cG.comp, file = "2_chr_cG_sample.RData")
