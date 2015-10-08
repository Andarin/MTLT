# -----------------------------------------------------------------------------
# Quantlet:      COAPMODEEntropyScore
# -----------------------------------------------------------------------------
# Description:   Based on the contingency tables for genome-genome comparisons
#                for different gene counting methodologies, this script
#                creates a databse with the log-p-values for chi-square test
#                and g-test (Kullback-Leibler)
# -----------------------------------------------------------------------------
# Inputs:        -
# -----------------------------------------------------------------------------
# Output:        Dataset for plotting in 
#                EntropyScoreComparisonsForGenomeAnalysis
# -----------------------------------------------------------------------------
# Keywords:      contingency table, discrete, distance, chi-square, 
#                Kullback-Leibler
# -----------------------------------------------------------------------------
# See also:      COAPMODEEntropyScore
# -----------------------------------------------------------------------------
# Author:        Lucas Tittmann, 2015-08-14
# -----------------------------------------------------------------------------
# Datafile:      COAPMODEAmniotaGood.nwk, COAPMODEHomologyTables/*.csv
# -----------------------------------------------------------------------------

# Load packages
require(reshape2)
require(ape)

# Load data
# At first, load phylogenetic tree
tree.amniota.good = read.tree("COAPMODEAmniotaGood.nwk")
amniota.good.dist = cophenetic.phylo(tree.amniota.good)
# Secondly, load contingency tables for genome-genome comparisons
setwd("COAPMODEHomologyTables/")
files = list.files(".")

# Define function
calc.inequality = function(disp.list, species) {
    disp.list.long = list()
    for (i in seq_along(disp.list)) {
        el = disp.list[[i]]
        tmp = melt(as.matrix(el))
        names(tmp) = c(species[2:1], names(disp.list)[i])
        if (i == 1) {
            disp.list.long = tmp
        } else {
            disp.list.long = merge(x = disp.list.long, 
                                   y = tmp, 
                                   by = names(tmp)[1:2])
        }
        
        for (j in 1:2) {
            HSagg = aggregate(x = tmp[[3]], 
                              by = list(tmp[[j]]), 
                              FUN = sum)
            HSagg = cbind(HSagg, 
                          sum(HSagg$x), 
                          HSagg[[2]]/sum(HSagg$x))
            names(HSagg) = c(names(tmp)[j], 
                             paste0("sum_", 
                                    names(tmp)[j], "_", names(disp.list)[i]), 
                             paste0("total_", names(disp.list)[i]), 
                             paste0("marginal_", names(tmp)[j], 
                  "_", names(disp.list)[i]))
            if (j == 2) {
                HSagg = HSagg[, -3]
            }
            disp.list.long = merge(x = disp.list.long, 
                                   y = HSagg, 
                                   by = names(tmp)[j])
        }
        disp.list.long$a = disp.list.long[, 
            paste0("marginal_", names(tmp)[1], "_", names(disp.list)[i])] * 
            disp.list.long[, paste0("marginal_", names(tmp)[2], 
                                    "_", names(disp.list)[i])]
        names(disp.list.long)[length(names(disp.list.long))] = 
            paste0("exp_prop_", names(disp.list)[i])
        
        disp.list.long$a = disp.list.long[, paste0("sum_", names(tmp)[1], "_", names(disp.list)[i])] * 
            disp.list.long[, paste0("sum_", names(tmp)[2], 
                                    "_", names(disp.list)[i])] * 
            1/disp.list.long[, paste0("total_", names(disp.list)[i])]
        names(disp.list.long)[length(names(disp.list.long))] = 
            paste0("exp_", names(disp.list)[i])
        
        disp.list.long$a = disp.list.long[, names(disp.list)[i]] * 
            1/disp.list.long[, paste0("total_", names(disp.list)[i])]
        names(disp.list.long)[length(names(disp.list.long))] = 
            paste0("real_prop_", names(disp.list)[i])
        
        tmp.exp = disp.list.long[, paste0("exp_", names(disp.list)[i])]
        tmp.real = disp.list.long[, names(disp.list)[i]]
        disp.list.long$a = (tmp.exp - tmp.real)^2 * 1/tmp.exp
        names(disp.list.long)[length(names(disp.list.long))] = 
          paste0("chi2_", names(disp.list)[i])
        
        disp.list.long$a = 2 * tmp.real * (log(tmp.real) - log(tmp.exp))
        names(disp.list.long)[length(names(disp.list.long))] = 
            paste0("G_", names(disp.list)[i])
    }
    disp.list.long = disp.list.long[, c(order(names(disp.list.long)[1:2]), 
                                        3:ncol(disp.list.long))]
    return(disp.list.long)
}

file.ind = 1
ineq.list = list()
while (file.ind <= length(files)) {
    disp.list = list()
    species = c(strsplit(files[file.ind], "-")[[1]][1], 
                strsplit(files[file.ind], 
        "-")[[1]][2])
    disp.list$alldup = read.table(files[file.ind], 
                                  header = T, 
                                  sep = ",", 
                                  check.names = F)
    disp.list$nodup = read.table(files[file.ind + 1], 
                                 header = T, 
                                 sep = ",", 
                                 check.names = F)
    disp.list$prob = read.table(files[file.ind + 2], 
                                header = T, 
                                sep = ",", 
                                check.names = F)
    
    ineq = calc.inequality(disp.list, species)
    ineq.s = ineq[, c(1, 2, 
                      grep("chi2_", names(ineq)), 
                      grep("G_", names(ineq)))]
    pseudo.cnt = 1
    ineq.pseudo = calc.inequality(lapply(disp.list, 
                                         FUN = function(x) 
                                           return(x + pseudo.cnt)), species)
    ineq.pseudo.s = ineq.pseudo[, c(1, 2, grep("chi2_", names(ineq.pseudo)), 
                                    grep("G_", names(ineq.pseudo)))]
    ineq.list = append(ineq.list, list(list(ineq.s, ineq.pseudo.s)))
    # Pseudocount makes the p-value a bit more conservative (bending to a 
    # uniform distribution)
    statistics = colSums(ineq.s[, -c(1:2)], na.rm = T)
    statistics.p = colSums(ineq.pseudo.s[, -c(1:2)], na.rm = T)
    dfree = (nrow(disp.list[[1]]) - 1) * (ncol(disp.list[[1]]) - 1)
    p.values = as.data.frame(matrix(pchisq(statistics, 
                                           dfree, 
                                           log.p = T, 
                                           lower.tail = F), 
                                           nrow = 1))
    p.values.p = as.data.frame(matrix(pchisq(statistics.p, 
                                             dfree, 
                                             log.p = T, 
                                             lower.tail = F), 
        nrow = 1))
    names(p.values.p) = names(p.values) = names(statistics)
    p.values$species = p.values.p$species = paste(species, collapse = "-")
    if (file.ind == 1) {
        p.values.all = p.values
        p.values.p.all = p.values.p
    } else {
        p.values.all = rbind(p.values.all, p.values)
        p.values.p.all = rbind(p.values.p.all, p.values.p)
    }
    file.ind = file.ind + 3
}

phyl.dist.long = melt(amniota.good.dist)
phyl.dist.long$species = gsub(".", "_", paste(phyl.dist.long[, 1], 
                                              phyl.dist.long[, 2], 
                                              sep = "-"), fixed = T)
phyl.dist.long$Var1 = NULL
phyl.dist.long$Var2 = NULL

p.values.all = merge(p.values.all, phyl.dist.long, by = "species")
names(p.values.all)[length(names(p.values.all))] = "phylo_distance"
p.values.all.long = melt(p.values.all, id.vars = c("species", 
                                                   "phylo_distance"))

p.values.p.all = merge(p.values.p.all, phyl.dist.long, by = "species")
names(p.values.p.all)[length(names(p.values.p.all))] = "phylo_distance"
p.values.p.all.long = melt(p.values.p.all, 
                           id.vars = c("species", "phylo_distance"))

p.values.all.all = merge(p.values.all.long, 
                         p.values.p.all.long[, names(p.values.p.all.long) != 
                   "phylo_distance"], by = c("species", "variable"))
names(p.values.all.all)[(ncol(p.values.all.all) - 1):
                          ncol(p.values.all.all)] = c("normal", "pseudocount")

df.entropy.score = p.values.all.all[grep("_nodup", p.values.all.all$variable),]
save(df.entropy.score, file = "4_df_entropy_score.RData") 
