# Basic metrics of DNA barcodes

Strategies of species detection, and even delimitation, using genetic characters are mostly build on the analysis of DNA markers called DNA barcodes. DNA barcoding have demonstrated Resolution of identification and delmitation both at local and regional if optimal referece libray are build. Notwithstanding, characterization of DNA barcode reference library depends on comprehensive analysis of K2P (i.e. Kimura 2-parameter model) distance matrices. Estimates of **Barcoding Gap** or **Neighbor species** are inferred from these kind of matrices. Therefore, data exploration of these matrices is a stepping tone process

Several Programs with graphical interfaz such as *MEGA* produce distance matrices as outcomes. Then, downstream procedures, however, do not rely   ,but rather outside of these programs interfaz are traditionally involved in the pipeline. Prone to errors

Neighbor or barcoding gap, it is not handled and it should invert too much code in R. Exploratory data. Data requiere to plot information. Little is know about tat program

R packages can represent a good represent a good option to handled sistematically large set of sequence information. Currently there are packages which provide information of barcodes (e.g _Spider_ package), however its demand you prior knowlegde in command R. The user experience will depend

, variability function is presented using base function and *Ape* package in R, using as input only fasta equences with formated names. The aim objective is obtain data quicklier to explore your barcodes and in turn, characterize your refernce library

## Installation and 

## Sequence name



Variability measures  of DNA barcodes sequences (i.e. Intraspecific and Interespecific) by using commands in R

```Rscript
variability <- function(barcodes){
        ##species name obtained from the self structure of sequence name
        spp = sapply(strsplit(names(barcodes), "\\|"), function(a){a[2]})
        ##species name
        uniq_spp = unique(spp)
##List storing intraspecific metrics by specie
intra_list = lapply(uniq_spp, function(b){
        g.tmp = grep(b, names(barcodes), value = T)
        intra_by_spp = barcodes[c(which(names(barcodes) %in% g.tmp))]
        dist.dna(intra_by_spp, "k80")})
##creating intraspecific data frame
dat = data.frame(intra_mean = sapply(intra_list, mean), 
                 intra_min = sapply(intra_list, min), 
                 intra_max = sapply(intra_list, max))   
row.names(dat) = uniq_spp
##as.matrix to boost the manipulation of row names
c = dist.dna(barcodes, model = "k80", as.matrix = T) ##whole distance matrix
d = list()
##with names to batch with them and to have as rows as length of unique(spp)
for(i in 1:length(spp)){
        d[[rownames(c)[i]]] = c[i,][!(rownames(c) %in% grep(spp[i], rownames(c), value = T))]
        }
##List storing interspecific metrics by specie
inter_list = lapply(uniq_spp, function(e){
        as.data.frame(d[names(d) %in% grep(e, names(d), value = T)])})

neighbor_each_column = sapply(sapply(inter_list, as.matrix), function(f){
        sapply(strsplit(rownames(which(f == min(f), arr.ind = T)), "\\|"), function(g){g[2]})
        })
neighbor_by_spp = sapply(neighbor_each_column, function(h){
        if(length(unique(h)) == 1){
                unique(h)
        }else{
                "Shared"}})
most_divergent_each_column = sapply(sapply(inter_list, as.matrix), function(j){
        sapply(strsplit(rownames(which(j == max(j), arr.ind = T)), "\\|"), function(k){k[2]})
        })
most_divergent_by_spp = sapply(most_divergent_each_column, function(l){
        if(length(unique(l)) == 1){
                unique(l)
        }else{
                "Shared"}})
##creating interspecific data frame
dat2 = data.frame(inter_mean = sapply(sapply(inter_list, as.matrix), mean),
           neighbor = neighbor_by_spp,
           inter_min = sapply(sapply(inter_list, as.matrix), min),
           most_divergent = most_divergent_by_spp,
           inter_max = sapply(sapply(inter_list, as.matrix), max))
row.names(dat2) = uniq_spp

main_list = list(Summary_table = cbind(dat, dat2),
Intraspecific_metrics = summary(unlist(intra_list)),
Interspecific_metrics = summary(unlist(sapply(d, as.matrix))),
Intraspecific_values = unlist(intra_list),
Interspecific_values = unlist(sapply(unname(d), as.numeric)),
Barcoding_Gap = min(unlist(sapply(unname(d), as.numeric))) - max(unlist(intra_list)))
return(main_list)
}
```
