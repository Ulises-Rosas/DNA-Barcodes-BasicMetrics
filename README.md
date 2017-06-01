# Basic metrics of DNA barcodes

Strategies of species detection, and even delimitation, using genetic characters are mostly build on the analysis of DNA markers called DNA barcodes. DNA barcoding have demonstrated high resolution in identifying and delimiting species either at local or regional scales if optimal referece library are constructed. Characterization of DNA barcode reference library depends on comprehensive analysis of K2P (i.e. Kimura 2-parameter model) distance matrices. For instance, estimates of **Barcoding Gap** or **Neighbor Species** are inferred from these kind of matrices. Therefore, stringent data exploration of these matrices is a stepping-stone towards characterizing of a DNA barcode reference library.

Several softwares with graphical user-friendly interfaces such as *MEGA* also produce K2P distance matrices. Notwithstanding, its outcomes provide metric by sequences instead of species. Accordingly downstream procedures do not come forward in the same program, but rather traditional pipelines turn to involve more than one program to analyse matrices. This is not only more time consuming, but also it is prone to errors. 

R programing represent an optimal way to sistematically handle large set of sequence information, included DNA barcodes. Upon working in a single console, errors due to file manipulation are eliminated. Currently there are R packages which provide information of barcodes (e.g _Spider_ package). Little is, however, known about direct outcomes of metrics by species (i.e. intraspecific and interspecific variability) using these packages. On top of this, retrieving the Barcoding Gap or Neighbor Species from these ones requiere prior knowledge in R programing and hence troubles emerge for whom do not know R programing.

In the following post, a simple function called `variability` is presented which conducts basic metrics of DNA barcodes by species using base commands and *Ape* package in R. To test that function, DNA sequences mined from the _GenBank_ repository were used. The aim objective is obtain directly foremost data to explore barcodes by species and, in turn, characterize a reference library.

## Input preparation

Before testing the variability function, the input data was prepared. To accomplish this, DNA sequences were download from GenBank. Then, since variability function only carry out its estimates with a single format of sequence names, names of mined sequences were restructured. Finally, only binomial system names were taked by sequence filtering. 

#### Data mining

the _Rentrez_ and _Ape_ packages were used to mine sequences from GenBank. All available sequences of species belonging to the family Sciaenidae were recruited.  

Settings are as the following lines:  

```Rscript
library(rentrez)
library(ape)
mining = entrez_search(db = 'nuccore', 
                       term = "Sciaenidae[Organism] AND (COI[Gene] OR COX[Gene]) AND 
                       (600[SLEN] : 650[SLEN])" ,
                       retmax = 560)
                       
```
```Rscript
library(rentrez)
library(ape)
mining = entrez_search(db = 'nuccore', 
                       term = "Sciaenidae[Organism] AND (COI[Gene] OR COX[Gene]) AND 
                       (600[SLEN] : 650[SLEN])" ,
                       retmax = 560)
                       
```


#### Name structure
```Rscript
seqs= paste(mining$ids, '|', attr(seqs_mining, 'species'), sep = "")
names(seqs_mining) = seqs
```
```
>ID|Name of the specie
TAAGTCAGCCCGGCGCACTTCTCGGAGATGACCAAGTTTA
TAACGTAATTGTTACGGCACATGCCTTCGTTATAATTTTC
TTTA...
```

#### Sequence filtering

```Rscript
sciae.spp = attr(seqs_mining, 'species')
ve = vector('character')
for(i in 1:length(sciae.spp)){
  if(length(strsplit(sciae.spp[i], '_')[[1]]) > 2){
    next
  }else{
    ve[i]=sciae.spp[i]
  }
}
ve2 = ve[!is.na(ve)]
```
```Rscript
new.names = unique(ve2)
new.names.formated = unlist(lapply(new.names, function(x){
  grep(x, names(seqs_mining), value = T)}))

seqs_filtrated = seqs_mining[c(which(names(seqs_mining) %in% new.names.formated))]
write.dna(seqs_filtrated, 'sciaenidae_mined.txt', format = 'fasta',nbcol=1, colw= 60)
```
[sequences](https://github.com/Ulises-Rosas/DNA-Barcodes-BasicMetrics/blob/master/sciaenidae_mined_linsi_gblocks.txt)

## Testing data

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
```
> sciaenidae$Summary_table
                               intra_mean   intra_min   intra_max inter_mean                   neighbor   inter_min          most_divergent inter_max
Pteroscion_peli              0.0056926611 0.005692661 0.005692661  0.2049753           Umbrina_bussingi 0.129972963         Johnius_carouna 0.3704298
Pseudotolithus_senegallus             NaN         Inf        -Inf  0.2004204       Pseudotolithus_typus 0.043292379         Johnius_carouna 0.3266601
Pseudotolithus_typus                  NaN         Inf        -Inf  0.2074564  Pseudotolithus_senegallus 0.043292379     Johnius_heterolepis 0.3214561
Pseudotolithus_brachygnathus 0.0967334894 0.000000000 0.241833723  0.2533595  Pseudotolithus_senegallus 0.112726809     Larimichthys_crocea 0.3652551
Umbrina_cirrosa                       NaN         Inf        -Inf  0.1900687        Umbrina_canariensis 0.075356336         Johnius_carouna 0.3529911
Umbrina_canariensis          0.0041919706 0.000000000 0.009544968  0.1951246            Umbrina_cirrosa 0.075356336      Johnius_borneensis 0.3426139
Sciaena_umbra                0.0018894681 0.001889468 0.001889468  0.1966161            Umbrina_cirrosa 0.088387300         Johnius_carouna 0.3479523
...
```
```
> sciaenidae$Intraspecific_metrics
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.000000 0.001890 0.003788 0.020809 0.015327 0.341685 
> sciaenidae$Interspecific_metrics
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.0000  0.1843  0.2048  0.2225  0.2562  0.4070 
> sciaenidae$Barcoding_Gap
-0.3416851
```
```Rscript
par(mfrow = c(1,2))
hist(sciae$Interspecific_values, col = 'gray', breaks = 150, border = 'gray', 
     xlab = 'K2P distances', main = NULL)
box()
hist(sciae$Intraspecific_values, col = 'green', breaks = 150, border = 'green', add=T)
legend('topleft', 
       legend=c("Intraspecific distances", "Interspecific distances"),
       col=c("green", "gray"), 
       pch = c(15,15), cex=0.95,
       bg='transparent',
       bty = 'n')
plot(sciae$Summary_table[,'intra_max'], sciae$Summary_table[,'inter_min'], 
     xlab = 'Maximun intraspecific distance', ylab = 'Nearest Neighbor distance')
abline(a=0,b =1, col ='gray', lwd=3, lty=5)
text(0.1253226, 0.1150362, ##position retrieved by locator(1) 
     '1:1 Line', col = 'gray', srt=60)
```
![Image of Ulises-Rosas](https://github.com/Ulises-Rosas/DNA-Barcodes-BasicMetrics/blob/master/plots.png)
