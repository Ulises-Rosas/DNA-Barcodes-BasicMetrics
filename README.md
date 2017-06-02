# Basic metrics of DNA barcodes

Strategies of species detection, and even delimitation, using genetic characters are mostly build on the analysis of DNA markers called DNA barcodes. DNA barcoding have demonstrated high resolution in identifying and delimiting species either at local or regional scales if optimal referece library are constructed. Characterization of DNA barcode reference library depends on comprehensive analysis of K2P (i.e. Kimura 2-parameter model) distance matrices. For instance, estimates of **Barcoding Gap** or **Neighbor Species** are inferred from these kind of matrices. Therefore, stringent data exploration of these matrices is a stepping-stone towards characterizing of a DNA barcode reference library.

Several softwares with graphical user-friendly interfaces such as *MEGA* also produce K2P distance matrices. Notwithstanding, its outcomes provide metrics by sequences instead of species. Accordingly downstream procedures do not come forward in the same program, but rather traditional pipelines turn to involve more than one program to analyse matrices. This is not only more time consuming, but also it is prone to errors. 

R programing represent an optimal way to sistematically handle large set of sequence information, included DNA barcodes. Upon working in a single console, errors due to file manipulation are eliminated. Currently there are R packages which provide information of barcodes (e.g _Spider_ package). Little is, however, known about direct outcomes of metrics by species (i.e. intraspecific and interspecific variability) using these packages. On top of this, retrieving the Barcoding Gap or Neighbor Species from these ones requiere prior knowledge in R programming and hence troubles emerge for whom do not know R programming.

In the following post, a simple function called `variability` is presented which conducts basic metrics of DNA barcodes by species using base commands and *Ape* package in R. To test that function, mined DNA sequences from the _GenBank_ repository were used. The aim objective is obtain directly foremost data to explore barcodes by species and, in turn, characterize a reference library.

## Input preparation
Before testing the variability function, input data was prepared. To accomplish this, DNA sequences were download from GenBank. Then, since variability function only carry out its estimates with a single format of sequence names, names of mined sequences were restructured. Finally, only binomial system names were taken by sequence filtering. 


#### Data mining
On this case, _Rentrez_ and _Ape_ packages were used to mine barcodes of species belonging to the family Sciaenidae from GenBank. All available ID's of COI or COX sequences, whose length is between 600-650 bp, were recruited. That interval of sequence lengths is due that barcode region have roughly 650 bp of length. Codes are shown in the following lines: 
```Rscript
library(rentrez)
library(ape)
ID_sciaenidae = entrez_search(db = 'nuccore', term = "Sciaenidae[Organism] AND (COI[Gene] OR COX[Gene]) AND 
                      (600[SLEN] : 650[SLEN])" ,
                       retmax = 560) ## only 560 were obtained following above conditions
                       
```
Thus, function _read.GenBank()_ download all those sequences whose IDs, stored in `ID_scaenidae` object, match with IDs within the GenBank repository:

```Rscript
seqs_sciaenidae = read.GenBank(ID_sciaenidae$ids)
```


#### Name structure
The variability function run with a single format of sequence names and it is like the following:
```
>ID|Name of the specie
TAAGTCAGCCCGGCGCACTTCTCGGAGATGACCAAGTTTA
TAACGTAATTGTTACGGCACATGCCTTCGTTATAATTTTC
TTTA...
```
In consequence, once obtained raw sequences, we need to change its name format. The forthcoming code is used to change names of `seqs_sciaenidae` object:
```Rscript
formatted_names = paste(ID_sciaenidae$ids, '|', attr(seqs_sciaenidae, 'species'), sep = "")
names(seqs_sciaenidae) = formatted_names
```


#### Sequence filtering

DNA barcoding detects signs of incomplete lineage sorting if correct identification of species by using morphological characters is conducted, its validation, however, needs additional markers. Thus, considering we are characterizing a DNA barcode reference library, in this example we will constraint us to analyse only sequences of species identified at species-level:

```Rscript
whole.spp = attr(seqs_sciaenidae, 'species')
binomial.spp = vector('character')
for(i in 1:length(whole.spp)){
        if(length(strsplit(whole.spp[i], '_')[[1]]) > 2){ 
                next  ##eliminate all names whose structure are composed by more than two words
        }else if(strsplit(whole.spp[i], '_')[[1]][2] == "sp"){
                next  ##eliminate all names composed with "sp" (i.e. just at genus-level identification)
        }else {
                binomial.spp[i]=whole.spp[i] ##recover binomal names only
        }
}
binomial.spp.whitout.na = binomial.spp[!is.na(binomial.spp)] ##erase all "NA" wrote by next() function
```
Then, unique names were held them as argument pattern in a _grep()_ function. It searches for matches within each element of the `names(seqs_sciaenidae)` vector. Matches were, in turn, used to extract sequences:

```Rscript
filtered.names = unique(binomial.spp.whitout.na) ##unique names of binomial.spp.whitout.na object

filtered.names.formated = unlist(lapply(filtered.names, function(x){
  grep(x, names(seqs_sciaenidae), value = T)})) ##it searches for matches  

##sequence extraction using names 
seqs_filtrated = seqs_sciaenidae[c(which(names(seqs_sciaenidae) %in% filtered.names.formated))]
```
Above command prepare the final version of sequences which will stay in our analyzes. On that account, sequences were exported in text file to align them:

```Rscript
write.dna(seqs_filtrated, 'sciaenidae_mined.txt', format = 'fasta',nbcol=1, colw= 60)
```
Alignments were performend with the program **MAFFT-L-INS-i** (i.e. Smith-Waterman algorithm). Those ends outside the alignment zone were cutted with the program **Gblocks**. Finally, we already have those sequences which will act as input data to variability function: [sequences](https://github.com/Ulises-Rosas/DNA-Barcodes-BasicMetrics/blob/master/sciaenidae_mined_linsi_gblocks.txt)

## Testing data
To test variability function we firstly must read our sequences :
```Rscript
>sciaenidae.barcodes = read.FASTA("sciaenidae_mined_linsi_gblocks.txt")
```
Then, we run the `variabilty` script:
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
calculate basic metric of DNA barcodes:

```
>metrics.sciaenidae = variability(sciaenidae.barcodes)
```
If it appears a warning message, it is due there are species represented by a single sequence. Summary table shows basic metrics by species of sequences. It include, **Neighbor Species** and **Most Divergent Species**:

```
> metrics.sciaenidae$Summary_table
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
Also, as part of their outcomes, variability function estimates the **Barcoding Gap**:
```
> metrics.sciaenidae$Intraspecific_metrics
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.000000 0.001890 0.003788 0.020809 0.015327 0.341685 
> metrics.sciaenidae$Interspecific_metrics
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.0000  0.1843  0.2048  0.2225  0.2562  0.4070 
> metrics.sciaenidae$Barcoding_Gap
-0.3416851
```
If its value is negative, it means that there are not a **Barcoding Gap** using the  whole reference library. We can also graphically demonstrate it running the following code:
```Rscript
par(mfrow = c(1,2))
hist(metrics.sciaenidae$Interspecific_values, col = 'gray', breaks = 150, border = 'gray', 
     xlab = 'K2P distances', main = NULL)
box()
hist(metrics.sciaenidae$Intraspecific_values, col = 'green', breaks = 150, border = 'green', add=T)
legend('topleft', 
       legend=c("Intraspecific distances", "Interspecific distances"),
       col=c("green", "gray"), 
       pch = c(15,15), cex=0.95,
       bg='transparent',
       bty = 'n')
plot(metrics.sciaenidae$Summary_table[,'intra_max'], metrics.sciaenidae$Summary_table[,'inter_min'], 
     xlab = 'Maximun intraspecific distance', ylab = 'Nearest Neighbor distance')
abline(a=0,b =1, col ='gray', lwd=3, lty=5)
text(0.1253226, 0.1150362, ##position retrieved by locator(1) 
     '1:1 Line', col = 'gray', srt=60)
```
These code used `metrics.sciaenidae$Interspecific_values`,` metrics.sciaenidae$Intraspecific_values` and` metrics.sciaenidae$Summary_table` objects

![Image of Ulises-Rosas](https://github.com/Ulises-Rosas/DNA-Barcodes-BasicMetrics/blob/master/plots.png)

Since we are dowloading sequences directly from the _GenBank_ repository rather than _BOLD_ repository, misannotation of sequences can occur. Nevertheless, the absence of a Barcoding Gap can be defined either in terms of sequence annotation as also effect of broad geographical scales geographical scales. 
