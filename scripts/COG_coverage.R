setwd("~/Documents/23-Pipeline_snakemake/final_pipeline/norm_coverage/")

require(edgeR)
require(plyr)
require(gplots)
require(RColorBrewer)
library(dendextend)

files = list.files(pattern="*COG.tsv")
myfiles = lapply(files, function(x) read.table(x, header=TRUE,dec='.',quote = '', comment.char = '',sep = '\t'))
data= data.frame(join_all(myfiles,by="cog_hit", type="full"))

df= data[,-c(1,2,3,4)]
df[is.na(df)] <- 0
group <- files
d <- DGEList(df)
d <- calcNormFactors(d,method = 'TMM')
cps <- cpm(d, normalized.lib.sizes=TRUE)
DATA_norm=data.frame(cps)
rownames(DATA_norm)<- paste(data$cog_hit,'[',data$cog_class,']', sep='')

LC75= c("3862_1325_1st_elution_TGTGAA_L002_covCOG.tsv","LC75","hydrothermal_chimney")
LC81= c("Serp_LC_H08_2_2_DCO_BRZ_GTGGCC_L002_covCOG.tsv","LC81","hydrothermal_chimney")
P27= c("Prony_27_ACACGA_L003_covCOG.tsv","P27","hydrothermal_chimney")
P28= c("Prony_28_GGATGT_L002_covCOG.tsv","P28","hydrothermal_chimney")
hydr.SG =c("ERR868087-ERR868115_covCOG.tsv","hydr.SG","hydrothermal_vents")
#hydr.GC=c("ERR868031-ERR868053_covCOG.tsv","hydr.GC","hydrothermal_vents")
Mrk113=c("ERR694207-ERR694213_covCOG.tsv","Mrk113","hydrothermal_vents")
Mrk33=c("ERR694199-ERR694206_covCOG.tsv","Mrk33","hydrothermal_vents")
TAT.2=c("SRR3961741_covCOG.tsv","TAT.2","hot_spring")
TAT.3=c("SRR3961742_covCOG.tsv","TAT.3","hot_spring")
TAT.4=c("SRR3961743_covCOG.tsv","TAT.4","hot_spring")
spring.O=c("SRR4030098_covCOG.tsv","spring.O","hot_spring")
spring.C=c("SRR4030102_covCOG.tsv","spring.C","hot_spring")
GOR.13=c("SRR1636513_covCOG.tsv","GOR.13","travertine")
GOR.12=c("SRR2058407_covCOG.tsv","GOR.12","travertine")
BR2.13=c("SRR1636510_covCOG.tsv","BR2.13","travertine")
BR2.12=c("SRR2058405_covCOG.tsv","BR2.12","travertine")
CdV=c("SRR1636515_covCOG.tsv","CdV","travertine")
CSW11=c("CR12Aug_CSW11AC_covCOG.tsv","CSW11","travertine")
QV11=c("CR12Aug_QV11A_covCOG.tsv","QV11","travertine")
CSW13=c("CR12Aug_CSW13A_covCOG.tsv","CSW13","travertine")
QV12=c("CR12Aug_QV12A_covCOG.tsv","QV12","travertine")
SE.9=c("SRR4101185_covCOG.tsv","SE.9","travertine")
samples=as.data.frame(t(cbind(LC75,LC81,P27,P28,hydr.SG,Mrk113,Mrk33,TAT.3,TAT.4,
               spring.O,spring.C,GOR.13,GOR.12,BR2.13,BR2.12,CdV,CSW11,CSW13,QV11,QV12,SE.9)))
colnames(samples)= c("full_name", "ID","type")
cols_type= c("orange","seagreen3","dodgerblue3","red4","khaki1")

heatmap.2(as.matrix(DATA_norm),
          labCol=as.character(samples[na.omit(match(files,samples$full_name)),2]),
          margins = c(5,15),labRow = '',
          ColSideColors=cols_type[as.numeric(samples[na.omit(match(files,samples$full_name)),3])],
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x,method="average"),
          symkey=FALSE,
          scale='row',
          trace='none',
          dendrogram = 'col',
          sepcolor="grey",
          cexRow = 1.2,
          col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
          keysize = 1,key=TRUE, main ="Clustering with correlation distance and average-linkage method 
           of normalized counts of COG annotations (edgeR)")

legend(y=0.9, x=-.10, xpd=TRUE,     
       legend = levels(samples$type),
       col =  cols_type,
       pch = 15,           
       cex=0.9,
       bty='n',
       title= "Environments"
       )

# with filtered data
data_filtered =DATA_norm[rowSums(DATA_norm==0)<=18&(rowSums(DATA_norm==0)<10| (rowSums(DATA_norm==0)>=10&rowSums(DATA_norm) >13.76)),]
# 13.76 = fisrt Qu.
colnames(data_filtered)= as.character(samples[na.omit(match(files,samples$full_name)),2])
saveRDS(data_filtered, "cog_data_09-17.rds")

# save non-normalized data
data_save= data[,-c(2,3,4)]
colnames(data_save)= c('COG',as.character(samples[na.omit(match(files,samples$full_name)),2]))
data_save[is.na(data_save)] <- 0
write.table(data_save, file = "samples_reads_percog.txt",sep = "\t",quote = FALSE, row.names = FALSE)

heatmap.2(as.matrix(data_filtered),
          labCol=as.character(samples[na.omit(match(files,samples$full_name)),2]),
          margins = c(5,15),labRow = '',
          ColSideColors=cols_type[as.numeric(samples[na.omit(match(files,samples$full_name)),3])],
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x,method="ward.D2"),
          scale='row',
          trace='none',
          dendrogram = 'col',
          sepcolor="grey",
          cexRow = 1.2,
          col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(15),
          keysize = 1,key=TRUE, main ="Clustering with correlation distance and Ward linkage method 
           of normalized counts of COG annotations (edgeR)",lhei = c(2, 8),symkey = FALSE,symbreaks= FALSE)

legend(y=0.9, x=-.10, xpd=TRUE,     
       legend = levels(samples$type),
       col =  cols_type,
       pch = 15,           
       cex=0.9,
       bty='n',
       title= "Environments"
)
# with high normalized counts
r1=DATA_norm[which(rowSums(DATA_norm>3000)>0),]
heatmap.2(as.matrix(r1),
          labCol=as.character(samples[na.omit(match(files,samples$full_name)),2]),
          margins = c(5,15),
          ColSideColors=cols_type[as.numeric(samples[na.omit(match(files,samples$full_name)),3])],
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x,method="average"),
          scale='row',
          trace='none',
          sepcolor="grey",
          cexRow = 1.2,
          col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
          keysize = 1,key=TRUE)

## select COG with high sd in Prony samples
d_scale = as.data.frame(scale(DATA_norm))
S= d_scale[which(((d_scale$hist_mean_sum_Prony_27_ACACGA_L003>=apply(d_scale,1, max)|
                   d_scale$hist_mean_sum_Prony_28_GGATGT_L002>=apply(d_scale,1, max))&
                  (d_scale$hist_mean_sum_Prony_27_ACACGA_L003>2 |
                   d_scale$hist_mean_sum_Prony_28_GGATGT_L002>2))),]
heatmap.2(as.matrix(S),
          labCol=as.character(samples[na.omit(match(files,samples$full_name)),2]),
          margins = c(6,16),
          ColSideColors=cols_type[as.numeric(samples[na.omit(match(files,samples$full_name)),3])],
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x,method="average"),
          scale='row',
          trace='none',
          sepcolor="grey",
          cexRow = 1.2,
          col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
          keysize = 1,key=TRUE)


## select COG with high sd in LC samples
d_scale = as.data.frame(t(scale(t(DATA_norm))))
#d_scale = as.data.frame(scale(DATA_norm))
S= d_scale[which(((d_scale$hist_mean_sum_3862_1325_1st_elution_TGTGAA_L002>=apply(d_scale,1, max)|
                  d_scale$hist_mean_sum_Serp_LC_H08_2_2_DCO_BRZ_GTGGCC_L002>=apply(d_scale,1, max))&
                  (d_scale$hist_mean_sum_3862_1325_1st_elution_TGTGAA_L002>2.5 &
                  d_scale$hist_mean_sum_Serp_LC_H08_2_2_DCO_BRZ_GTGGCC_L002>2.5))),]
heatmap.2(as.matrix(S),
          labCol=as.character(samples[na.omit(match(files,samples$full_name)),2]),
          margins = c(6,16),
          ColSideColors=cols_type[as.numeric(samples[na.omit(match(files,samples$full_name)),3])],
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x,method="average"),
          scale='row',
          trace='none',
          sepcolor="grey",
          cexRow = 1.2,
          col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
          keysize = 1,key=TRUE)
# 
S3= d_scale[which(d_scale$hist_mean_sum_ERR694199.ERR694206>0.5 &
                 d_scale$hist_mean_sum_ERR694207.ERR694213>0.5 &
                 d_scale$hist_mean_sum_ERR868031.ERR868053>0.5 &
                 d_scale$hist_mean_sum_ERR868087.ERR868115>0.5 &
                 d_scale$hist_mean_sum_3862_1325_1st_elution_TGTGAA_L002>0.5 &
                 d_scale$hist_mean_sum_Serp_LC_H08_2_2_DCO_BRZ_GTGGCC_L002>0.5),]

## Comparing two dendrograms
# Compute distance matrix
res.dist <- as.dist(1-cor((DATA_norm)))
# Compute 2 hierarchical clusterings
hc1 <- hclust(res.dist, method = "average")
hc2 <- hclust(res.dist, method = "ward.D2")
# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)
# Create a list of dendrograms
dend_list <- dendlist(dend1, dend2)
#A lower entanglement coefficient corresponds to a good alignment
tanglegram(dend1, dend2,  common_subtrees_color_branches = TRUE, # Color common branches ,
           margin_inner=20,
           main = paste("entanglement =", round(entanglement(dend_list), 2)))
