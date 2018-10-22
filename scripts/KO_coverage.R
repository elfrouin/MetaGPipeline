#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                          Preparation                                  ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# set wd
setwd("~/Documents/22-LostProny/clean_version/scripts_normalization/ko_data")
#libraries
require(edgeR)
require(plyr)
require(gplots)
require(RColorBrewer)
require(dendextend)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                           Import data                                 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## import template data
ko_name= read.table("K_complete.tab",sep='\t', comment.char = "",quote = "",fill = NA) 
rownames(ko_name)= ko_name$V1

## import results data
# Get the files names
files = list.files(pattern="*reads_perko.csv$")
# First apply read.csv, table of frequecies ,then rbind
myfiles = lapply(files, function(x) read.table(x, header=TRUE,dec='.',quote = '', comment.char = '',sep = '\t'))

# concatenate results data
data= data.frame(join_all(myfiles,by="ko_hit", type="full"))
# remove sample hydr.GC (because it's related to ultramafic rocks)
data= subset(data, select = -hist_mean_sum_ERR868031.ERR868053)
data= subset(data, select = -hist_mean_sum_SRR4101184)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                      Sampling information                             ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
LC75= c("3862_1325_1st_elution_TGTGAA_L002_reads_perko.csv","LC75","hydrothermal_chimney")
LC81= c("Serp_LC_H08_2_2_DCO_BRZ_GTGGCC_L002_reads_perko.csv","LC81","hydrothermal_chimney")
P27= c("Prony_27_ACACGA_L003_reads_perko.csv","P27","hydrothermal_chimney")
P28= c("Prony_28_GGATGT_L002_reads_perko.csv","P28","hydrothermal_chimney")
hydr.SG =c("ERR868087-ERR868115_reads_perko.csv","hydr.SG","hydrothermal_vents")
hydr.GC=c("ERR868031-ERR868053_reads_perko.csv","hydr.GC","hydrothermal_vents")
Mrk113=c("ERR694207-ERR694213_reads_perko.csv","Mrk113","hydrothermal_vents")
Mrk33=c("ERR694199-ERR694206_reads_perko.csv","Mrk33","hydrothermal_vents")
TAT.2=c("SRR3961741_reads_perko.csv","TAT.2","hot_spring")
TAT.3=c("SRR3961742_reads_perko.csv","TAT.3","hot_spring")
TAT.4=c("SRR3961743_reads_perko.csv","TAT.4","hot_spring")
spring.O=c("SRR4030098_reads_perko.csv","spring.O","hot_spring")
spring.C=c("SRR4030102_reads_perko.csv","spring.C","hot_spring")
GOR.13=c("SRR1636513_reads_perko.csv","GOR.13","travertine")
GOR.12=c("SRR2058407_reads_perko.csv","GOR.12","travertine")
BR2.13=c("SRR1636510_reads_perko.csv","BR2.13","travertine")
BR2.12=c("SRR2058405_reads_perko.csv","BR2.12","travertine")
CdV=c("SRR1636515_reads_perko.csv","CdV","travertine")
CSW11=c("CR12Aug_CSW11AC_reads_perko.csv","CSW11","travertine")
QV11=c("CR12Aug_QV11A_reads_perko.csv","QV11","travertine")
CSW13=c("CR12Aug_CSW13A_reads_perko.csv","CSW13","travertine")
QV12=c("CR12Aug_QV12A_reads_perko.csv","QV12","travertine")
SE.9=c("SRR4101185_reads_perko.csv","SE.9","travertine")
                            spring.O,spring.C,GOR.13,GOR.12,BR2.13,BR2.12,CdV,CSW11,CSW13,QV11,QV12)))
samples=as.data.frame(t(cbind(LC75,LC81,P27,P28,hydr.SG,Mrk113,Mrk33,TAT.3,TAT.4,
                              spring.O,spring.C,GOR.13,GOR.12,BR2.13,BR2.12,CdV,CSW11,CSW13,QV11,QV12,SE.9)))
colnames(samples)= c("full_name", "ID","type")
cols_type= c("orange","seagreen3","dodgerblue3","red4","khaki1")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                          Normalization                                ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# concatenate data
df=data[,-1]
df[ is.na(df) ] <- 0 

group <- files
d <- DGEList(df)
d <- calcNormFactors(d,method = 'TMM')
cps <- cpm(d, normalized.lib.sizes=TRUE)
DATA_norm=data.frame(cps)
rownames(DATA_norm)<- data$ko_hit

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                      Function for heatmap                             ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
plot_heatmap <- function (data,dendo,title,labrow=rownames(data))
  {
  heatmap.2(as.matrix(data),
          labCol=as.character(samples[na.omit(match(files,samples$full_name)),2]),
          margins = c(6,6),
          labRow = labrow,
          ColSideColors=cols_type[as.numeric(samples[na.omit(match(files,samples$full_name)),3])],
          distfun=function(x) as.dist(1-cor(t(x))),
          hclustfun=function(x) hclust(x,method="ward.D2"),
          scale='row',
          trace='none',
          dendrogram = dendo,
          sepcolor="grey",
          col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(15),
          keysize = 1,key=TRUE,main=title,lhei = c(2, 8),symkey = FALSE,symbreaks= FALSE)
  legend(y=0.9, x=-.10, xpd=TRUE,     
       legend = levels(samples$type),
       col =  cols_type,
       pch = 15,           
       cex=0.9,
       bty='n',
       title= "Environments"
       )
  }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                    Function for little heatmap                        ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
plot_little_heatmap <- function (data,dendo,scale,title)
{
  heatmap.2(as.matrix(data),
            labCol=as.character(samples[na.omit(match(files,samples$full_name)),2]),
            margins = c(6,26),
            ColSideColors=cols_type[as.numeric(samples[na.omit(match(files,samples$full_name)),3])],
            distfun=function(x) as.dist(1-cor(t(x))),
            hclustfun=function(x) hclust(x,method="average"),
            scale=scale,
            trace='none',
            dendrogram = dendo,
            sepcolor="grey",
            col=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
            keysize = 1,key=TRUE,main=title)
  legend(y=0.9, x=-.10, xpd=TRUE,     
         legend = levels(samples$type),
         col =  cols_type,
         pch = 15,           
         cex=0.9,
         bty='n',
         title= "Environments"
  )
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
####                                 Plots                                 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# total plot
plot_heatmap(DATA_norm,'col',"Clustering with correlation distance and average-linkage method 
                              of normalized counts of Ko annotations (edgeR)",'')


#with Ko filetred
data_filtered =DATA_norm[rowSums(DATA_norm==0)<=18&(rowSums(DATA_norm==0)<10| (rowSums(DATA_norm==0)>=10&rowSums(DATA_norm) >20)),]
plot_heatmap(data_filtered,'col',"Clustering with correlation distance and Ward-linkage method 
                              of normalized counts of Ko annotations (edgeR)",'')
#colnames(data_filtered)=as.character(samples[na.omit(match(files,samples$full_name)),2])
# Save the data_filtered object
colnames(data_filtered)= as.character(samples[na.omit(match(files,samples$full_name)),2])
saveRDS(data_filtered, "ko_filtered_clean_09-17.rds")
