# Guillem Ylla, Ph.D. 2020
# Extavour Lab, Harvard University
library(rtracklayer)
library("circlize")
library(ggtree)
library("viridis")


## Read scaffold lengths from file
Scaff_len<-read.delim("Gbi_Contig_len.txt", sep="\t",stringsAsFactors = FALSE,header = TRUE)
Scaff_len_sorted<-Scaff_len[order(Scaff_len$End, decreasing = TRUE),]
seqlengths<-Scaff_len[,3]
names(seqlengths)<-Scaff_len[,1]


## Set Genome window size
genomewidows = tileGenome(seqlengths, tilewidth=100000, cut.last.tile.in.chrom=T)# split genome 10 000 

############ Genes
## read gff and format names
GbiGffgenes0<-read.delim("Gb_V3.onlyGenes.GBI.3.gff", header = FALSE, sep="\t", stringsAsFactors = FALSE, skip = 1)
GbiGffgenes_2<-GbiGffgenes0[GbiGffgenes0$V3=="gene",c(1,4,5,9)]
GbiGffgenes_2$name0<-sapply(strsplit(GbiGffgenes_2$V9, ";"), "[[", 1)
GbiGffgenes_2$name<-sapply(strsplit(GbiGffgenes_2$name0, "="), "[[", 2)
GbiGffgenes<-GbiGffgenes_2[,c(1,2,3,6)]
#rename columns
colnames(GbiGffgenes)<-c("scaff",  "start"  ,"end",  "name")


gff3genes<-import.gff3("Gb_V3.onlyGenes.GBI.3.gff")
#select genes
gff3genes<-subset(gff3genes, type =="mRNA"  )## only genes with mRNA annotation 
gff3genes


seqlengths<-Scaff_len[,3]
names(seqlengths)<-Scaff_len[,1]
seqlengths

# count genes inside genome windows
Overlapsgene<-findOverlaps(genomewidows, gff3genes)
UniqOverlapsgene<-Overlapsgene[!duplicated(subjectHits(Overlapsgene))]

# calcualte gene density within windows
genedensity=genomewidows
genedensity$totgenes<-0
genedensity[unique(queryHits(UniqOverlapsgene))]$totgenes<-table(queryHits(UniqOverlapsgene))
genedensity_df<-as.data.frame(genedensity)


genedensity_df%>%
  filter( totgenes >0)%>%
  group_by(seqnames )%>%
  summarize(sum(totgenes))%>%
  arrange(desc( `sum(totgenes)`)) %>%
  print(n = 20)

############ Gene Family: PPks

## Get PPk genes from the gene tree file
PPkTree <- read.nexus( "PPK_7OGS")
Gbi_PPkgenes<-PPkTree$tip.label[grep("Gbi",PPkTree$tip.label)]
length(Gbi_PPkgenes)
Gbi_PPkgenes_names0<-sapply(strsplit(Gbi_PPkgenes, "Gbi_"), "[[", 2)
Gbi_PPkgenes_names<-sapply(strsplit(Gbi_PPkgenes_names0, "-"), "[[", 1)

Ppk_coordinates<-GbiGffgenes[GbiGffgenes$name %in%Gbi_PPkgenes_names,]
rownames(Ppk_coordinates)<-NULL
Ppk_coordinates$value<-rep(0, nrow(Ppk_coordinates))

Ppk_coordinatesGR<-makeGRangesFromDataFrame(Ppk_coordinates, seqnames.field="scaff")
Ppk_coordinatesGR

# Find overlaps with genome windows
OverlapsPpk<-findOverlaps(genomewidows, Ppk_coordinatesGR)
UniqOverlapsPpk<-OverlapsPpk[!duplicated(subjectHits(OverlapsPpk))]
Ppkgenes_density=genomewidows
Ppkgenes_density$totgenes<-0
Ppkgenes_density[unique(queryHits(UniqOverlapsPpk))]$totgenes<-table(queryHits(UniqOverlapsPpk))
Ppkgenes_density_df<-as.data.frame(Ppkgenes_density)

############ CpGoe
## Get CpGoe per CDS
CpGoe_Gbi<-read.delim("Gryllus_bimaculatus_longest_CDS_CpGoe.csv")
CpGoe_Gbi$name<-sapply(strsplit(as.character(CpGoe_Gbi$ID), '-'), "[[", 1)
CpGoe_Gbi_list<-CpGoe_Gbi[,c("name","CpGoe")]
head(CpGoe_Gbi_list)
dim(CpGoe_Gbi_list)
CpGoe_Gbi_coordinates<-merge(GbiGffgenes, CpGoe_Gbi_list, all.y=TRUE, by.x="name", by.y="name")
dim(CpGoe_Gbi_coordinates)
head(CpGoe_Gbi_coordinates)


############ TEs
gff3ALL<-import.gff3("Gb_V3.all.GBI.3.gff")
gff3ALL
#select genes
gff3Rep<-subset(gff3ALL, source =="repeatmasker" & type =="match" )
gff3Rep$Name2<-sapply(strsplit(gff3Rep$Name, "genus:" ), `[`, 2)
gff3Rep$Class<-sapply(strsplit(gff3Rep$Name2, '/' ), `[`, 1)
table(gff3Rep$Class)

# calculate overlaps withe genome windows
OverlapsRep<-findOverlaps(genomewidows, gff3Rep)
UniqOverlapsRep<-OverlapsRep[!duplicated(subjectHits(OverlapsRep))]

RepDensity=genomewidows
RepDensity$RepCounts<-0
RepDensity[unique(queryHits(UniqOverlapsRep))]$RepCounts<-table(queryHits(UniqOverlapsRep))
RepDensity_df<-as.data.frame(RepDensity)



GffOthers=gff3Rep[gff3Rep$Class =="DNA?"  | gff3Rep$Class =="RC" |gff3Rep$Class =="rRNA" |gff3Rep$Class =="SINE?" |gff3Rep$Class =="snRNA" |gff3Rep$Class =="Unknown" ]
OthersDensity=genomewidows
OthersDensity$counts = countOverlaps(genomewidows, GffOthers)
OthersDensity_df<-as.data.frame(OthersDensity)



head(RepDensity_df[order(RepDensity_df$RepCounts, decreasing = TRUE),])

RepMeandensity_chr<-aggregate(RepDensity_df$RepCounts, by=list(Category=RepDensity_df$seqnames), FUN=mean)
head(RepMeandensity_chr[order(RepMeandensity_chr$x, decreasing = TRUE),])




library(RColorBrewer)
Selectedcolors<-brewer.pal(n = 8, name = 'Dark2')


##########################
##  Circos Plot
###########################

png("Circos_Plot_Gbi.png",  bg = "transparent", res=400, height = 6, width = 6, units = "in")

# select gene sin scaffolds to plot
GbiGffgenes_selectedregion<-GbiGffgenes[GbiGffgenes$scaff %in% Scaff_len_sorted$Seq, ]


# Initialize Circus
circos.clear()
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 0, gap.after=0)
circos.genomicInitialize(Scaff_lenLOTS, major.by=1000000000, tickLabelsStartFromZero=FALSE, labels.cex=0.00000000001) #major.by avery 1Gb5


# Gene density track
genedensity_df_region<-genedensity_df[genedensity_df$seqnames %in%Scaff_lenLOTS$Seq,c(1,2,3,6)]

circos.genomicTrackPlotRegion(genedensity_df_region ,  numeric.column = 4,
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, type = "l", col =Selectedcolors[2], ...)
                              }, track.height = 0.2,bg.lty=0,bg.lwd=0.00001, bg.border = NA)




PPkgenes_density_region<-Ppkgenes_density_df[Ppkgenes_density_df$seqnames %in%Scaff_lenLOTS$Seq,c(1,2,3,6)]
PPkgenes_density_region$gene<-"ppk"

Genestoplot<-rbind(PPkgenes_density_region )


circos.genomicTrack(Genestoplot, numeric.column = 4,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = value[[1]]/2,col= ifelse(value[[2]] =="ppk", Selectedcolors[5], Selectedcolors[6]), ...)
                    }, track.height = 0.05,bg.lty=0,bg.lwd=0.00001, bg.border = NA, track.margin=c(0,0.03))



## CpGoe
Gbi_cpg_threshold=0.749
CpGoe_High<-CpGoe_Gbi_coordinates[CpGoe_Gbi_coordinates$CpGoe>Gbi_cpg_threshold ,2:5]
CpGoe_Low<-CpGoe_Gbi_coordinates[CpGoe_Gbi_coordinates$CpGoe<=Gbi_cpg_threshold, 2:5]

circos.genomicRainfall(list(CpGoe_High,CpGoe_Low ), pch =16,stack = FALSE,
                       cex =0.5, col = c("#E69F00", "#56B4E9"), track.height = 0.2,
                       bg.lty=0,bg.lwd=0.00001, bg.border = NA)




# rep density
RepDensity_df_region<-RepDensity_df[RepDensity_df$seqnames %in%Scaff_lenLOTS$Seq,c(1,2,3,6)]


circos.genomicTrackPlotRegion(RepDensity_df_region ,  numeric.column = 4,
                              panel.fun = function(region, value, ...) {
                                #i = getI(...)
                                circos.genomicLines(region, value, type = "l", col =Selectedcolors[1], ...)
                              } , track.height = 0.2,bg.lty=0,bg.lwd=0.00001, bg.border = NA)




# add color track scaff >N50 and >N90
circos.track(ylim = c(0, 1),
             bg.col = c(rep(Selectedcolors[3], 71),rep(Selectedcolors[4], 307-71),rep(Selectedcolors[8], nrow(Scaff_lenLOTS)-307) ), # add N50 scafffolds colored
             bg.border = NA, track.height = 0.05)


dev.off()      

