### Plot depths and heterozygosities for a pair of strains for each scaffold

options(scipen=999)

data<-read.table("300_R1_val.versus.Pero_masked_genome_V1.fasta.pileup.bins.csv", header=T)

par(mfrow=c(3,4))

for (scaffold.name in c('Scaffold_00002', 'Scaffold_00003', 'Scaffold_00010', 'Scaffold_01882', 'Scaffold_01981', 'Scaffold_02656', 'Scaffold_02555',  'Scaffold_02600', 'Scaffold_01474',
'Scaffold_02202', 'Scaffold_00514', 'Scaffold_02408'  ) ) {
   
#filename <- paste(scaffold.name, 'png', sep='.')
#png(filename=filename, height=900, width=1800)


### Plot depth of coverage
plot(data$midpoint[data$scaffold==scaffold.name],
     data$mean_depth[data$scaffold==scaffold.name]/median(data$mean_depth), type='l', col='darkgreen', xlab='', ylab='', main=paste(scaffold.name, sep=' '), log='y', lwd=3 )

abline(h=1, lty=6)
abline(h=0.5,lty=3)
abline(h=2,lty=3)


### Plot frequency of bin50 heterozygous sites
lines(data$midpoint[data$scaffold==scaffold.name],
     data$bin50[data$scaffold==scaffold.name], type='l', col='black' 
)

#dev.off()

}

### Heatmap 
loh<-read.table("LOH.ramorum1.csv", header=T)
attach(loh)
names(loh)
library(RColorBrewer)
hmcol<-colorRampPalette(brewer.pal(10,"Blues"))(256)
library(pheatmap)
x<-loh[,2:21]
x<-data.matrix(x)
row.names(x)<-loh$scaffold
png(filename="heatmap-heterozygosity.png", width=1200, height=1600)
pheatmap(x, col=hmcol, display_numbers=T, fontsiz_row=6, fontsize_col=24, cellheight=20)
dev.off()








### Plot depths and heterozygosities for a pair of strains for each scaffold

strain1<-read.table("2275.versus.ramorum1.fasta.masked.rmdup.pileup.csv", header=F)
strain2<-read.table("CC1011.versus.ramorum1.fasta.masked.rmdup.pileup.csv", header=F)



scaffold <- 'scaffold_7'

filename <- paste('CC2275', 'CC1011', scaffold, 'png', sep='.')
#start <- 320000
#stop <- 380000

start <- 310000
stop <- 340000

png(filename=filename, height=900, width=1800)
par(mfrow=c(2,1))

plot(strain1$V2[strain1$V1==scaffold], strain1$V3[strain1$V1==scaffold]/median(strain1$V3), type='l',  col='darkgreen', xlab='', ylab='', main=paste('CC2275 (green) and CC1011 (purple) Depth',scaffold,sep=' '), log='y',  lwd=3, xlim=c(start,stop) )

lines(strain2$V2[strain2$V1==scaffold], strain2$V3[strain2$V1==scaffold]/median(strain2$V3), type='l', col='purple', xlab='', ylab='', lwd=3)

abline(h=1, lty=6)
#abline(h=0.5,lty=3 )

plot(strain1$V2[strain1$V1==scaffold], strain1$V4[strain1$V1==scaffold]/10, type='l',  col='darkgreen' , xlab ='', , ylab='', main=paste('CC2275 (green) and CC1011 (purple) Heterozygosity (50:50)',scaffold,sep=' '), lwd=3,   xlim=c(start,stop) )
lines(strain2$V2[strain2$V1==scaffold], strain2$V4[strain2$V1==scaffold]/10, type='l',  col='purple', lwd=3    )

dev.off()

