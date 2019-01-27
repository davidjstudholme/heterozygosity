### Load the date from .csv files (generated using get_freqs_of_first_and_second_allele.pl)
phals<-read.table("300bp.versus.Phals_V2_repeatmasked_Nuclear_genome.fasta.rmdup.pileup.freq-first-and-second-most-abundant.csv", header=F)
peb<-read.table("300_R1_val.versus.Pero_masked_genome_V1.fasta.pileup.freq-first-and-second-most-abundant.csv", header=F)
emoy2<-read.table("Emoy2.combined.versus.GCA_000173235.2_HyaAraEmoy2_2.0_genomic.aln.sorted.bam.pileup.freq-first-and-second-most-abundant.csv", header=F)
infestans<-read.table("ERR248813.versus.P_infestans.GCA_000142945.1_ASM14294v1_genomic.pileup.freq-first-and-second-most-abundant.csv", header=F)
kerno<-read.table("Kerno238_431.versus.GCA_000785735.2_NZFS2646v2_genomic.pileup.freq-first-and-second-most-abundant.csv", header=F)
lateralis<-read.table("MPF6.versus.GCA_000318465.2_MPF4_v2.0_genomic.pileup.freq-first-and-second-most-abundant.csv", header=F)
ramorum<-read.table("SOD_22_12_GTAGAG_L008_R1_001_val.versus.ramorum1.fasta.pileup.freq-first-and-second-most-abundant.csv",header=F)
parasitica<-read.table("SRR058717.versus.S_parasitica.GCA_000151545.2_ASM15154v2_genomic.pileup.freq-first-and-second-most-abundant.csv", header=F)
sojae<-read.table("SRR1046799.versus.P_sojae.GCA_000149755.2_P.sojae_V3.0_genomic.pileup.freq-first-and-second-most-abundant#.csv", header=F)
albugo<-read.table("SRR1811471.versus.A_candida.CAIW01000001.1.pileup.freq-first-and-second-most-abundant.csv", header=F)
ptab<-read.table("SRR2146895.versus.Ptab_968-S26_published_assembly.fasta.rename.fasta.pileup.freq-first-and-second-most-abundant.csv", header=F)
ultimum<-read.table("SRR235484.versus.P_ultimum.AKYB02000001.1.pileup.freq-first-and-second-most-abundant.csv", header=F)
vexans<-read.table("SRR235485.versus.P_vexans.AKYC02000001.1.pileup.freq-first-and-second-most-abundant.csv", header=F)
cubensis<-read.table("SRR412826.versus.P_cubensis.AHJF01000001.1.pileup.freq-first-and-second-most-abundant.csv", header=F)
capsici<-read.table("SRR945695.versus.P_capsici.ADVJ01000001.1.pileup.freq-first-and-second-most-abundant.csv", header=F)

### Define function
plot.heterozygosity.histogram <- function(genome, genome.title, upper.lim, lower.lim, num.classes, min.depth, show.legend){
most.abundant <- genome$V3[genome$V3 >= lower.lim & genome$V3 <= upper.lim & !is.na(genome$V3) & genome$V5 >= min.depth ]
second.most.abundant <- genome$V4[ genome$V4 >= lower.lim & genome$V4 <= upper.lim & !is.na(genome$V4) & genome$V5 >= min.depth ]

if (show.legend == T) {
   	hist(second.most.abundant,
		#nclass=num.classes,
		main = genome.title,
		xlab='Proportion', 
		ylab='Density',
		ylim=c(0,8),	
		xlim=c(lower.lim, upper.lim),
		col = rgb(1,0,0,0.4),
		breaks = num.classes,
		prob=T
		)
	} else {
	 hist(second.most.abundant,
		#nclass=num.classes,
		main = genome.title,
		xlab='', 
		ylab='',
		ylim=c(0,8),
		xlim=c(lower.lim, upper.lim),
		col = rgb(1,0,0,0.4),
		breaks = num.classes,
		prob=T
		)
	}

hist(most.abundant,
	add=T,
	#nclass=num.classes,
	col = rgb(0,0,1,0.4),
	breaks = num.classes,	
	prob=T
)

adjust <- 1

lines( density( most.abundant, adjust=adjust), col = rgb(0,0,1,0.4), lwd=4 )
lines( density( second.most.abundant, adjust=adjust), col = rgb(1,0,0,0.4), lwd=4 )

return(1)

}

### Do stuff
upper.lim <- 0.90
lower.lim <- 0.10
num.classes <- 16
min.depth <- 10

pdf('~/Dropbox/Peb.heterozygosity.pdf')

### Multiple graphs on the page
op<-par(mfrow=c(4,4))

#dummy <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#hist(dummy, main="", xlab='', ylab='', yaxt='n', xaxt='n')
#legend('left',c('Second-most','Most'), fill = rgb(1:0,0,0:1,0.4), bty = 'n', border = NA, inset=0)
                                
plot.heterozygosity.histogram(peb, "Pe. belbahrii", upper.lim, lower.lim, num.classes, min.depth, T)
plot.heterozygosity.histogram(ptab, "Pe. tabacina", upper.lim, lower.lim, num.classes, min.depth, F)

plot.heterozygosity.histogram(cubensis, "Ps. cubensis", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(emoy2, "H. arabidopsidis", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(phals, "Pl. halstedii", upper.lim, lower.lim, num.classes, min.depth, F)

plot.heterozygosity.histogram(capsici, "Ph. capsici", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(infestans, "Ph. infestans", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(kerno, "Ph. kernoviae", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(lateralis, "Ph. lateralis", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(parasitica, "Ph. parasitica", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(ramorum, "Ph. ramorum", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(sojae, "Ph. sojae", upper.lim, lower.lim, num.classes, min.depth, F)

plot.heterozygosity.histogram(ultimum, "Py. ultimum", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(vexans, "Py. vexans", upper.lim, lower.lim, num.classes, min.depth, F)
plot.heterozygosity.histogram(albugo, "A. candida", upper.lim, lower.lim, num.classes, min.depth, F)

dummy <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
hist(dummy, main="", xlab='', ylab='', yaxt='n', xaxt='n')
legend('left',c('Second-most','Most frequent'), fill = rgb(1:0,0,0:1,0.4), bty = 'n', border = F, inset=0, cex=0.85)

dev.off()

### Revert par
par(op)
