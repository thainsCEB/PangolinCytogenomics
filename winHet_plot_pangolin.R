library(plyr)
library(gtools)
library(showtext)

winsize=1e6
stepsize=1e6

plotwinhet=function(searchstring, myname){
  window=sprintf("%.0f",winsize)
  stepsize=sprintf("%.0f",stepsize)
  ending=paste("_",window,"win_",stepsize,"step.txt",sep="")
	hetfiles=list.files(pattern=ending)
	hetfiles=mixedsort(hetfiles)
#	hetfiles=hetfiles[c(1:4,6:30)]
	allhet=ldply(hetfiles, read.table, header=TRUE, sep="\t")
	temp=allhet[which(allhet[,4]>=(0.5*winsize)),]
	
	meanhet=sum(temp[,5])/sum(temp[,4])
	plotname=paste("mean het. = ", sprintf("%.4f", meanhet), sep="")
	
	pos=as.numeric(rownames(unique(data.frame(temp$chrom)[1])))
	pos=append(pos,length(temp$chrom))
	numpos=NULL
	for (i in 1:length(pos)-1){numpos[i]=(pos[i]+pos[i+1])/2}

	mycols=NULL
	for (i in (seq(1,length(numpos), by=2))){mycols[i]="#005a32"}
	for (i in (seq(2,length(numpos), by=2))){mycols[i]="#238b45"}
	
  par(mar=c(6,5,4,1))
	b=barplot(temp[,5]/temp[,4], ylim=c(0,.015), border=mycols[as.factor(temp$chrom)], col=mycols[as.factor(temp$chrom)], ylab="Heterozygosity", main=myname,las=1,cex.axi=0.9,font.main=3,family="Helvetica")
	axis(side=1, at=b[pos], labels=F)
	axis(side=1, at=b[numpos], tick=F, labels=mixedsort(as.character(unique(temp$chrom))), las=1, line=-.5)
	mtext(plotname,cex=1)
}

setwd("/Users/taylorhains/PhaTri1_Het/")
searchstring="PhaTri1"
myname="Phataginus tricuspsis"
pdf(paste("winHet_1Mbwin_1Mbstep_",searchstring,".pdf", sep=""), width=15, height=5)
plotwinhet(searchstring, myname)
dev.off()

setwd("/Users/taylorhains/ManJav1_Het/")
searchstring="ManJav1"
myname="Manis javanica"
pdf(paste("winHet_1Mbwin_1Mbstep_",searchstring,".pdf", sep=""), width=9, height=3)
plotwinhet(searchstring, myname)
dev.off()

setwd("/Users/taylorhains/ManPen1_Het/")
searchstring="ManPen1"
myname="Manis pentadactyla"
pdf(paste("winHet_1Mbwin_1Mbstep_",searchstring,".pdf", sep=""), width=15, height=3)
plotwinhet(searchstring, myname)
dev.off()

