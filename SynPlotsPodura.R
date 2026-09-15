library(data.table)

## synteny plots for T. podura phased genomes

## read synteny dat
#########################################################################
## first compare with cristinae, green hap 1 from Hwy 154
##############
dat_pod_cris<-fread("out_cactusTpod1_TcrGUSH1.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod_cris)

## target = pod1
## query = gush1 cristinae

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
crisCh<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% crisCh)
subDfdat<-dfdat[keep,]## 13 cristinae, 14 podura

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tc_gush1<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,52,4)])
tpod_1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,111,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,112,8)])## sizes for podura


## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura
}
pdf("SynTcrisStripe_podura_GUSH1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. cristinae (H GUS1)",ylab="T. podura (1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,13,length.out=13)/13,tc_gush1,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## chrom number, gs,gsr1, gsh1, gush1, tpod
chtab<-matrix(c(1,8483,12,13,22,7,
        2,14640,6,5,23,4,
        3,42935,2,3,16,3,
        3,42935,2,3,16,6,
        3,42935,2,3,16,9,
        4,42912,1,1,64,6,
        4,42912,1,1,64,12,
        5,18722,7,12,5,14,
        6,9928,8,4,11,5,
        7,10660,10,10,54,13,
        8,7748,11,11,7,2,
        9,16151,5,8,46,8,
        10,14160,4,7,15,3,
        11,12033,9,9,2,11,
        12,12380,13,6,1,10,
        13,14101,3,2,36,1),nrow=16,ncol=6,byrow=TRUE)


pdf("AlnPlotTcris_podura_GUSH1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:16){
	tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	tcr_gush1<-grep(x=subDfdat[,10],pattern=paste("6S_",chtab[i,5],"_",sep=""))
	cc<-tpo_1[tpo_1 %in% tcr_gush1]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura",ylab="T. cristinae")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tpod_1==chtab[i,6])]-subd[j,16:17])
		}
	}
}
dev.off()

#########################################################################
## now compar podura haps
##############
dat_pod1_pod2<-fread("out_cactusTpod1_TcrPod2.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(unlist(strsplit(x=rownames(tab),split="_"))[seq(2,111,8)])
tpod_1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,111,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,112,8)])## sizes for podura1


## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}
## clean as hell

pdf("SynTcrisStripe_podura1_podura2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (2)",ylab="T. podura (1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,7,9,
        2,14640,6,5,23,4,6,
        3,42935,2,3,16,9,8,
        4,42912,1,1,64,12,2,
        5,18722,7,12,5,14,14,
        6,9928,8,4,11,5,4,
        7,10660,10,10,54,13,11,
        8,7748,11,11,7,2,12,
        9,16151,5,8,46,8,13,
        10,14160,4,7,15,3,1,
        11,12033,9,9,2,11,10,
        12,12380,13,6,1,10,7,
        13,14101,3,2,36,1,3,
	14,NA,NA,NA,NA,6,5),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnPlotTcris_pod1_pod2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	tpo_2<-grep(x=subDfdat[,10],pattern=paste("old_",chtab[i,7],"_",sep="")) 
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura 1",ylab="T. podura 2")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tpod_1==chtab[i,6])]-subd[j,16:17])
		}
	}
}
dev.off()

#########################################################################
## now compare to newer podura haps
############## M x M
dat_pod1_pod2<-fread("out_syn_TpodH1_TpodE240140H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,111,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,112,8)])## sizes for podura1


## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodH1_TpodE240140H1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240140H1)",ylab="T. podura (H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,7,12,
        2,14640,6,5,23,4,9,
        3,42935,2,3,16,9,1,
        4,42912,1,1,64,12,6,
        5,18722,7,12,5,14,2,
        6,9928,8,4,11,5,8,
        7,10660,10,10,54,13,14,
        8,7748,11,11,7,2,7,
        9,16151,5,8,46,8,10,
        10,14160,4,7,15,3,5,
        11,12033,9,9,2,11,13,
        12,12380,13,6,1,10,11,
        13,14101,3,2,36,1,4,
	14,NA,NA,NA,NA,6,3),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodH1_TpodE240140H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura 1",ylab="T. podura E240140H1")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tpod_1==chtab[i,6])]-subd[j,16:17])
		}
	}
}
dev.off()

############## M x M
dat_pod1_pod2<-fread("out_syn_TpodH1_TpodE240140H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,111,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,112,8)])## sizes for podura1


## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodH1_TpodE240140H2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240140H2)",ylab="T. podura (H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,7,12,
        2,14640,6,5,23,4,9,
        3,42935,2,3,16,9,1,
        4,42912,1,1,64,12,6,
        5,18722,7,12,5,14,2,
        6,9928,8,4,11,5,8,
        7,10660,10,10,54,13,14,
        8,7748,11,11,7,2,7,
        9,16151,5,8,46,8,10,
        10,14160,4,7,15,3,5,
        11,12033,9,9,2,11,13,
        12,12380,13,6,1,10,11,
        13,14101,3,2,36,1,4,
	14,NA,NA,NA,NA,6,3),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodH1_TpodE240140H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura 1",ylab="T. podura E240140H2")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tpod_1==chtab[i,6])]-subd[j,16:17])
		}
	}
}
dev.off()

############## M x G
dat_pod1_pod2<-fread("out_syn_TpodH1_TpodE240154H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,111,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,112,8)])## sizes for podura1


## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodH1_TpodE240154H1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240154H1)",ylab="T. podura (H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,7,13,
        2,14640,6,5,23,4,9,
        3,42935,2,3,16,9,2,
        4,42912,1,1,64,12,5,
        5,18722,7,12,5,14,7,
        6,9928,8,4,11,5,8,
        7,10660,10,10,54,13,14,
        8,7748,11,11,7,2,6,
        9,16151,5,8,46,8,11,
        10,14160,4,7,15,3,4,
        11,12033,9,9,2,11,12,
        12,12380,13,6,1,10,10,
        13,14101,3,2,36,1,3,
	14,NA,NA,NA,NA,6,1),nrow=14,ncol=7,byrow=TRUE)

pdf("AlnTpod_TpodH1_TpodE240154H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura 1",ylab="T. podura E240154H1")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tpod_1==chtab[i,6])]-subd[j,16:17])
		}
	}
}
dev.off()

############## M x G
dat_pod1_pod2<-fread("out_syn_TpodH1_TpodE240154H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(2,111,8)])
sizes<-as.numeric(unlist(strsplit(x=colnames(tab),split="_"))[seq(8,112,8)])## sizes for podura1


## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodH1_TpodE240154H2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240154H2)",ylab="T. podura (H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,7,13,
        2,14640,6,5,23,4,9,
        3,42935,2,3,16,9,2,
        4,42912,1,1,64,12,5,
        5,18722,7,12,5,14,7,
        6,9928,8,4,11,5,8,
        7,10660,10,10,54,13,14,
        8,7748,11,11,7,2,6,
        9,16151,5,8,46,8,11,
        10,14160,4,7,15,3,4,
        11,12033,9,9,2,11,12,
        12,12380,13,6,1,10,10,
        13,14101,3,2,36,1,3,
	14,NA,NA,NA,NA,6,1),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodH1_TpodE240154H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura 1",ylab="T. podura E240154H2")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],sizes[which(tpod_1==chtab[i,6])]-subd[j,16:17])
		}
	}
}
dev.off()

## now Edinburgh genomes against each other

############## M x M
dat_pod1_pod2<-fread("out_syn_TpodE240140H1_TpodE240140H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
## sizes approximated below

## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodE240140H1_TpodE240140H2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240140H2)",ylab="T. podura (E240140H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,7,12,
        2,14640,6,5,23,4,9,
        3,42935,2,3,16,9,1,
        4,42912,1,1,64,12,6,
        5,18722,7,12,5,14,2,
        6,9928,8,4,11,5,8,
        7,10660,10,10,54,13,14,
        8,7748,11,11,7,2,7,
        9,16151,5,8,46,8,10,
        10,14160,4,7,15,3,5,
        11,12033,9,9,2,11,13,
        12,12380,13,6,1,10,11,
        13,14101,3,2,36,1,4,
	14,NA,NA,NA,NA,6,3),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodE240140H1_podE240140H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	#tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,7] > 9){
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr",chtab[i,7],sep="")) 
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr0",chtab[i,7],sep="")) 
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura E240140H1",ylab="T. podura E240140H2")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

############## G x G
dat_pod1_pod2<-fread("out_syn_TpodE240154H1_TpodE240154H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
## sizes approximated below

## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodE240154H1_TpodE240154H2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240154H2)",ylab="T. podura (E240154H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,7,13,
        2,14640,6,5,23,4,9,
        3,42935,2,3,16,9,2,
        4,42912,1,1,64,12,5,
        5,18722,7,12,5,14,7,
        6,9928,8,4,11,5,8,
        7,10660,10,10,54,13,14,
        8,7748,11,11,7,2,6,
        9,16151,5,8,46,8,11,
        10,14160,4,7,15,3,4,
        11,12033,9,9,2,11,12,
        12,12380,13,6,1,10,10,
        13,14101,3,2,36,1,3,
	14,NA,NA,NA,NA,6,1),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodE240154H1_podE240154H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	#tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,7] > 9){
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr",chtab[i,7],sep="")) 
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr0",chtab[i,7],sep="")) 
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura E240154H1",ylab="T. podura E240154H2")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

############## M x G
dat_pod1_pod2<-fread("out_syn_TpodE240140H1_TpodE240154H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
## sizes approximated below

## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodE240140H1_TpodE240154H1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240154H1)",ylab="T. podura (E240140H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,12,13,
        2,14640,6,5,23,9,9,
        3,42935,2,3,16,1,2,
        4,42912,1,1,64,6,5,
        5,18722,7,12,5,2,7,
        6,9928,8,4,11,8,8,
        7,10660,10,10,54,14,14,
        8,7748,11,11,7,7,6,
        9,16151,5,8,46,10,11,
        10,14160,4,7,15,5,4,
        11,12033,9,9,2,13,12,
        12,12380,13,6,1,11,10,
        13,14101,3,2,36,4,3,
	14,NA,NA,NA,NA,3,1),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodE240140H1_podE240154H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	#tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,6] > 9){
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr",chtab[i,6],sep="")) 
	} else{
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr0",chtab[i,6],sep="")) 
	}
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura E240140H1",ylab="T. podura E240154H1")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

############## M x G
dat_pod1_pod2<-fread("out_syn_TpodE240140H1_TpodE240154H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
## sizes approximated below

## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodE240140H1_TpodE240154H2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240154H2)",ylab="T. podura (E240140H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,12,13,
        2,14640,6,5,23,9,9,
        3,42935,2,3,16,1,2,
        4,42912,1,1,64,6,5,
        5,18722,7,12,5,2,7,
        6,9928,8,4,11,8,8,
        7,10660,10,10,54,14,14,
        8,7748,11,11,7,7,6,
        9,16151,5,8,46,10,11,
        10,14160,4,7,15,5,4,
        11,12033,9,9,2,13,12,
        12,12380,13,6,1,11,10,
        13,14101,3,2,36,4,3,
	14,NA,NA,NA,NA,3,1),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodE240140H1_podE240154H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	#tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,6] > 9){
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr",chtab[i,6],sep="")) 
	} else{
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr0",chtab[i,6],sep="")) 
	}
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura E240140H1",ylab="T. podura E240154H2")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

############## M x G
dat_pod1_pod2<-fread("out_syn_TpodE240140H2_TpodE240154H1.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
## sizes approximated below

## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodE240140H2_TpodE240154H1.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240154H2)",ylab="T. podura (E240140H1)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,12,13,
        2,14640,6,5,23,9,9,
        3,42935,2,3,16,1,2,
        4,42912,1,1,64,6,5,
        5,18722,7,12,5,2,7,
        6,9928,8,4,11,8,8,
        7,10660,10,10,54,14,14,
        8,7748,11,11,7,7,6,
        9,16151,5,8,46,10,11,
        10,14160,4,7,15,5,4,
        11,12033,9,9,2,13,12,
        12,12380,13,6,1,11,10,
        13,14101,3,2,36,4,3,
	14,NA,NA,NA,NA,3,1),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodE240140H2_podE240154H1.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	#tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,6] > 9){
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr",chtab[i,6],sep="")) 
	} else{
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr0",chtab[i,6],sep="")) 
	}
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura E240140H2",ylab="T. podura E240154H1")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()

############## M x G
dat_pod1_pod2<-fread("out_syn_TpodE240140H2_TpodE240154H2.psl",header=FALSE)
dfdat<-as.data.frame(dat_pod1_pod2)

## target = pod1
## query = pod2

## verify/identify large scaffolds (all should be) 
xx<-table(dfdat[,14])
pod1Ch<-names(xx)[xx>500]
xx<-table(dfdat[,10])
pod2Ch<-names(xx)[xx>500]
keep<-(dfdat[,14] %in% pod1Ch) & (dfdat[,10] %in% pod2Ch)
subDfdat<-dfdat[keep,]## 14 chroms each

tab<-tapply(X=subDfdat[,1],INDEX=list(qg=subDfdat[,10],tg=subDfdat[,14]),sum)
tpod_2<-as.numeric(gsub(x=unlist(strsplit(x=rownames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
tpod_1<-as.numeric(gsub(x=unlist(strsplit(x=colnames(tab),split="_"))[seq(2,28,2)],pattern="Chr",replacement=""))
## sizes approximated below

## normalize  
ntab<-tab
ntab[is.na(ntab)]<-0
for(i in 1:14){
	ntab[,i]<-ntab[,i]/sum(ntab[,i]) ## relative to podura1
}

pdf("SynTpod_TpodE240140H2_TpodE240154H2.pdf",width=6,height=6)
par(mar=c(5,5,1,1))
image(ntab,axes=FALSE,xlab="T. podura (E240154H2)",ylab="T. podura (E240140H2)",cex.lab=1.4)
axis(2,at=seq(0,14,length.out=14)/14,tpod_1,las=2)
axis(1,at=seq(0,14,length.out=14)/14,tpod_2,las=2)
box()
dev.off()

## colinearity plots for all homologous chromsomes
## simplified for podura focus
## chrom number, gs,gsr1, gsh1, gush1, tpod1, tpod2
chtab<-matrix(c(1,8483,12,13,22,12,13,
        2,14640,6,5,23,9,9,
        3,42935,2,3,16,1,2,
        4,42912,1,1,64,6,5,
        5,18722,7,12,5,2,7,
        6,9928,8,4,11,8,8,
        7,10660,10,10,54,14,14,
        8,7748,11,11,7,7,6,
        9,16151,5,8,46,10,11,
        10,14160,4,7,15,5,4,
        11,12033,9,9,2,13,12,
        12,12380,13,6,1,11,10,
        13,14101,3,2,36,4,3,
	14,NA,NA,NA,NA,3,1),nrow=14,ncol=7,byrow=TRUE)


pdf("AlnTpod_TpodE240140H2_podE240154H2.pdf",width=11,height=11)
par(mfrow=c(4,4))
par(mar=c(4.5,5.5,2.5,1.5))
for(i in 1:14){
	#tpo_1<-grep(x=subDfdat[,14],pattern=paste("old_",chtab[i,6],"_",sep="")) 
	if(chtab[i,6] > 9){
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr",chtab[i,6],sep="")) 
	} else{
		tpo_1<-grep(x=subDfdat[,14],pattern=paste("Chr0",chtab[i,6],sep="")) 
	}
	if(chtab[i,7] > 9){
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr",chtab[i,7],sep="")) 
	} else{
		tpo_2<-grep(x=subDfdat[,10],pattern=paste("Chr0",chtab[i,7],sep="")) 
	}
	cc<-tpo_1[tpo_1 %in% tpo_2]
	subd<-subDfdat[cc,]
	xub<-max(subd[,13]);yub<-max(subd[,17])	

	plot(as.numeric(subd[1,12:13]),as.numeric(subd[1,16:17]),type='n',xlim=c(0,xub),ylim=c(0,yub),cex.lab=1.4,xlab="T. podura E240140H2",ylab="T. podura E240154H2")
	title(main=paste("Chromosome",chtab[i,1]),cex.main=1.4)
	N<-dim(subd)[1]
	for(j in 2:N){
		if(subd[j,9]=="++"){
			lines(subd[j,12:13],subd[j,16:17])
		}
		else{
			lines(subd[j,12:13],yub-subd[j,16:17])
		}
	}
}
dev.off()
