#import functions
#load('~/Dropbox/Libelle2/mammut/data/illumina16/CleanPub/Normfunctions.R')
source('LinRegRNAseqData_11.R')
source('numerize_11.R')
source('calc_absSE.R')
source('linregsuperflex2.R')

#reformatting is of course, necessary
lababsmunits<-read.csv(file="milliUnitsAct_ed.csv")#this file has milliUnits enzyme activity / mgFW already calculated
mUEnzyme<-lababsmunits[order(lababsmunits[,'Fraction']),]
num_mUEnzyme<-t(mUEnzyme[,6:length(mUEnzyme[1,])])

rename<-c()
for (i in 1:120){
	rename<-c(rename,paste(mUEnzyme[i,2],mUEnzyme[i,1],mUEnzyme[i,3],sep='_'))
}
colnames(num_mUEnzyme)<-rename
#TODO collapse technical replicates before calculating significance please! manu
enzact <- t(num_mUEnzyme)
enzact2 <- c()
for (i in 0:59){
        enzact2<-rbind(enzact2,colMeans(enzact[(i*2+1:2),]))
}
rownames(enzact2)<-colnames(num_mUEnzyme)[0:59*2+1]#apply(enzact[0:59*2+1,1:5],1,function(x) paste(x[1:3],collapse='_'))
enzact2<-t(enzact2)

#calculate tempoary variables for normalization (gene):
mx2<-enzact2['PEPC',]
bx2<-enzact2['NADPME',]
overm<-rbind(fraction_ize(rbind(bx2,mx2,enzact2)))
rownames(overm)[1:2]<-c('BS_markers','M_markers')
overm<-t(t(overm)/overm['M_markers',])
my_segsums<-flexi_segs_sums(enzact2,reps=4,spacing=20)
#colnames(my_segsums)<-colnames(my_segsums)[5:1

#save gene
#enznormout<-custom_logstats(overm,colsam=rep(c(rep(1,8),rep(2,8),rep(3,8),rep(4,8),rep(5,8)),3))
enznormout<-custoutlie_logstats(overm,colsam=rep(c(rep(1,4),rep(2,4),rep(3,4),rep(4,4),rep(5,4)),3))
#enznormout<-custom_logstats(overm,reps=8,spacing=40)
normenz_BSM<-cbind(enznormout[3:length(enznormout[,1]),1:5]*my_segsums,(1-enznormout[3:length(enznormout[,1]),1:5])*my_segsums)
normenz_p<-enznormout[3:length(enznormout[,1]),6:10]
normenz_se<-enznormout[3:length(enznormout[,1]),11:15]
colnames(normenz_p)<-paste('S',1:5,sep='')
colnames(normenz_se)<-paste('S',1:5,sep='')
colnames(normenz_BSM)<-c('BS1','BS2','BS3','BS4','BS5','M1','M2','M3','M4','M5')
#activity level
enzsums<-absSEsep(enzact2,nslice=5,ntissue=3,nrep=4)
colnames(enzsums[[1]])<-paste('S',1:5,sep='')
colnames(enzsums[[2]])<-paste('S',1:5,sep='')

#quick overview plot
barplot(normenz_BSM[1:7,10:1],col=tmp,ylab='deconvoluted mU enzyme activity')
legend('topleft',legend=rownames(normenz_BSM)[1:7],fill=tmp)
dev.copy2eps(file='barenzymes.R')
save(list=c('normenz_BSM','normenz_p'),file='DeconEnzymes.Robject')
save(num_mUEnzyme,file='MarkEnzymes.Robject')
save(list=c('enzsums','num_mUEnzyme','normenz_BSM','normenz_p','normenz_se'),file='enzyme_sum.Robject')



