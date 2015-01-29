#this assumes one is working with the exact same format that fastoldmets/gradoldmets has
#grouped first by tissue, then by slice, then by rep
#adds up BS + Int + M tissues, then takes the mean and the standard error of the replicates
absSEsep<-function(x,nslice=2,ntissue=3,nrep=4){
	if (nslice*ntissue*nrep != length(x[1,])){
		print('Trouble on deck! check your format')
		print(length(x[1,]))
	}
	step=nslice*nrep
	mean.out<-c()
	se.out<-c()
	for (slice in 1:nslice){
		slicesum<-c()
		sadd <- slice - 1
		for(rep in 1:nrep){
			coi<-c(0,step,step*2) + rep + (sadd * nrep)
			slicesum <- cbind(slicesum,apply(x[,coi],1,sum))*2/ntissue
		}
		#print(head(slicesum))
		slicemean<-apply(slicesum,1,mean)
		slicese<-apply(slicesum,1,function(y) sd(y)/sqrt(nrep))
		mean.out<-cbind(mean.out,slicemean)
		se.out<-cbind(se.out,slicese)
	}
	return(list(mean.out,se.out))
}
