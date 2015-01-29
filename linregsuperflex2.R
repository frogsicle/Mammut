### this function deconvolutes 'yourdata' in different slices/subsets (denoted by colsamp, and named by coln)
#set colsamp to 0 to skip a replicate
#note, obviously mishanding of global variables like overm... but don't want to mess without testing, see pubnorm.R
custoutlie_logstats<-function(yourdata,colsamp=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5),coln=c("tip","mid-tip","mid","base-mid","base"))
        {
	#outliers to NAs
	spread<-as.vector(yourdata)[which(as.vector(overm)!=0)]
	mysd<-sd(log(spread))
	#toremove<-log(spread)[which(abs(log(spread)-mean(log(spread)))>3*mysd)]
	print('removing outliers > 3sd')
	print(log(yourdata[which(abs(log(yourdata)-mean(log(spread)))>3*mysd)]))
	yourdata[which(abs(log(yourdata)-mean(log(spread)))>3*mysd)]<-0
	#prep for deconvoluting
        counting=0
        pdata<-matrix(NA,length(yourdata[,1]),max(colsamp))
	sdata<-matrix(NA,length(yourdata[,1]),max(colsamp))
        edata<-matrix(NA,length(yourdata[,1]),max(colsamp))
        row.names(pdata)<-row.names(yourdata)
        colnames(pdata)<-coln
        for (i in seq(from=1,to=max(colsamp),by=1))
                {
                coi<-which(colsamp==i)
                out<-apply(yourdata,1,function(x) my_sublmslope_p(my_x=x,my_coi=coi,BSoM=yourdata["BS_markers",]))
                pdata[,i]<-unlist(out[1,])
		sdata[,i]<-unlist(out[2,])
                edata[,i]<-unlist(out[3,])
                }
        return(cbind(pdata,sdata,edata))
        }



