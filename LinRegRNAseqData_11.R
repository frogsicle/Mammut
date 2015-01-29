###funtion to take whole linear regressions and spit out intercepts (fraction in M or BS) for each gene across whole leaf

ill.normalization<-function(yourdata,start=3){
	end<-length(yourdata[1,])
	regressions=matrix(data=NA,nrow=4,end)
	colnames(regressions)<-colnames(yourdata)
	row.names(regressions)<-c("intercept(M)","intercept(BS)","p-value","stdev")
	for (i in start:end)
	{
		if (!(any(is.na(yourdata[,i]))))
		{
			linregdelta<-lm(yourdata[,i] ~ yourdata[,(start-1)])
			linregsig<-summary(lm((yourdata[,i]-0.5) ~ yourdata[,(start-1)]))
			regressions[1,i]<-linregdelta[[1]][[1]]
			regressions[2,i]<-1-regressions[1,i]
			regressions[3,i]<-linregsig[[4]][[7]]
			regressions[4,i]<-linregsig[[4]][[3]]

		}
	}	
	return (regressions)
}

###function to take linear regressions of each replicate, and return just the intercept, representing the fraction in the mesophyll
##note, this is very much set up for M S T being the most external sort function
repill.normalization<-function(yourdata,start=3,spacing=15){
        end<-length(yourdata[1,])
        regressions=matrix(data=NA,nrow=spacing,end)
        colnames(regressions)<-colnames(yourdata)
        row.names(regressions)<-c("AD_14","AD_19","AD_24","AD_15","AD_20","AD_25","AD_16","AD_21","AD_26","AD_17","AD_22","AD_27","AD_18","AD_23","AD_28")
        for (i in start:end)
        {
		for (j in 1:spacing)
		{
			k=j+spacing
			l=k+spacing
	
        	        if (!(any(is.na(yourdata[c(j,k,l),i]))))
                	{
                        	linregdelta<-lm(yourdata[c(j,k,l),i] ~ yourdata[c(j,k,l),(start-1)])
                        	linregsig<-summary(lm((yourdata[c(j,k,l),i]-0.5) ~ yourdata[c(j,k,l),(start-1)]))
                        	regressions[j,i]<-linregdelta[[1]][[1]]
                	}
		}	
        }
        return (regressions)
}



###function to take linear regressions of each segment(with all bioreps), and return just the intercept, representing the fraction in the mesophyll
##note, this is very much set up for M S T being the most external sort function, and then segments, and then reps
my_sublmslope<-function(my_x,my_coi,BSoM)
	{
	if (length(which(my_x[my_coi]>0))>0)
		{
	        return(lm(log(my_x[my_coi][which(my_x[my_coi]>0)])~log(BSoM[my_coi][which(my_x[my_coi]>0)]))$coefficients[[2]])
		}
	else
		{
		return(NA)
		}
	}

#adjusted to spit out p-val as well 
my_sublmslope_p<-function(my_x,my_coi,BSoM)
        {
        if (length(which(my_x[my_coi]>0))>0)
                {
		#offset by *0.5 to adjust slope against which the p.val is tested to 0.5
		tmp.lm<-lm(log(my_x[my_coi][which(my_x[my_coi]>0)])~log(BSoM[my_coi][which(my_x[my_coi]>0)]),offset=log(BSoM[my_coi][which(my_x[my_coi]>0)])*0.5)
		#due to offset, the slope needs to be corrected by 0.5
		tmp.slope<-tmp.lm$coefficients[[2]]+0.5
		#null hypothesis of pval is now slope = 0.5
		tmp.pval<-summary(tmp.lm)$coefficients[[8]]
		#the standard error for the sake of plotting
		tmp.se<-summary(tmp.lm)$coefficients[[4]]

                return(c(tmp.slope,tmp.pval,tmp.se))
                }
        else
                {
                return(c(NA,NA,NA))
                }
        }

#attempting to trouble shoot whether I have a marker or a raw data difference
my_sublmslope_t<-function(my_x,my_coi,BSoM)
        {
        if (length(which(my_x[my_coi]>0))>0)
                {
		#offset by *0.5 to adjust slope against which the p.val is tested to 0.5
		test<-length(my_x[my_coi][which(my_x[my_coi]>0)])
		tmp.lm<-lm(log(my_x[my_coi][which(my_x[my_coi]>0)])~jitter(rep(1,test)))
		#due to offset, the slope needs to be corrected by 0.5
		tmp.slope<-tmp.lm$coefficients[[2]]+0.5
		#null hypothesis of pval is now slope = 0.5
		tmp.pval<-summary(tmp.lm)$coefficients[[8]]
                return(c(tmp.slope,tmp.pval))
                }
        else
                {
                return(c(NA,NA))
                }
        }



#major normalization of fraction in BS function
flexi_segs_logi<-function(yourdata,reps=3,spacing=15,coln=c("tip","mid-tip","mid","base-mid","base"))
	{
	counting=0
	pdata<-matrix(NA,length(yourdata[,1]),spacing/reps)
	row.names(pdata)<-row.names(yourdata)
	colnames(pdata)<-coln
	for (i in seq(from=1,to=spacing,by=reps))
		{
		counting=counting+1
		extra=reps-1
		j<-i+spacing
		k<-i+spacing*2
		coi<-c(c(i:(i+extra)),c(j:(j+extra)),c(k:(k+extra)))
		out<-apply(yourdata,1,function(x) my_sublmslope(my_x=x,my_coi=coi,BSoM=yourdata["BS_markers",]))
		pdata[,counting]<-unlist(out)
		}
	return(pdata)
	}
#major norm now with significance too
flexi_segs_logstats<-function(yourdata,reps=3,spacing=15,coln=c("tip","mid-tip","mid","base-mid","base"))
        {
        counting=0
        pdata<-matrix(NA,length(yourdata[,1]),spacing/reps)
	sdata<-matrix(NA,length(yourdata[,1]),spacing/reps)
        row.names(pdata)<-row.names(yourdata)
        colnames(pdata)<-coln
        for (i in seq(from=1,to=spacing,by=reps))
                {
                counting=counting+1
                extra=reps-1
                j<-i+spacing
                k<-i+spacing*2
                coi<-c(c(i:(i+extra)),c(j:(j+extra)),c(k:(k+extra)))
		print(coi)
                out<-apply(yourdata,1,function(x) my_sublmslope_p(my_x=x,my_coi=coi,BSoM=yourdata["BS_markers",]))
		print(dim(out))
		print(dim(pdata))
                pdata[,counting]<-unlist(out[1,])
		sdata[,counting]<-unlist(out[2,])
                }
        return(cbind(pdata,sdata))
        }
###major norm w/o specific samples
#enter a 0 in colsamp for any sample you wish to skip
custom_logstats<-function(yourdata,colsamp=c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5),coln=c("S1","S2","S3","S4","S5"))
        {
        counting=0
        pdata<-matrix(NA,length(yourdata[,1]),max(colsamp))
	sdata<-matrix(NA,length(yourdata[,1]),max(colsamp))
        edata<-matrix(NA,length(yourdata[,1]),max(colsamp))
        row.names(pdata)<-row.names(yourdata)
        colnames(pdata)<-coln
        for (i in seq(from=1,to=max(colsamp),by=1))
                {
                coi<-which(colsamp==i)
		#print(coi)
                out<-apply(yourdata,1,function(x) my_sublmslope_p(my_x=x,my_coi=coi,BSoM=yourdata["BS_markers",]))
		#print(dim(out))
		print(head(out))
		#print(dim(pdata))
                pdata[,i]<-unlist(out[1,])
		sdata[,i]<-unlist(out[2,])
		edata[,i]<-unlist(out[3,])
                }
        return(cbind(pdata,sdata,edata))
        }


#norm outdated
segill.normalization<-function(yourdata,start=3,spacing=15){
        end<-length(yourdata[1,])
        regressions=matrix(data=NA,nrow=spacing/3,end)
        colnames(regressions)<-colnames(yourdata)
        row.names(regressions)<-c("tip","mid-tip","mid","Base-mid","base")
        for (i in start:end)
        {
		counting=0
                for (j in seq(from=1,to=spacing,by=3))
                {
			counting=counting+1
                        k=j+spacing
                        l=k+spacing
			coi<-c(c(j:(j+2)),c(k:(k+2)),c(l:(l+2)))
                        if (!(any(is.na(yourdata[coi,i]))))
                        {
                                linregdelta<-lm(yourdata[coi,i] ~ yourdata[coi,(start-1)])
                                linregsig<-summary(lm((yourdata[coi,i]-0.5) ~ yourdata[coi,(start-1)]))
                                regressions[counting,i]<-linregdelta[[1]][[1]]
                        }
                }
        }
        return (regressions)
}

