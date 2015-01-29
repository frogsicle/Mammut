
##get the fractions of whatever you dump in
subfrac<-function(x)
	{
	if(min(x)>0)
		{
		my_frac<-x/sum(x)	
		}
	else
		{
		my_frac<-x*0
		}
	return(my_frac)
	}
###get the sums of just the segments used for normalization

fraction_ize<-function(x){
	y<-x
	third<-length(x[1,])/3
	for (i in 1:third)
	{
		z<-c(i,(i+third),(i+third*2))
		y[,z]<-t(apply(x[,z],1,subfrac))
	}	
	return (y)
}

seg_sums<-function(x)
	{
	y<-matrix(NA,length(x[,1]),15)
	row.names(y)<-row.names(x)
	colnames(y)<-c("AD_14","AD_19","AD_24","AD_15","AD_20","AD_25","AD_16","AD_21","AD_26","AD_17","AD_22","AD_27","AD_18","AD_23","AD_28")
        for (i in 1:15)
        	{
		for (j in 1:length(x[,1]))
			{
			y[j,i]<-sum(na.omit(x[j,c(i,(i+15),(i+30))]))
			}
        	}
        return (y)
	}



#a more flexible summation, that already matches the format out from flexi_segs_illnorm
flexi_segs_sums<-function(yourdata,reps=3,spacing=15)
        {
        counting=0
        pdata<-matrix(NA,length(yourdata[,1]),spacing/reps)
        row.names(pdata)<-row.names(yourdata)
        colnames(pdata)<-c("tip","mid-tip","mid","base-mid","base")
        for (i in seq(from=1,to=spacing,by=reps))
                {
                counting=counting+1
                extra=reps-1
                j<-i+spacing
                k<-i+spacing*2
                coi<-c(c(i:(i+extra)),c(j:(j+extra)),c(k:(k+extra)))
                out<-apply(yourdata[,coi],1,sum)/3
                pdata[,counting]<-unlist(out)
                }
        return(pdata)
        }




