#######################################################################################
#
# Calculate 8 single filter statistics and 3 summary statistics for Methylation data (beta values)
#
# 1) SD
# 2) SD/SDmax
# 3) mad
# 4) distance from theoretical center after beta transformation
# 5) quantile difference from theoretical CDF distribution
# 6) quantile difference from theoretical PDF distribution
# 7) dip statistics for unimordality
# 8) using SD.M rank (logit.x <- log(x.adj/(1-x.adj), base=2)
# 9) using best rank out of the above 8 ranks as raw value
# 10) using average rank of the above 8 ranks as raw value
# 11) using different weighted average of the above 8 ranks as raw value
# 12) combine half SD-b + half TM-GOF top probes
#
# by Xinhui Wang, Kim Siegmund
#
#######################################################################################

library(diptest)

###############################################################################
## 
## start with non-missing methylation data set x (beta values) 
##        
###############################################################################

nr.x=nrow(x) #samples
nc.x=ncol(x) #probes


## get means of Xj for each gene probe j (column)
mean.x <- apply(x,2,mean,na.rm=T)


## 
## Calculate SD-b and MAD
##
SD.x <- apply(x,2,sd,na.rm=T)
mad.x <- apply(x,2,mad,na.rm=T)


## 
## Calculate SDob/SDmax (1/precision)
##
SD.x.max <- sqrt(mean.x*(1-mean.x))
Ratio.x <- SD.x/SD.x.max
Ratio.x[1:20]


##
## calculate SD-m
##
# are there boundary values? if yes, shrinkage them to be between (0,1)
sum(x==0,na.rm=T)
sum(x==1,na.rm=T)
x.adj = x
x.adj = ifelse(x.adj==0,1e-8,x.adj)
x.adj = ifelse(x.adj==1,1-1e-8,x.adj)

logit.x <- log(x.adj/(1-x.adj), base=2)


# get SD of logit transformed Xj for each gene probe j (column)
SD.logit.x <- apply(logit.x,2,sd,na.rm=T)



##
##  Do beta transformation and calculate three statistics (TM-GOF, BQ-GOF, TQ-GOF)
##  
# Using beta distribution parameters caculated from mean and variance to transform the data 
alpha_hat.x <- mean.x**2*(1-mean.x)/SD.x**2-mean.x
beta_hat.x <- mean.x*(1-mean.x)**2/SD.x**2-1+mean.x
para.x <- rbind(alpha_hat.x, beta_hat.x)

phi.x <- alpha_hat.x+beta_hat.x


#
# Use grand (alpha, beta) to transform dataset x to y
#
cx=rbind(x,para.x)
my.y=apply(cx,2,function(x){
 			pbeta(x[1:nr.x],x[nr.x+1],x[nr.x+2],lower.tail=TRUE,log.p=FALSE)
 			}	
 		)
mean.y <- apply(my.y,2,mean,na.rm=T)
SD.y <- apply(my.y,2,sd,na.rm=T)

# calculate the variance of mean.y and SD.y
SD.mean.y <- sd(mean.y, na.rm = TRUE)
SD.SD.y <- sd(SD.y, na.rm = TRUE)


#
# calculate the standardized distance of mean.y and SD.y to the theoretical center point (TM-GOF)
#
dist.y <- sqrt(((mean.y-0.5)/SD.mean.y)**2+((SD.y-1/sqrt(12))/SD.SD.y)**2)


#
# calculate sum of quantiles using pdf (BQ-GOF)
#
# calculate 26 quantiles of estimated beta distribution for each probe 
observed.x.qts.pdf <- apply(x,2,function(fx){
		quantile(fx, probs = seq(0, 1, 0.04),na.rm = TRUE, names = TRUE)}
)
expected.x.qts.pdf <- apply(para.x,2,function(fx){
	qbeta(p=seq(0, 1, 0.04),fx[1],fx[2], ncp = 0, lower.tail = TRUE, log.p = FALSE)}
)

pdf.quantile.sum <- apply(abs(observed.x.qts.pdf-expected.x.qts.pdf),2,sum,na.rm=T)


#
# calculate sum of quantiles using cdf (TQ-GOF)
#
# calculate 26 quantiles of transformed beta distribution for each probe 
observed.x.qts.cdf <- apply(my.y,2,function(fx){
		quantile(fx, probs = seq(0, 1, 0.04),na.rm = TRUE, names = TRUE)}
)
expected.x.qts.cdf <- array(seq(0,1,0.04),dim(observed.x.qts))

cdf.quantile.sum <- apply(abs(observed.x.qts.cdf-expected.x.qts.cdf),2,sum,na.rm=T)


##
## calculate dip for each probe, bigger means more informative (DIP)
##
dip.x <- apply(x,2,dip)

single.method=cbind.data.frame(mean.x,SD.x,SD.logit.x,mad.x,dip.x,1/phi.x,
				pdf.quantile.sum,dist.y,cdf.quantile.sum)



###############################################################################
##
##  ADD SUMMARY FILTERS
##
###############################################################################

rank.single.methods=apply(single.method,2,rank)
rank.eight.filts=rank.single.methods[,c("SD.x","SD.logit.x","mad.x","dip.x","1/phi.x",
					"pdf.quantile.sum","dist.y","cdf.quantile.sum")]

o.rank.eight.filts=rank.eight.filts
colnames(o.rank.eight.filts)=paste("r",1:8,sep="")
for (i in 1:nrow(o.rank.eight.filts)) o.rank.eight.filts[i,]=rev(sort(rank.eight.filts[i,]))
head(o.rank.eight.filts)


#
# calculate best rank (BR), average rank (AR), weighted average rank (WAR) statistics
#
BR=apply(rank.eight.filts,1,max)
AR=apply(o.rank.eight.filts[,1:2],1,mean)
WAR=o.rank.eight.filts[,1:4]%*%as.matrix(4:1)/10

sum.filts=cbind.data.frame(BR,AR,WAR)
rank.sum.filts=apply(sum.filts,2,rank)
rmean.x=rank(single.method$mean.x)

## matrix includes mean in last position
rank.eleven.filts=cbind.data.frame(rank.eight.filts,rank.sum.filts,rmean.x)
colnames(rank.eleven.filts)=c("SDb","SDm","MAD","DIP","InvPr","BQ-GOF",
					"TM-GOF","TQ-GOF","BR","AR","WR","Mean")
head(rank.eleven.filts)



#
# filter data set using top N probes: row=probe; col=sample 
# 

nc.x=ncol(x)
npick=1000

dat.SDb.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="SDb")]>(nc.x-npick)]
dat.SDm.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="SDm")]>(nc.x-npick)]
dat.mad.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="MAD")]>(nc.x-npick)]
dat.dip.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="DIP")]>(nc.x-npick)]
dat.InvPr.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="InvPr")]>(nc.x-npick)]
dat.BQ.GOF.filt		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="BQ-GOF")]>(nc.x-npick)]
dat.TM.GOF.filt 	<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="TM-GOF")]>(nc.x-npick)]
dat.TQ.GOF.filt 	<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="TQ-GOF")]>(nc.x-npick)]
dat.BR.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="BR")]>(nc.x-npick)]
dat.AR.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="AR")]>(nc.x-npick)]
dat.WR.filt 		<- x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="WR")]>(nc.x-npick)]
dat.SDb.TM.GOF.filt 	<- cbind(x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="SDb")]>(nc.x-npick/2)],
				x[, rank.eleven.filts[,which(colnames(rank.eleven.filts)=="TM-GOF")]>(nc.x-npick/2)])


