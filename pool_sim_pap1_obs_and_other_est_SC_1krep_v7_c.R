#script_name: pool_sim_pap1_obs_and_other_est_SC_1krep_v7_c.r
#project: 4-way decomp: paper 1
#script author: Teri North
#script purpose: pool estimates across simulation repeats by
#                 -taking the mean betahat & SE of betahats (to generate MC 95% CI for betahat) 
#                 -take the mean SE and the SD of betahats
#                 -calculate power, type i error and coverage where applicable
#date created: 15/01/2019
#last edited: 18/01/2019
#notes:

setwd('') #Folder 1

install.packages('car')
library('car')

#number of repeats in each sim
repeats=100

#headers
headers=c('mediator_coeff','interac_coeff', 'F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation')

for (nval in c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,500000,1000000)){
  
  meanF=c(1:25)
  for (i in c(1:25)){
    meanF[i]=NA
  }
  
  meanz1F=c(1:25)
  for (i in c(1:25)){
    meanz1F[i]=NA
  }
  

  #calculating the mean F statistics
  first=1
  

  for (seedval in c(7821897,8376154,649384402,238140535,170379645,312006101,713795870,169378934,456561608,28335714)){
    
    for (rep in c(1:repeats)){
      
      if (first==1){
        
        data=read.table(file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_F_res.txt",sep=''),sep='\t',header=TRUE)
        true_vals=data.frame(data$x_coeff_m,data$xm_coeff_y)
        first=0
        
      } else if (first==0){
        
        new=read.table(file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_F_res.txt",sep=''),sep='\t',header=TRUE)
        data=data+new
        
        
      }
      
      
      
    }
    
  }
  
  #remove true values
  data_est=data.frame(data$F_2SLS,data$F_Z1_2SLS)
  
  #mean betas
  mean_denom=repeats*10 #no. rep within seeds * no. seeds
  
  mean_F_stats=data_est/mean_denom

  #add in the true params
  mean_F_stats=cbind(true_vals,mean_F_stats)
  mean_F_stats=rbind(headers,mean_F_stats)

  write.table(mean_F_stats,file=paste(nval,'_mean_F_stats.txt',sep=''),sep='\t',row.names=FALSE,quote=FALSE,col.names=FALSE)

}



#now read back in, de-duplicate and analyse across sample sizes
F_10K=read.table(file="10000_mean_F_stats.txt",sep='\t',header=TRUE)
F_20K=read.table(file="20000_mean_F_stats.txt",sep='\t',header=TRUE)
F_30K=read.table(file="30000_mean_F_stats.txt",sep='\t',header=TRUE)
F_40K=read.table(file="40000_mean_F_stats.txt",sep='\t',header=TRUE)
F_50K=read.table(file="50000_mean_F_stats.txt",sep='\t',header=TRUE)
F_60K=read.table(file="60000_mean_F_stats.txt",sep='\t',header=TRUE)
F_70K=read.table(file="70000_mean_F_stats.txt",sep='\t',header=TRUE)
F_80K=read.table(file="80000_mean_F_stats.txt",sep='\t',header=TRUE)
F_90K=read.table(file="90000_mean_F_stats.txt",sep='\t',header=TRUE)
F_100K=read.table(file="1e+05_mean_F_stats.txt",sep='\t',header=TRUE)
F_500K=read.table(file="5e+05_mean_F_stats.txt",sep='\t',header=TRUE)
F_1000K=read.table(file="1e+06_mean_F_stats.txt",sep='\t',header=TRUE)

F_10K=data.frame(F_10K$mediator_coeff,F_10K$F_statistic_2sls_assuming_no_mediation,F_10K$F_statistic_2sls_assuming_mediation)
F_20K=data.frame(F_20K$mediator_coeff,F_20K$F_statistic_2sls_assuming_no_mediation,F_20K$F_statistic_2sls_assuming_mediation)
F_30K=data.frame(F_30K$mediator_coeff,F_30K$F_statistic_2sls_assuming_no_mediation,F_30K$F_statistic_2sls_assuming_mediation)
F_40K=data.frame(F_40K$mediator_coeff,F_40K$F_statistic_2sls_assuming_no_mediation,F_40K$F_statistic_2sls_assuming_mediation)
F_50K=data.frame(F_50K$mediator_coeff,F_50K$F_statistic_2sls_assuming_no_mediation,F_50K$F_statistic_2sls_assuming_mediation)
F_60K=data.frame(F_60K$mediator_coeff,F_60K$F_statistic_2sls_assuming_no_mediation,F_60K$F_statistic_2sls_assuming_mediation)
F_70K=data.frame(F_70K$mediator_coeff,F_70K$F_statistic_2sls_assuming_no_mediation,F_70K$F_statistic_2sls_assuming_mediation)
F_80K=data.frame(F_80K$mediator_coeff,F_80K$F_statistic_2sls_assuming_no_mediation,F_80K$F_statistic_2sls_assuming_mediation)
F_90K=data.frame(F_90K$mediator_coeff,F_90K$F_statistic_2sls_assuming_no_mediation,F_90K$F_statistic_2sls_assuming_mediation)
F_100K=data.frame(F_100K$mediator_coeff,F_100K$F_statistic_2sls_assuming_no_mediation,F_100K$F_statistic_2sls_assuming_mediation)
F_500K=data.frame(F_500K$mediator_coeff,F_500K$F_statistic_2sls_assuming_no_mediation,F_500K$F_statistic_2sls_assuming_mediation)
F_1000K=data.frame(F_1000K$mediator_coeff,F_1000K$F_statistic_2sls_assuming_no_mediation,F_1000K$F_statistic_2sls_assuming_mediation)

F_10K=unique(F_10K)
F_20K=unique(F_20K)
F_30K=unique(F_30K)
F_40K=unique(F_40K)
F_50K=unique(F_50K)
F_60K=unique(F_60K)
F_70K=unique(F_70K)
F_80K=unique(F_80K)
F_90K=unique(F_90K)
F_100K=unique(F_100K)
F_500K=unique(F_500K)
F_1000K=unique(F_1000K)

n_10k=c(10000,10000,10000,10000,10000)
n_20k=c(20000,20000,20000,20000,20000)
n_30k=c(30000,30000,30000,30000,30000)
n_40k=c(40000,40000,40000,40000,40000)
n_50k=c(50000,50000,50000,50000,50000)
n_60k=c(60000,60000,60000,60000,60000)
n_70k=c(70000,70000,70000,70000,70000)
n_80k=c(80000,80000,80000,80000,80000)
n_90k=c(90000,90000,90000,90000,90000)
n_100k=c(100000,100000,100000,100000,100000)
n_500k=c(500000,500000,500000,500000,500000)
n_1000k=c(1000000,1000000,1000000,1000000,1000000)

F_10K=cbind(F_10K,n_10k)
F_20K=cbind(F_20K,n_20k)
F_30K=cbind(F_30K,n_30k)
F_40K=cbind(F_40K,n_40k)
F_50K=cbind(F_50K,n_50k)
F_60K=cbind(F_60K,n_60k)
F_70K=cbind(F_70K,n_70k)
F_80K=cbind(F_80K,n_80k)
F_90K=cbind(F_90K,n_90k)
F_100K=cbind(F_100K,n_100k)
F_500K=cbind(F_500K,n_500k)
F_1000K=cbind(F_1000K,n_1000k)


colnames(F_10K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_20K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_30K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_40K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_50K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_60K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_70K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_80K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_90K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_100K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_500K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')
colnames(F_1000K)=c('mediator_coeff','F_statistic_2sls_assuming_no_mediation','F_statistic_2sls_assuming_mediation','N')

ALLDAT=rbind(F_10K,F_20K,F_30K,F_40K,F_50K,F_60K,F_70K,F_80K,F_90K,F_100K,F_500K,F_1000K)

ALLDAT=ALLDAT[order(ALLDAT$mediator_coeff,ALLDAT$N),]
#write table to file
write.table(ALLDAT,file="F_stats_final_table.txt",sep='\t',row.names=FALSE,quote=FALSE,col.names=TRUE)

#plots
setwd('')

tiff('SWF_nomed.tif',width=7,height=3.5,units='in',res=400)
scatterplot(ALLDAT$N,ALLDAT$F_statistic_2sls_assuming_no_mediation,groups=ALLDAT$mediator_coeff,xlab='sample size',
            ylab='S-W F Statistic, Z=(Z1,Z2,Z1Z2)',legend=list(title='alpha',cex=0.5),cex.lab=0.5,cex.axis=0.5)
abline(h=10)
dev.off()

tiff('SWF_med.tif',width=7,height=3.5,units='in',res=400)
scatterplot(ALLDAT$N,ALLDAT$F_statistic_2sls_assuming_mediation,groups=ALLDAT$mediator_coeff,xlab='sample size',
            ylab='S-W F Statistic, Z=(Z1,Z2,Z1Z2,Z1Z1)',legend=list(title='alpha',cex=0.5),cex.lab=0.5,cex.axis=0.5)
abline(h=10)
dev.off()

smallsamp=ALLDAT[ALLDAT$N<=100000,]

tiff('SWF_nomed_smallsamp.tif',width=7,height=3.5,units='in',res=400)
scatterplot(smallsamp$N,smallsamp$F_statistic_2sls_assuming_no_mediation,groups=smallsamp$mediator_coeff,xlab='sample size',
            ylab='S-W F Statistic, Z=(Z1,Z2,Z1Z2)',legend=list(title='alpha',cex=0.5),cex.lab=0.5,cex.axis=0.5)
abline(h=10)
dev.off()

tiff('SWF_med_smallsamp.tif',width=7,height=3.5,units='in',res=400)
scatterplot(smallsamp$N,smallsamp$F_statistic_2sls_assuming_mediation,groups=smallsamp$mediator_coeff,xlab='sample size',
            ylab='S-W F Statistic, Z=(Z1,Z2,Z1Z2,Z1Z1)',legend=list(title='alpha',cex=0.5),cex.lab=0.5,cex.axis=0.5)
abline(h=10)
dev.off()
