#script_name: Pool_meta_sensitivity_v6.r
#project: 4-way decomp: paper 1
#script author: Teri North
#script purpose: pool estimates across simulation repeats by
#                 -taking the mean betahat & SE of betahats (to generate MC 95% CI for betahat) 
#                 -take the mean SE and the SD of betahats
#                 -calculate power, type i error and coverage where applicable
#date created: 09/08/2018
#last edited: 11/10/2018
#notes:

setwd('') #Folder 1

#number of repeats in each sim
repeats=100
nval=500000

#tracker for erroneous calls - will get two numbers one for model 1 and the other for model 2

n_z1problem=c(1:2)
for (j in c(1:2)){n_z1problem[j]=0}
n_z1prob_track=1



for (model in c(1:2)){
  
  xm_z1_2sls_detec=c(1:25)
  for (i in c(1:25)){
    xm_z1_2sls_detec[i]=0
  }
  
  
  z1_coverage=c(1:25)
  for (i in c(1:25)){
    z1_coverage[i]=0
  }
  
  
  #reality check
  #how many times is the interaction detected (p<0.05), but the estimate is in the opposite direction to true effect?

  z1problem=c(1:25) 
  for (i in c(1:25)){
    z1problem[i]=0
  }
  
  
  
  fmr_y_int=c(1:25)  # counter for # times interaction detected factorial approach (Wald test 5%)
  for (i in c(1:25)){
    fmr_y_int[i]=0
  }
  
  
  #calculating the mean betas   
  first=1
  
  
  for (seedval in c(520160447,267639401,37905828,750891730,435580371,945959183,141153971,456264979,86129334,119011473)){
    
    for (rep in c(1:repeats)){
      
      if (first==1){
        
        data=read.table(file=paste(seedval,'_rep',rep,'_model',model,"_pleio_res.txt",sep=''),sep='\t',header=TRUE)
        true_vals=data.frame(data$x_coeff_m,data$x_coeff_y,data$m_coeff_y, data$xm_coeff_y)
        first=0
        
        ll=data.frame(
          xm_2sls_ll=data$xm_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_2sls_se
        )
        
        ul=data.frame(
          xm_2sls_ul=data$xm_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_2sls_se
        )
        
        
        for (i in c(1:25)){
          if (ll$xm_2sls_ll[i]>0 | ul$xm_2sls_ul[i]<0){xm_z1_2sls_detec[i]=1}
        }

        for (i in c(1:25)){
          if ((ll$xm_2sls_ll[i]<data$xm_coeff_y[i]) & (ul$xm_2sls_ul[i]>data$xm_coeff_y[i])){z1_coverage[i]=1}
        }  
        

        for (i in c(1:25)){
          if (((ll$xm_2sls_ll[i]>0) & (data$xm_coeff_y[i]<0))|((ul$xm_2sls_ul[i]<0) & (data$xm_coeff_y[i]>0))) {z1problem[i]=1}#if interac detec, but coeff wrong direc
        } 
        
 
        for (i in c(1:25)){
          if (data$fmr_interac_p[i]<0.05){fmr_y_int[i]=1}
        }
        
        
        
      } else if (first==0){
        
        new=read.table(file=paste(seedval,'_rep',rep,'_model',model,"_pleio_res.txt",sep=''),sep='\t',header=TRUE)
        data=data+new
        
        ll_new=data.frame(
          xm_2sls_ll=new$xm_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_2sls_se
        )
        
        ul_new=data.frame(
          xm_2sls_ul=new$xm_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_2sls_se
        )
        
        
        for (i in c(1:25)){
          if (ll_new$xm_2sls_ll[i]>0 | ul_new$xm_2sls_ul[i]<0){xm_z1_2sls_detec[i]=xm_z1_2sls_detec[i]+1}
        }
        

        for (i in c(1:25)){
          if ((ll_new$xm_2sls_ll[i]<new$xm_coeff_y[i]) & (ul_new$xm_2sls_ul[i]>new$xm_coeff_y[i])){z1_coverage[i]=z1_coverage[i]+1}
        }  
        

        for (i in c(1:25)){
          if (((ll_new$xm_2sls_ll[i]>0) & (new$xm_coeff_y[i]<0))|((ul_new$xm_2sls_ul[i]<0) & (new$xm_coeff_y[i]>0))) {z1problem[i]=z1problem[i]+1}#if interac detec, but coeff wrong direc
        } 

        
        for (i in c(1:25)){
          if (new$fmr_interac_p[i]<0.05){fmr_y_int[i]=fmr_y_int[i]+1}
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  
  
  #remove true values
  data_est=data.frame(data$xm_2sls,data$xm_2sls_se)
  
  
  
  #mean betas
  mean_denom=repeats*10 #no. rep within seeds * no. seeds
  data_mean=data_est/mean_denom #gives mean beta and mean se
  #add in the true params
  mean_betas=cbind(true_vals,data_mean)
  
  
  #############################################################################################################################################################################

  #now for the standard error
  checker=1
  
  for (seeds in  c(520160447,267639401,37905828,750891730,435580371,945959183,141153971,456264979,86129334,119011473)){
    
    for (rep in c(1:repeats)){
      
      if (checker==1){
        
        newdata=read.table(file=paste(seedval,'_rep',rep,'_model',model,"_pleio_res.txt",sep=''),sep='\t',header=TRUE)
        newdata=(newdata-(data/mean_denom))^2
        checker=0
        
      } else if (checker==0){
        
        newer=read.table(file=paste(seedval,'_rep',rep,'_model',model,"_pleio_res.txt",sep=''),sep='\t',header=TRUE)
        newdata=newdata+(newer-(data/mean_denom))^2
        newdata=data.frame(newdata)
      }
      
    }
  }
  
  
  newdata_est=data.frame(newdata$xm_2sls)
  
  
  
  #divide by n-1 to get s^2
  s2=newdata_est/(repeats*10-1)
  se=sqrt(s2/(repeats*10))
  
  
  ##################################################################################################################################################################
  

  xm_z1_2sls_detec=data.frame(xm_z1_2sls_detec)
  z1_coverage=data.frame(z1_coverage)
  fmr_y_int=data.frame(fmr_y_int)
  
  #results table 
  res=data.frame(mean_betas$data.x_coeff_m,
                 mean_betas$data.x_coeff_y,
                 mean_betas$data.m_coeff_y,
                 mean_betas$data.xm_coeff_y,
                 mean_betas$data.xm_2sls,
                 se$newdata.xm_2sls,
                 mean_betas$data.xm_2sls_se,
                 xm_z1_2sls_detec$xm_z1_2sls_detec,
                 z1_coverage$z1_coverage,
                 s2$newdata.xm_2sls,
                 fmr_y_int$fmr_y_int
  )
  
  
  
  
  
  
  write.table(res,file=paste('500000_EXTRA_final_res_model_',model,'.txt',sep=''),sep='\t',row.names=FALSE)
  
  
  #how many times across repeat sims is an interaction detected in the incorrect direction? 
  n_z1problem[n_z1prob_track]=sum(z1problem)
  n_z1prob_track=n_z1prob_track+1
  
  
  
}

write(n_z1problem, file='z1problem.txt',append=FALSE, sep = "\n")



for (model in c(1:2)){

t50=data.frame(read.table(file=paste('500000_EXTRA_final_res_model_',model,'.txt',sep=''),header=TRUE))

all=t50

headers=c('mediator_coeff','\t','interac_coeff', '\t', 'mean_est','\t','sd(est)','\t','mean(se(est))','\t','se(est)','\t',
          'power','\t','type_i','\t','coverage')

res_l_0=all[round(all$mean_betas.data.xm_coeff_y,3)==0.000,]
res_l_m3=all[round(all$mean_betas.data.xm_coeff_y,3)==-0.111,]
res_l_3=all[round(all$mean_betas.data.xm_coeff_y,3)==0.111,]
res_l_5=all[round(all$mean_betas.data.xm_coeff_y,3)==0.167,]
res_l_1=all[round(all$mean_betas.data.xm_coeff_y,3)==0.333,]

res_l_0=res_l_0[order(res_l_0$mean_betas.data.x_coeff_m),]
res_l_m3=res_l_m3[order(res_l_m3$mean_betas.data.x_coeff_m),]
res_l_3=res_l_3[order(res_l_3$mean_betas.data.x_coeff_m),]
res_l_5=res_l_5[order(res_l_5$mean_betas.data.x_coeff_m),]
res_l_1=res_l_1[order(res_l_1$mean_betas.data.x_coeff_m),]

blank=c(1:5)
for (i in c(1:5)){blank[i]='NA'}

##################################INTERACTION COEFFICIENT#####################################################################################################################

#REMEMBER THAT THE VARIANCE NEEDS TO BE SQRT'D TO CONVERT TO SD
#POWER, TYPE I AND COVERAGE NEED TO BE DIVIDED BY 10 TO CONVERT TO %

###################
#Z=Z1+Z2+Z1Z2+Z1Z1#
###################

editZ1_res_l_0=data.frame(res_l_0$mean_betas.data.x_coeff_m,res_l_0$mean_betas.data.xm_coeff_y,res_l_0$mean_betas.data.xm_2sls,
                          sqrt(res_l_0$s2.newdata.xm_2sls),res_l_0$mean_betas.data.xm_2sls_se,res_l_0$se.newdata.xm_2sls,blank,
                          (res_l_0$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,(res_l_0$z1_coverage.z1_coverage)/10)

editZ1_res_l_m3=data.frame(res_l_m3$mean_betas.data.x_coeff_m,res_l_m3$mean_betas.data.xm_coeff_y,res_l_m3$mean_betas.data.xm_2sls,
                           sqrt(res_l_m3$s2.newdata.xm_2sls),res_l_m3$mean_betas.data.xm_2sls_se,res_l_m3$se.newdata.xm_2sls,(res_l_m3$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                           blank, (res_l_m3$z1_coverage.z1_coverage)/10)

editZ1_res_l_3=data.frame(res_l_3$mean_betas.data.x_coeff_m,res_l_3$mean_betas.data.xm_coeff_y,res_l_3$mean_betas.data.xm_2sls,
                          sqrt(res_l_3$s2.newdata.xm_2sls),res_l_3$mean_betas.data.xm_2sls_se,res_l_3$se.newdata.xm_2sls,(res_l_3$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                          blank, (res_l_3$z1_coverage.z1_coverage)/10)

editZ1_res_l_5=data.frame(res_l_5$mean_betas.data.x_coeff_m,res_l_5$mean_betas.data.xm_coeff_y,res_l_5$mean_betas.data.xm_2sls,
                          sqrt(res_l_5$s2.newdata.xm_2sls),res_l_5$mean_betas.data.xm_2sls_se,res_l_5$se.newdata.xm_2sls,(res_l_5$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                          blank, (res_l_5$z1_coverage.z1_coverage)/10)

editZ1_res_l_1=data.frame(res_l_1$mean_betas.data.x_coeff_m,res_l_1$mean_betas.data.xm_coeff_y,res_l_1$mean_betas.data.xm_2sls,
                          sqrt(res_l_1$s2.newdata.xm_2sls),res_l_1$mean_betas.data.xm_2sls_se,res_l_1$se.newdata.xm_2sls,(res_l_1$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                          blank, (res_l_1$z1_coverage.z1_coverage)/10)


#interaction coefficient=0
write.table(headers, file=paste('TSLS_MED_L0_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('TSLS_MED_L0_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editZ1_res_l_0, file=paste('TSLS_MED_L0_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=m3
write.table(headers, file=paste('TSLS_MED_LM3_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('TSLS_MED_LM3_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editZ1_res_l_m3, file=paste('TSLS_MED_LM3_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=3
write.table(headers, file=paste('TSLS_MED_L3_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('TSLS_MED_L3_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editZ1_res_l_3, file=paste('TSLS_MED_L3_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=5
write.table(headers, file=paste('TSLS_MED_L5_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('TSLS_MED_L5_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editZ1_res_l_5, file=paste('TSLS_MED_L5_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=1
write.table(headers, file=paste('TSLS_MED_L1_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('TSLS_MED_L1_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editZ1_res_l_1, file=paste('TSLS_MED_L1_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)



#####
#FMR#
#####

headers2=c('mediator_coeff','\t','interac_coeff', '\t', 'power','\t','type_i')


editfmr_res_l_0=data.frame(res_l_0$mean_betas.data.x_coeff_m,res_l_0$mean_betas.data.xm_coeff_y,blank,(res_l_0$fmr_y_int.fmr_y_int)/10)
editfmr_res_l_m3=data.frame(res_l_m3$mean_betas.data.x_coeff_m,res_l_m3$mean_betas.data.xm_coeff_y,(res_l_m3$fmr_y_int.fmr_y_int)/10,blank)
editfmr_res_l_3=data.frame(res_l_3$mean_betas.data.x_coeff_m,res_l_3$mean_betas.data.xm_coeff_y,(res_l_3$fmr_y_int.fmr_y_int)/10,blank)
editfmr_res_l_5=data.frame(res_l_5$mean_betas.data.x_coeff_m,res_l_5$mean_betas.data.xm_coeff_y,(res_l_5$fmr_y_int.fmr_y_int)/10,blank)
editfmr_res_l_1=data.frame(res_l_1$mean_betas.data.x_coeff_m,res_l_1$mean_betas.data.xm_coeff_y,(res_l_1$fmr_y_int.fmr_y_int)/10,blank)


#interaction coefficient=0
write.table(headers2, file=paste('fmr_L0_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('fmr_L0_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editfmr_res_l_0, file=paste('fmr_L0_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=m3
write.table(headers2, file=paste('fmr_LM3_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('fmr_LM3_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editfmr_res_l_m3, file=paste('fmr_LM3_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=3
write.table(headers2, file=paste('fmr_L3_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('fmr_L3_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editfmr_res_l_3, file=paste('fmr_L3_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=5
write.table(headers2, file=paste('fmr_L5_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('fmr_L5_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editfmr_res_l_5, file=paste('fmr_L5_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=1
write.table(headers2, file=paste('fmr_L1_',model,'_.txt'),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file=paste('fmr_L1_',model,'_.txt'),append=TRUE, sep = "\n")
write.table(editfmr_res_l_1, file=paste('fmr_L1_',model,'_.txt'),append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)

}
