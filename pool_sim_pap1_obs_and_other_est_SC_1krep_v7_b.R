#script_name: pool_sim_pap1_obs_and_other_est_SC_1krep_v7_b.r
#project: 4-way decomp: paper 1
#script author: Teri North
#script purpose: pool estimates across simulation repeats by
#                 -taking the mean betahat & SE of betahats (to generate MC 95% CI for betahat) 
#                 -take the mean SE and the SD of betahats
#                 -calculate power, type i error and coverage where applicable
#date created: 03/01/2019
#last edited: 17/01/2019
#notes:

setwd('') #Folder 1

#number of repeats in each sim
repeats=100

n_z1problem=c(1:2)
for (j in c(1:2)){n_z1problem[j]=0}
n_z1prob_track=1

for (nval in c(50000,500000)){
  
  xm_z1_2sls_detec=c(1:25)
  for (i in c(1:25)){
    xm_z1_2sls_detec[i]=0
  }
  xm_obs_detec=c(1:25) # counter for # times interaction is detected (null is rejected) at 5% level
  for (i in c(1:25)){
    xm_obs_detec[i]=0
  }  
  
  
  z1_coverage=c(1:25)
  for (i in c(1:25)){
    z1_coverage[i]=0
  }
  obscoverage=c(1:25)
  for (i in c(1:25)){
    obscoverage[i]=0
  }
  
  #reality check
  #how many times is the interaction detected (p<0.05), but the estimate is in the opposite direction to true effect?
  z1problem=c(1:25) 
  for (i in c(1:25)){
    z1problem[i]=0
  }


  #calculating the mean betas   
  first=1
  

  for (seedval in c(7821897,8376154,649384402,238140535,170379645,312006101,713795870,169378934,456561608,28335714)){
    
    for (rep in c(1:repeats)){
      
      if (first==1){
        
        data=read.table(file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_ROBUST_res.txt",sep=''),sep='\t',header=TRUE)
        true_vals=data.frame(data$x_coeff_m,data$x_coeff_y,data$m_coeff_y, data$xm_coeff_y)
        first=0
        
        ll=data.frame(
          xm_z1_2sls_ll=data$xm_z1_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_z1_2sls_se,
          xm_obs_ll=data$xm_obs-(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_obs_se
        )
        
        ul=data.frame(
          xm_z1_2sls_ul=data$xm_z1_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_z1_2sls_se,
          xm_obs_ul=data$xm_obs+(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_obs_se
        )
        
        
        for (i in c(1:25)){
          if (ll$xm_z1_2sls_ll[i]>0 | ul$xm_z1_2sls_ul[i]<0){xm_z1_2sls_detec[i]=1}
        } 
        
        
        for (i in c(1:25)){
          if (ll$xm_obs_ll[i]>0 | ul$xm_obs_ul[i]<0){xm_obs_detec[i]=1}
        }
        
        
        for (i in c(1:25)){
          if ((ll$xm_z1_2sls_ll[i]<data$xm_coeff_y[i]) & (ul$xm_z1_2sls_ul[i]>data$xm_coeff_y[i])){z1_coverage[i]=1}
        }
        

        for (i in c(1:25)){
          if ((ll$xm_obs_ll[i]<data$xm_coeff_y[i]) & (ul$xm_obs_ul[i]>data$xm_coeff_y[i])){obscoverage[i]=1}
        }
        

        for (i in c(1:25)){
          if (((ll$xm_z1_2sls_ll[i]>0) & (data$xm_coeff_y[i]<0))|((ul$xm_z1_2sls_ul[i]<0) & (data$xm_coeff_y[i]>0))) {z1problem[i]=1}#if interac detec, but coeff wrong direc
        }
        


        
      } else if (first==0){
        
        new=read.table(file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_ROBUST_res.txt",sep=''),sep='\t',header=TRUE)
        data=data+new
        
        ll_new=data.frame(
          xm_z1_2sls_ll=new$xm_z1_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_z1_2sls_se,
          xm_obs_ll=new$xm_obs-(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_obs_se
        )
        
        ul_new=data.frame(
          xm_z1_2sls_ul=new$xm_z1_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_z1_2sls_se,
          xm_obs_ul=new$xm_obs+(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_obs_se
        )
        
        

        for (i in c(1:25)){
          if (ll_new$xm_z1_2sls_ll[i]>0 | ul_new$xm_z1_2sls_ul[i]<0){xm_z1_2sls_detec[i]=xm_z1_2sls_detec[i]+1}
        }
        

        for (i in c(1:25)){
          if (ll_new$xm_obs_ll[i]>0 | ul_new$xm_obs_ul[i]<0){xm_obs_detec[i]=xm_obs_detec[i]+1}
        }
        

        
        for (i in c(1:25)){
          if ((ll_new$xm_z1_2sls_ll[i]<new$xm_coeff_y[i]) & (ul_new$xm_z1_2sls_ul[i]>new$xm_coeff_y[i])){z1_coverage[i]=z1_coverage[i]+1}
        }
        

        for (i in c(1:25)){
          if ((ll_new$xm_obs_ll[i]<new$xm_coeff_y[i]) & (ul_new$xm_obs_ul[i]>new$xm_coeff_y[i])){obscoverage[i]=obscoverage[i]+1}
        }
        
        
        

        for (i in c(1:25)){
          if (((ll_new$xm_z1_2sls_ll[i]>0) & (new$xm_coeff_y[i]<0))|((ul_new$xm_z1_2sls_ul[i]<0) & (new$xm_coeff_y[i]>0))) {z1problem[i]=z1problem[i]+1}#if interac detec, but coeff wrong direc
        } 
        


        
      }
      
      
      
    }
    
  }
  
  
  
  #remove true values
  data_est=data.frame(data$x_obs,data$x_obs_se,
                      data$m_obs,data$m_obs_se,
                      data$xm_obs,data$xm_obs_se,
                      data$x_2sls,data$x_2sls_se,
                      data$x_z1_2sls,data$x_z1_2sls_se,
                      data$m_2sls,data$m_2sls_se,
                      data$m_z1_2sls,data$m_z1_2sls_se,
                      data$xm_2sls,data$xm_2sls_se,
                      data$xm_z1_2sls,data$xm_z1_2sls_se
  )
  
  
  
  #mean betas
  mean_denom=repeats*10 #no. rep within seeds * no. seeds
  data_mean=data_est/mean_denom #gives mean beta and mean se
  #add in the true params
  mean_betas=cbind(true_vals,data_mean)
  
  
  #############################################################################################################################################################################
  
  #now for the standard error
  checker=1
  
  for (seeds in c(7821897,8376154,649384402,238140535,170379645,312006101,713795870,169378934,456561608,28335714)){
    
    for (rep in c(1:repeats)){
      
      if (checker==1){
        
        newdata=read.table(file=paste(seeds,'_rep',rep,'_samp',nval,"_FMR_ROBUST_res.txt",sep=''),sep='\t',header=TRUE)
        newdata=(newdata-(data/mean_denom))^2
        checker=0
        
      } else if (checker==0){
        
        
        newer=read.table(file=paste(seeds,'_rep',rep,'_samp',nval,"_FMR_ROBUST_res.txt",sep=''),sep='\t',header=TRUE)
        newdata=newdata+(newer-(data/mean_denom))^2
        newdata=data.frame(newdata)
        
        
      }
      
    }
  }
  
  
  newdata_est=data.frame(newdata$x_obs,
                         newdata$m_obs,
                         newdata$xm_obs,
                         newdata$x_2sls,
                         newdata$x_z1_2sls,
                         newdata$m_2sls,
                         newdata$m_z1_2sls,
                         newdata$xm_2sls,
                         newdata$xm_z1_2sls
  )
  
  
  
  #divide by n-1 to get s^2
  s2=newdata_est/(repeats*10-1)
  se=sqrt(s2/(repeats*10))
  
  
  ##################################################################################################################################################################
  
  xm_z1_2sls_detec=data.frame(xm_z1_2sls_detec)
  xm_obs_detec=data.frame(xm_obs_detec)
  z1_coverage=data.frame(z1_coverage)
  obscoverage=data.frame(obscoverage)

  #results table 
  res=data.frame(mean_betas$data.x_coeff_m,
                 mean_betas$data.x_coeff_y,
                 mean_betas$data.m_coeff_y,
                 mean_betas$data.xm_coeff_y,
                 mean_betas$data.xm_obs,
                 se$newdata.xm_obs,
                 mean_betas$data.xm_z1_2sls,
                 se$newdata.xm_z1_2sls,
                 mean_betas$data.xm_obs_se,
                 mean_betas$data.xm_z1_2sls_se,
                 xm_obs_detec$xm_obs_detec,
                 xm_z1_2sls_detec$xm_z1_2sls_detec,
                 z1_coverage$z1_coverage,
                 obscoverage$obscoverage,
                 s2$newdata.xm_obs,
                 s2$newdata.xm_z1_2sls)
  
  
  
  
  
  
  write.table(res,file=paste(nval,'_EXTRA_final_ROBUST_res.txt',sep=''),sep='\t',row.names=FALSE)
  
  
  #how many times across all models and repeat sims is an interaction detected in the incorrect direction? 
  n_z1problem[n_z1prob_track]=sum(z1problem)
  n_z1prob_track=n_z1prob_track+1
  

  
}

write(n_z1problem, file='z1_ROBUST_problem.txt',append=FALSE, sep = "\n")


#Now read back in so that we have all the data across all sample sizes


n_t5=c(1:25)
for (i in c(1:25)){
  n_t5[i]=50000
}
n_t50=c(1:25)
for (i in c(1:25)){
  n_t50[i]=500000
}



sampsize=n_t5
t5=data.frame(sampsize,read.table(file=paste(50000,'_EXTRA_final_ROBUST_res.txt',sep=''),header=TRUE))
sampsize=n_t50
t50=data.frame(sampsize,read.table(file=paste(500000,'_EXTRA_final_ROBUST_res.txt',sep=''),header=TRUE))

all=rbind(t5,t50)

headers=c('mediator_coeff','\t','interac_coeff', '\t', 'sample_size','\t','mean_est','\t',
          'sd(est)','\t','mean(se(est))','\t','se(est)','\t',
          'power/type_i','\t','coverage')

res_N_50=all[all$sampsize==50000,]
res_N_500=all[all$sampsize==500000,]

res_N_50=res_N_50[order(res_N_50$mean_betas.data.x_coeff_m,res_N_50$mean_betas.data.xm_coeff_y),]
res_N_500=res_N_500[order(res_N_500$mean_betas.data.x_coeff_m,res_N_500$mean_betas.data.xm_coeff_y),]

##################################INTERACTION COEFFICIENT#####################################################################################################################

#REMEMBER THAT THE VARIANCE NEEDS TO BE SQRT'D TO CONVERT TO SD
#POWER, TYPE I AND COVERAGE NEED TO BE DIVIDED BY 10 TO CONVERT TO %

###################
#Z=Z1+Z2+Z1Z2+Z1Z1#
###################

editZ1_res_N_50=data.frame(res_N_50$mean_betas.data.x_coeff_m,res_N_50$mean_betas.data.xm_coeff_y,res_N_50$sampsize,res_N_50$mean_betas.data.xm_z1_2sls,
                          sqrt(res_N_50$s2.newdata.xm_z1_2sls),res_N_50$mean_betas.data.xm_z1_2sls_se,res_N_50$se.newdata.xm_z1_2sls, 
                          (res_N_50$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,(res_N_50$z1_coverage.z1_coverage)/10)
editZ1_res_N_500=data.frame(res_N_500$mean_betas.data.x_coeff_m,res_N_500$mean_betas.data.xm_coeff_y,res_N_500$sampsize,res_N_500$mean_betas.data.xm_z1_2sls,
                           sqrt(res_N_500$s2.newdata.xm_z1_2sls),res_N_500$mean_betas.data.xm_z1_2sls_se,res_N_500$se.newdata.xm_z1_2sls,(res_N_500$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                           (res_N_500$z1_coverage.z1_coverage)/10)



#sample size=50k
write.table(headers, file='TSLS_MED_SS50K_ROBSE.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_MED_SS50K_ROBSE.txt',append=TRUE, sep = "\n")
write.table(editZ1_res_N_50, file='TSLS_MED_SS50K_ROBSE.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#sample size=500k
write.table(headers, file='TSLS_MED_SS500K_ROBSE.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_MED_SS500K_ROBSE.txt',append=TRUE, sep = "\n")
write.table(editZ1_res_N_500, file='TSLS_MED_SS500K_ROBSE.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)



###############
#OBSERVATIONAL#
###############

edit_obs_res_N_50=data.frame(res_N_50$mean_betas.data.x_coeff_m,res_N_50$mean_betas.data.xm_coeff_y,res_N_50$sampsize,res_N_50$mean_betas.data.xm_obs,
                           sqrt(res_N_50$s2.newdata.xm_obs),res_N_50$mean_betas.data.xm_obs_se,res_N_50$se.newdata.xm_obs,
                           (res_N_50$xm_obs_detec.xm_obs_detec)/10,(res_N_50$obscoverage.obscoverage)/10)

edit_obs_res_N_500=data.frame(res_N_500$mean_betas.data.x_coeff_m,res_N_500$mean_betas.data.xm_coeff_y,res_N_500$sampsize,res_N_500$mean_betas.data.xm_obs,
                            sqrt(res_N_500$s2.newdata.xm_obs),res_N_500$mean_betas.data.xm_obs_se,res_N_500$se.newdata.xm_obs,(res_N_500$xm_obs_detec.xm_obs_detec)/10,
                            (res_N_500$obscoverage.obscoverage)/10)

#sample size=50k
write.table(headers, file='OBS_SS50K_ROBSE.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='OBS_SS50K_ROBSE.txt',append=TRUE, sep = "\n")
write.table(edit_obs_res_N_50, file='OBS_SS50K_ROBSE.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#sample size=500k
write.table(headers, file='OBS_SS500K_ROBSE.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='OBS_SS500K_ROBSE.txt',append=TRUE, sep = "\n")
write.table(edit_obs_res_N_500, file='OBS_SS500K_ROBSE.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)

