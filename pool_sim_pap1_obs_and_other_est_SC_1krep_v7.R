#script_name: pool_sim_pap1_obs_and_other_est_SC_1krep_v7.r
#project: 4-way decomp: paper 1
#script author: Teri North
#script purpose: pool estimates across simulation repeats by
#                 -taking the mean betahat & SE of betahats (to generate MC 95% CI for betahat) 
#                 -take the mean SE and the SD of betahats
#                 -calculate power, type i error and coverage where applicable
#date created: 03/01/2019
#last edited: 03/01/2019
#notes:

setwd('') #Folder 1

#number of repeats in each sim
repeats=100

#tracker for erroneous calls
n_problem=c(1:12)
for (j in c(1:12)){n_problem[j]=0}
n_prob_track=1

n_z1problem=c(1:12)
for (j in c(1:12)){n_z1problem[j]=0}
n_z1prob_track=1

for (nval in c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,500000,1000000)){
  
  xm_2sls_detec=c(1:25) # counter for # times interaction is detected (null is rejected) at 5% level
  for (i in c(1:25)){
    xm_2sls_detec[i]=0
  }
  xm_z1_2sls_detec=c(1:25)
  for (i in c(1:25)){
    xm_z1_2sls_detec[i]=0
  }
  xm_obs_detec=c(1:25) # counter for # times interaction is detected (null is rejected) at 5% level
  for (i in c(1:25)){
    xm_obs_detec[i]=0
  }  
  
  
  coverage=c(1:25) # counter for # times 95% CI contains true value
  for (i in c(1:25)){
    coverage[i]=0
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
  problem=c(1:25) 
  for (i in c(1:25)){
    problem[i]=0
  }
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
  

  for (seedval in c(7821897,8376154,649384402,238140535,170379645,312006101,713795870,169378934,456561608,28335714)){
    
    for (rep in c(1:repeats)){
      
      if (first==1){
        
        data=read.table(file=paste(seedval,'_rep',rep,'_samp',nval,"_FMRres.txt",sep=''),sep='\t',header=TRUE)
        true_vals=data.frame(data$x_coeff_m,data$x_coeff_y,data$m_coeff_y, data$xm_coeff_y)
        first=0
        
        ll=data.frame(
          xm_2sls_ll=data$xm_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_2sls_se,
          xm_z1_2sls_ll=data$xm_z1_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_z1_2sls_se,
          xm_obs_ll=data$xm_obs-(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_obs_se
        )
        
        ul=data.frame(
          xm_2sls_ul=data$xm_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_2sls_se,
          xm_z1_2sls_ul=data$xm_z1_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_z1_2sls_se,
          xm_obs_ul=data$xm_obs+(qt(0.025,(nval-4),lower.tail=FALSE))*data$xm_obs_se
        )
        
        
        for (i in c(1:25)){
          if (is.na(ll$xm_2sls_ll[i])|is.na(ul$xm_2sls_ul[i])){
            xm_2sls_detec[i]=NA
          } else if (ll$xm_2sls_ll[i]>0 | ul$xm_2sls_ul[i]<0){xm_2sls_detec[i]=1}
        }
        
        for (i in c(1:25)){
          if (ll$xm_z1_2sls_ll[i]>0 | ul$xm_z1_2sls_ul[i]<0){xm_z1_2sls_detec[i]=1}
        } 
        
        
        for (i in c(1:25)){
          if (ll$xm_obs_ll[i]>0 | ul$xm_obs_ul[i]<0){xm_obs_detec[i]=1}
        }
        
        
        for (i in c(1:25)){
          if (is.na(ll$xm_2sls_ll[i])|is.na(ul$xm_2sls_ul[i])){
            coverage[i]=NA
          } else if ((ll$xm_2sls_ll[i]<data$xm_coeff_y[i]) & (ul$xm_2sls_ul[i]>data$xm_coeff_y[i])){coverage[i]=1}
        }  
        
        for (i in c(1:25)){
          if ((ll$xm_z1_2sls_ll[i]<data$xm_coeff_y[i]) & (ul$xm_z1_2sls_ul[i]>data$xm_coeff_y[i])){z1_coverage[i]=1}
        }
        

        for (i in c(1:25)){
          if ((ll$xm_obs_ll[i]<data$xm_coeff_y[i]) & (ul$xm_obs_ul[i]>data$xm_coeff_y[i])){obscoverage[i]=1}
        }
        
        
        for (i in c(1:25)){
          if (is.na(ll$xm_2sls_ll[i])|is.na(ul$xm_2sls_ul[i])){
            problem[i]=NA
          } else if (((ll$xm_2sls_ll[i]>0) & (data$xm_coeff_y[i]<0))|((ul$xm_2sls_ul[i]<0) & (data$xm_coeff_y[i]>0))) {problem[i]=1}#if interac detec, but coeff wrong direc
        } 
        
        for (i in c(1:25)){
          if (((ll$xm_z1_2sls_ll[i]>0) & (data$xm_coeff_y[i]<0))|((ul$xm_z1_2sls_ul[i]<0) & (data$xm_coeff_y[i]>0))) {z1problem[i]=1}#if interac detec, but coeff wrong direc
        }
        

        for (i in c(1:25)){
          if (data$fmr_interac_p[i]<0.05){fmr_y_int[i]=1}
        }
        
        
        
      } else if (first==0){
        
        new=read.table(file=paste(seedval,'_rep',rep,'_samp',nval,"_FMRres.txt",sep=''),sep='\t',header=TRUE)
        data=data+new
        
        ll_new=data.frame(
          xm_2sls_ll=new$xm_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_2sls_se,
          xm_z1_2sls_ll=new$xm_z1_2sls-(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_z1_2sls_se,
          xm_obs_ll=new$xm_obs-(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_obs_se
        )
        
        ul_new=data.frame(
          xm_2sls_ul=new$xm_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_2sls_se,
          xm_z1_2sls_ul=new$xm_z1_2sls+(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_z1_2sls_se,
          xm_obs_ul=new$xm_obs+(qt(0.025,(nval-4),lower.tail=FALSE))*new$xm_obs_se
        )
        
        
        for (i in c(1:25)){
          if (is.na(ll_new$xm_2sls_ll[i])|is.na(ul_new$xm_2sls_ul[i])){
            xm_2sls_detec[i]=NA
          } else if (ll_new$xm_2sls_ll[i]>0 | ul_new$xm_2sls_ul[i]<0){xm_2sls_detec[i]=xm_2sls_detec[i]+1}
        }
        
        for (i in c(1:25)){
          if (ll_new$xm_z1_2sls_ll[i]>0 | ul_new$xm_z1_2sls_ul[i]<0){xm_z1_2sls_detec[i]=xm_z1_2sls_detec[i]+1}
        }
        

        for (i in c(1:25)){
          if (ll_new$xm_obs_ll[i]>0 | ul_new$xm_obs_ul[i]<0){xm_obs_detec[i]=xm_obs_detec[i]+1}
        }
        
        
        
        for (i in c(1:25)){
          if (is.na(ll_new$xm_2sls_ll[i])|is.na(ul_new$xm_2sls_ul[i])){
            coverage[i]=NA
          } else if ((ll_new$xm_2sls_ll[i]<new$xm_coeff_y[i]) & (ul_new$xm_2sls_ul[i]>new$xm_coeff_y[i])){coverage[i]=coverage[i]+1}
        }  
        
        for (i in c(1:25)){
          if ((ll_new$xm_z1_2sls_ll[i]<new$xm_coeff_y[i]) & (ul_new$xm_z1_2sls_ul[i]>new$xm_coeff_y[i])){z1_coverage[i]=z1_coverage[i]+1}
        }
        

        for (i in c(1:25)){
          if ((ll_new$xm_obs_ll[i]<new$xm_coeff_y[i]) & (ul_new$xm_obs_ul[i]>new$xm_coeff_y[i])){obscoverage[i]=obscoverage[i]+1}
        }
        
        
        
        
        for (i in c(1:25)){
          if (is.na(ll_new$xm_2sls_ll[i])|is.na(ul_new$xm_2sls_ul[i])){
            problem[i]=NA
          } else if (((ll_new$xm_2sls_ll[i]>0) & (new$xm_coeff_y[i]<0))|((ul_new$xm_2sls_ul[i]<0) & (new$xm_coeff_y[i]>0))) {problem[i]=problem[i]+1}#if interac detec, but coeff wrong direc
        } 
        
        for (i in c(1:25)){
          if (((ll_new$xm_z1_2sls_ll[i]>0) & (new$xm_coeff_y[i]<0))|((ul_new$xm_z1_2sls_ul[i]<0) & (new$xm_coeff_y[i]>0))) {z1problem[i]=z1problem[i]+1}#if interac detec, but coeff wrong direc
        } 
        

        for (i in c(1:25)){
          if (new$fmr_interac_p[i]<0.05){fmr_y_int[i]=fmr_y_int[i]+1}
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
        
        newdata=read.table(file=paste(seeds,'_rep',rep,'_samp',nval,"_FMRres.txt",sep=''),sep='\t',header=TRUE)
        newdata=(newdata-(data/mean_denom))^2
        checker=0
        
      } else if (checker==0){
        
        
        newer=read.table(file=paste(seeds,'_rep',rep,'_samp',nval,"_FMRres.txt",sep=''),sep='\t',header=TRUE)
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
  
  xm_2sls_detec=data.frame(xm_2sls_detec)
  xm_z1_2sls_detec=data.frame(xm_z1_2sls_detec)
  xm_obs_detec=data.frame(xm_obs_detec)
  coverage=data.frame(coverage)
  z1_coverage=data.frame(z1_coverage)
  obscoverage=data.frame(obscoverage)
  fmr_y_int=data.frame(fmr_y_int)
  
  #results table 
  res=data.frame(mean_betas$data.x_coeff_m,
                 mean_betas$data.x_coeff_y,
                 mean_betas$data.m_coeff_y,
                 mean_betas$data.xm_coeff_y,
                 mean_betas$data.x_obs,
                 se$newdata.x_obs,
                 mean_betas$data.m_obs,
                 se$newdata.m_obs,
                 mean_betas$data.xm_obs,
                 se$newdata.xm_obs,
                 mean_betas$data.x_2sls,
                 se$newdata.x_2sls,
                 mean_betas$data.x_z1_2sls,
                 se$newdata.x_z1_2sls,
                 mean_betas$data.m_2sls,
                 se$newdata.m_2sls,
                 mean_betas$data.m_z1_2sls,
                 se$newdata.m_z1_2sls,
                 mean_betas$data.xm_2sls,
                 se$newdata.xm_2sls,
                 mean_betas$data.xm_z1_2sls,
                 se$newdata.xm_z1_2sls,
                 mean_betas$data.xm_obs_se,
                 mean_betas$data.xm_2sls_se,
                 mean_betas$data.xm_z1_2sls_se,
                 xm_obs_detec$xm_obs_detec,
                 xm_2sls_detec$xm_2sls_detec,
                 xm_z1_2sls_detec$xm_z1_2sls_detec,
                 coverage$coverage,
                 z1_coverage$z1_coverage,
                 obscoverage$obscoverage,
                 s2$newdata.xm_obs,
                 s2$newdata.xm_2sls,
                 s2$newdata.xm_z1_2sls,
                 fmr_y_int$fmr_y_int
  )
  
  
  
  
  
  
  write.table(res,file=paste(nval,'_EXTRA_final_res.txt',sep=''),sep='\t',row.names=FALSE)
  
  
  #how many times across all models and repeat sims is an interaction detected in the incorrect direction? 
  n_problem[n_prob_track]=sum(problem)
  n_prob_track=n_prob_track+1
  
  n_z1problem[n_z1prob_track]=sum(z1problem)
  n_z1prob_track=n_z1prob_track+1
  

  
}

write(n_problem, file='problem.txt',append=FALSE, sep = "\n")
write(n_z1problem, file='z1problem.txt',append=FALSE, sep = "\n")


#Now read back in so that we have all the data across all sample sizes

n_t1=c(1:25)
for (i in c(1:25)){
  n_t1[i]=10000
}
n_t2=c(1:25)
for (i in c(1:25)){
  n_t2[i]=20000
}
n_t3=c(1:25)
for (i in c(1:25)){
  n_t3[i]=30000
}
n_t4=c(1:25)
for (i in c(1:25)){
  n_t4[i]=40000
}
n_t5=c(1:25)
for (i in c(1:25)){
  n_t5[i]=50000
}
n_t6=c(1:25)
for (i in c(1:25)){
  n_t6[i]=60000
}
n_t7=c(1:25)
for (i in c(1:25)){
  n_t7[i]=70000
}
n_t8=c(1:25)
for (i in c(1:25)){
  n_t8[i]=80000
}
n_t9=c(1:25)
for (i in c(1:25)){
  n_t9[i]=90000
}
n_t10=c(1:25)
for (i in c(1:25)){
  n_t10[i]=100000
}
n_t50=c(1:25)
for (i in c(1:25)){
  n_t50[i]=500000
}
n_t100=c(1:25)
for (i in c(1:25)){
  n_t100[i]=1000000
}


sampsize=n_t1
t1=data.frame(sampsize,read.table(file=paste(10000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t2
t2=data.frame(sampsize,read.table(file=paste(20000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t3
t3=data.frame(sampsize,read.table(file=paste(30000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t4
t4=data.frame(sampsize,read.table(file=paste(40000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t5
t5=data.frame(sampsize,read.table(file=paste(50000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t6
t6=data.frame(sampsize,read.table(file=paste(60000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t7
t7=data.frame(sampsize,read.table(file=paste(70000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t8
t8=data.frame(sampsize,read.table(file=paste(80000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t9
t9=data.frame(sampsize,read.table(file=paste(90000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t10
t10=data.frame(sampsize,read.table(file=paste(100000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t50
t50=data.frame(sampsize,read.table(file=paste(500000,'_EXTRA_final_res.txt',sep=''),header=TRUE))
sampsize=n_t100
t100=data.frame(sampsize,read.table(file=paste(1000000,'_EXTRA_final_res.txt',sep=''),header=TRUE))

all=rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t50,t100)

headers=c('mediator_coeff','\t','interac_coeff', '\t', 'sample_size','\t','mean_est','\t',
          'sd(est)','\t','mean(se(est))','\t','se(est)','\t',
          'power','\t','type_i','\t','coverage')

res_l_0=all[round(all$mean_betas.data.xm_coeff_y,3)==0.000,]
res_l_m3=all[round(all$mean_betas.data.xm_coeff_y,3)==-0.111,]
res_l_3=all[round(all$mean_betas.data.xm_coeff_y,3)==0.111,]
res_l_5=all[round(all$mean_betas.data.xm_coeff_y,3)==0.167,]
res_l_1=all[round(all$mean_betas.data.xm_coeff_y,3)==0.333,]

res_l_0=res_l_0[order(res_l_0$mean_betas.data.x_coeff_m,res_l_0$sampsize),]
res_l_m3=res_l_m3[order(res_l_m3$mean_betas.data.x_coeff_m,res_l_m3$sampsize),]
res_l_3=res_l_3[order(res_l_3$mean_betas.data.x_coeff_m,res_l_3$sampsize),]
res_l_5=res_l_5[order(res_l_5$mean_betas.data.x_coeff_m,res_l_5$sampsize),]
res_l_1=res_l_1[order(res_l_1$mean_betas.data.x_coeff_m,res_l_1$sampsize),]

blank=c(1:60)
for (i in c(1:60)){blank[i]='NA'}

##################################INTERACTION COEFFICIENT#####################################################################################################################

#REMEMBER THAT THE VARIANCE NEEDS TO BE SQRT'D TO CONVERT TO SD
#POWER, TYPE I AND COVERAGE NEED TO BE DIVIDED BY 10 TO CONVERT TO %
##############
#Z=Z1+Z2+Z1Z2#
##############
edit_res_l_0=data.frame(res_l_0$mean_betas.data.x_coeff_m,res_l_0$mean_betas.data.xm_coeff_y,res_l_0$sampsize,res_l_0$mean_betas.data.xm_2sls,
                        sqrt(res_l_0$s2.newdata.xm_2sls),res_l_0$mean_betas.data.xm_2sls_se,res_l_0$se.newdata.xm_2sls,blank,
                        (res_l_0$xm_2sls_detec.xm_2sls_detec)/10,(res_l_0$coverage.coverage)/10)
edit_res_l_m3=data.frame(res_l_m3$mean_betas.data.x_coeff_m,res_l_m3$mean_betas.data.xm_coeff_y,res_l_m3$sampsize,res_l_m3$mean_betas.data.xm_2sls,
                         sqrt(res_l_m3$s2.newdata.xm_2sls),res_l_m3$mean_betas.data.xm_2sls_se, res_l_m3$se.newdata.xm_2sls,(res_l_m3$xm_2sls_detec.xm_2sls_detec)/10,
                         blank, (res_l_m3$coverage.coverage)/10)
edit_res_l_3=data.frame(res_l_3$mean_betas.data.x_coeff_m,res_l_3$mean_betas.data.xm_coeff_y,res_l_3$sampsize,res_l_3$mean_betas.data.xm_2sls,
                        sqrt(res_l_3$s2.newdata.xm_2sls),res_l_3$mean_betas.data.xm_2sls_se,res_l_3$se.newdata.xm_2sls,(res_l_3$xm_2sls_detec.xm_2sls_detec)/10,
                        blank, (res_l_3$coverage.coverage)/10)
edit_res_l_5=data.frame(res_l_5$mean_betas.data.x_coeff_m,res_l_5$mean_betas.data.xm_coeff_y,res_l_5$sampsize,res_l_5$mean_betas.data.xm_2sls,
                        sqrt(res_l_5$s2.newdata.xm_2sls),res_l_5$mean_betas.data.xm_2sls_se,res_l_5$se.newdata.xm_2sls,(res_l_5$xm_2sls_detec.xm_2sls_detec)/10,
                        blank, (res_l_5$coverage.coverage)/10)
edit_res_l_1=data.frame(res_l_1$mean_betas.data.x_coeff_m,res_l_1$mean_betas.data.xm_coeff_y,res_l_1$sampsize,res_l_1$mean_betas.data.xm_2sls,
                        sqrt(res_l_1$s2.newdata.xm_2sls),res_l_1$mean_betas.data.xm_2sls_se,res_l_1$se.newdata.xm_2sls,(res_l_1$xm_2sls_detec.xm_2sls_detec)/10,
                        blank, (res_l_1$coverage.coverage)/10)


#interaction coefficient=0
write.table(headers, file='TSLS_NOMED_L0.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_NOMED_L0.txt',append=TRUE, sep = "\n")
write.table(edit_res_l_0, file='TSLS_NOMED_L0.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=m3
write.table(headers, file='TSLS_NOMED_LM3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_NOMED_LM3.txt',append=TRUE, sep = "\n")
write.table(edit_res_l_m3, file='TSLS_NOMED_LM3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=3
write.table(headers, file='TSLS_NOMED_L3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_NOMED_L3.txt',append=TRUE, sep = "\n")
write.table(edit_res_l_3, file='TSLS_NOMED_L3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=5
write.table(headers, file='TSLS_NOMED_L5.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_NOMED_L5.txt',append=TRUE, sep = "\n")
write.table(edit_res_l_5, file='TSLS_NOMED_L5.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=1
write.table(headers, file='TSLS_NOMED_L1.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_NOMED_L1.txt',append=TRUE, sep = "\n")
write.table(edit_res_l_1, file='TSLS_NOMED_L1.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)

###################
#Z=Z1+Z2+Z1Z2+Z1Z1#
###################

editZ1_res_l_0=data.frame(res_l_0$mean_betas.data.x_coeff_m,res_l_0$mean_betas.data.xm_coeff_y,res_l_0$sampsize,res_l_0$mean_betas.data.xm_z1_2sls,
                          sqrt(res_l_0$s2.newdata.xm_z1_2sls),res_l_0$mean_betas.data.xm_z1_2sls_se,res_l_0$se.newdata.xm_z1_2sls, blank,
                          (res_l_0$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,(res_l_0$z1_coverage.z1_coverage)/10)
editZ1_res_l_m3=data.frame(res_l_m3$mean_betas.data.x_coeff_m,res_l_m3$mean_betas.data.xm_coeff_y,res_l_m3$sampsize,res_l_m3$mean_betas.data.xm_z1_2sls,
                           sqrt(res_l_m3$s2.newdata.xm_z1_2sls),res_l_m3$mean_betas.data.xm_z1_2sls_se,res_l_m3$se.newdata.xm_z1_2sls,(res_l_m3$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                           blank, (res_l_m3$z1_coverage.z1_coverage)/10)
editZ1_res_l_3=data.frame(res_l_3$mean_betas.data.x_coeff_m,res_l_3$mean_betas.data.xm_coeff_y,res_l_3$sampsize,res_l_3$mean_betas.data.xm_z1_2sls,
                          sqrt(res_l_3$s2.newdata.xm_z1_2sls),res_l_3$mean_betas.data.xm_z1_2sls_se,res_l_3$se.newdata.xm_z1_2sls,(res_l_3$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                          blank, (res_l_3$z1_coverage.z1_coverage)/10)
editZ1_res_l_5=data.frame(res_l_5$mean_betas.data.x_coeff_m,res_l_5$mean_betas.data.xm_coeff_y,res_l_5$sampsize,res_l_5$mean_betas.data.xm_z1_2sls,
                          sqrt(res_l_5$s2.newdata.xm_z1_2sls),res_l_5$mean_betas.data.xm_z1_2sls_se,res_l_5$se.newdata.xm_z1_2sls,(res_l_5$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                          blank, (res_l_5$z1_coverage.z1_coverage)/10)
editZ1_res_l_1=data.frame(res_l_1$mean_betas.data.x_coeff_m,res_l_1$mean_betas.data.xm_coeff_y,res_l_1$sampsize,res_l_1$mean_betas.data.xm_z1_2sls,
                          sqrt(res_l_1$s2.newdata.xm_z1_2sls),res_l_1$mean_betas.data.xm_z1_2sls_se,res_l_1$se.newdata.xm_z1_2sls,(res_l_1$xm_z1_2sls_detec.xm_z1_2sls_detec)/10,
                          blank, (res_l_1$z1_coverage.z1_coverage)/10)


#interaction coefficient=0
write.table(headers, file='TSLS_MED_L0.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_MED_L0.txt',append=TRUE, sep = "\n")
write.table(editZ1_res_l_0, file='TSLS_MED_L0.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=m3
write.table(headers, file='TSLS_MED_LM3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_MED_LM3.txt',append=TRUE, sep = "\n")
write.table(editZ1_res_l_m3, file='TSLS_MED_LM3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=3
write.table(headers, file='TSLS_MED_L3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_MED_L3.txt',append=TRUE, sep = "\n")
write.table(editZ1_res_l_3, file='TSLS_MED_L3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=5
write.table(headers, file='TSLS_MED_L5.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_MED_L5.txt',append=TRUE, sep = "\n")
write.table(editZ1_res_l_5, file='TSLS_MED_L5.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=1
write.table(headers, file='TSLS_MED_L1.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='TSLS_MED_L1.txt',append=TRUE, sep = "\n")
write.table(editZ1_res_l_1, file='TSLS_MED_L1.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)



#####
#FMR#
#####

headers2=c('mediator_coeff','\t','interac_coeff', '\t', 'sample_size','\t','power','\t','type_i')


editfmr_res_l_0=data.frame(res_l_0$mean_betas.data.x_coeff_m,res_l_0$mean_betas.data.xm_coeff_y,res_l_0$sampsize,blank,(res_l_0$fmr_y_int.fmr_y_int)/10)
editfmr_res_l_m3=data.frame(res_l_m3$mean_betas.data.x_coeff_m,res_l_m3$mean_betas.data.xm_coeff_y,res_l_m3$sampsize,(res_l_m3$fmr_y_int.fmr_y_int)/10,blank)
editfmr_res_l_3=data.frame(res_l_3$mean_betas.data.x_coeff_m,res_l_3$mean_betas.data.xm_coeff_y,res_l_3$sampsize,(res_l_3$fmr_y_int.fmr_y_int)/10,blank)
editfmr_res_l_5=data.frame(res_l_5$mean_betas.data.x_coeff_m,res_l_5$mean_betas.data.xm_coeff_y,res_l_5$sampsize,(res_l_5$fmr_y_int.fmr_y_int)/10,blank)
editfmr_res_l_1=data.frame(res_l_1$mean_betas.data.x_coeff_m,res_l_1$mean_betas.data.xm_coeff_y,res_l_1$sampsize,(res_l_1$fmr_y_int.fmr_y_int)/10,blank)


#interaction coefficient=0
write.table(headers2, file='fmr_L0.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='fmr_L0.txt',append=TRUE, sep = "\n")
write.table(editfmr_res_l_0, file='fmr_L0.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=m3
write.table(headers2, file='fmr_LM3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='fmr_LM3.txt',append=TRUE, sep = "\n")
write.table(editfmr_res_l_m3, file='fmr_LM3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=3
write.table(headers2, file='fmr_L3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='fmr_L3.txt',append=TRUE, sep = "\n")
write.table(editfmr_res_l_3, file='fmr_L3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=5
write.table(headers2, file='fmr_L5.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='fmr_L5.txt',append=TRUE, sep = "\n")
write.table(editfmr_res_l_5, file='fmr_L5.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=1
write.table(headers2, file='fmr_L1.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='fmr_L1.txt',append=TRUE, sep = "\n")
write.table(editfmr_res_l_1, file='fmr_L1.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)



###############
#OBSERVATIONAL#
###############

edit_obsres_l_0=data.frame(res_l_0$mean_betas.data.x_coeff_m,res_l_0$mean_betas.data.xm_coeff_y,res_l_0$sampsize,res_l_0$mean_betas.data.xm_obs,
                           sqrt(res_l_0$s2.newdata.xm_obs),res_l_0$mean_betas.data.xm_obs_se,res_l_0$se.newdata.xm_obs,blank,
                           (res_l_0$xm_obs_detec.xm_obs_detec)/10,(res_l_0$obscoverage.obscoverage)/10)

edit_obsres_l_m3=data.frame(res_l_m3$mean_betas.data.x_coeff_m,res_l_m3$mean_betas.data.xm_coeff_y,res_l_m3$sampsize,res_l_m3$mean_betas.data.xm_obs,
                            sqrt(res_l_m3$s2.newdata.xm_obs),res_l_m3$mean_betas.data.xm_obs_se,res_l_m3$se.newdata.xm_obs,(res_l_m3$xm_obs_detec.xm_obs_detec)/10,
                            blank, (res_l_m3$obscoverage.obscoverage)/10)

edit_obsres_l_3=data.frame(res_l_3$mean_betas.data.x_coeff_m,res_l_3$mean_betas.data.xm_coeff_y,res_l_3$sampsize,res_l_3$mean_betas.data.xm_obs,
                           sqrt(res_l_3$s2.newdata.xm_obs),res_l_3$mean_betas.data.xm_obs_se,res_l_3$se.newdata.xm_obs,(res_l_3$xm_obs_detec.xm_obs_detec)/10,
                           blank, (res_l_3$obscoverage.obscoverage)/10)

edit_obsres_l_5=data.frame(res_l_5$mean_betas.data.x_coeff_m,res_l_5$mean_betas.data.xm_coeff_y,res_l_5$sampsize,res_l_5$mean_betas.data.xm_obs,
                           sqrt(res_l_5$s2.newdata.xm_obs),res_l_5$mean_betas.data.xm_obs_se,res_l_5$se.newdata.xm_obs,(res_l_5$xm_obs_detec.xm_obs_detec)/10,
                           blank, (res_l_5$obscoverage.obscoverage)/10)

edit_obsres_l_1=data.frame(res_l_1$mean_betas.data.x_coeff_m,res_l_1$mean_betas.data.xm_coeff_y,res_l_1$sampsize,res_l_1$mean_betas.data.xm_obs,
                           sqrt(res_l_1$s2.newdata.xm_obs),res_l_1$mean_betas.data.xm_obs_se,res_l_1$se.newdata.xm_obs,(res_l_1$xm_obs_detec.xm_obs_detec)/10,
                           blank, (res_l_1$obscoverage.obscoverage)/10)



#interaction coefficient=0
write.table(headers, file='OBS_L0.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='OBS_L0.txt',append=TRUE, sep = "\n")
write.table(edit_obsres_l_0, file='OBS_L0.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=m3
write.table(headers, file='OBS_LM3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='OBS_LM3.txt',append=TRUE, sep = "\n")
write.table(edit_obsres_l_m3, file='OBS_LM3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=3
write.table(headers, file='OBS_L3.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='OBS_L3.txt',append=TRUE, sep = "\n")
write.table(edit_obsres_l_3, file='OBS_L3.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=5
write.table(headers, file='OBS_L5.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='OBS_L5.txt',append=TRUE, sep = "\n")
write.table(edit_obsres_l_5, file='OBS_L5.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
#interaction coefficient=1
write.table(headers, file='OBS_L1.txt',append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "",eol="")
write("", file='OBS_L1.txt',append=TRUE, sep = "\n")
write.table(edit_obsres_l_1, file='OBS_L1.txt',append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
