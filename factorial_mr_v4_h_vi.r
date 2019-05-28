#script_name: factorial_mr_v4_h_vi.r
#project: 4-way decomp
#script author: Teri North
#script purpose: simulation to examine performance of Ference (2015) factorial MR approach - same as v2 script but with coeffs to impose 1% variance for z1 and z2
#date created: 12/04/2017
#last edited: 11/1/2019
#notes:
#https://cran.r-project.org/web/packages/systemfit/systemfit.pdf
#https://stat.ethz.ch/R-manual/R-devel/library/utils/html/write.table.html
#https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html

#read in seed
seedval=commandArgs(6)
print(seedval)
seedval=as.numeric(seedval)

library("AER")
library("ivpack")

sessionInfo()

#define headers
headers=c("x_coeff_m","x_coeff_y","m_coeff_y","xm_coeff_y",
          "x_obs", "x_obs_se",
          "m_obs", "m_obs_se",
          "xm_obs", "xm_obs_se",
          "f1", "f1_se",
          "f2", "f2_se",
          "f3","f3_se",
          "fmr_interac_p",
          "x_2sls", "x_2sls_se", 
          "m_2sls", "m_2sls_se",  
          "xm_2sls", "xm_2sls_se", 
          "x_z1_2sls", "x_z1_2sls_se", 
          "m_z1_2sls", "m_z1_2sls_se",  
          "xm_z1_2sls", "xm_z1_2sls_se")


rob_headers=c("x_coeff_m","x_coeff_y","m_coeff_y","xm_coeff_y",
              "x_obs", "x_obs_se",
              "m_obs", "m_obs_se",
              "xm_obs", "xm_obs_se",
              "f1", "f1_se",
              "f2", "f2_se",
              "f3","f3_se",
              "fmr_interac_p",
              "x_2sls", "x_2sls_se", 
              "m_2sls", "m_2sls_se",  
              "xm_2sls", "xm_2sls_se", 
              "x_z1_2sls", "x_z1_2sls_se", 
              "m_z1_2sls", "m_z1_2sls_se",  
              "xm_z1_2sls", "xm_z1_2sls_se")

F_headers=c("x_coeff_m","x_coeff_y","m_coeff_y","xm_coeff_y",
            "F_2SLS", "F_Z1_2SLS")



set.seed(seedval)

for (nval in c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,500000,1000000)){
  for (rep in c(1:100)){
    #log session for each seed, sample size and sim repeat
    sink(file=paste(seedval,'_rep',rep,'_samp',nval,"_FMRlog.txt",sep=""))
    
    #write file headers 
    #begin with append=FALSE to overwrite current file
    write.table(headers, file=paste(seedval,'_rep',rep,'_samp',nval,"_FMRres.txt",sep=""),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
    write.table(rob_headers, file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_ROBUST_res.txt",sep=""),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
    write.table(F_headers, file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_F_res.txt",sep=""),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
    
    #generate data for system
    samp_size=nval
    #error terms
    u=rnorm(samp_size)
    v=rnorm(samp_size)
    e=rnorm(samp_size)
    #confounders
    c=rnorm(samp_size)
    #instruments
    z1=rnorm(samp_size)
    z2=rnorm(samp_size)
    
    
    #begin looping over coefficient combinations
    
    for (i in c(0,0.333,0.5,1,-0.333)) {
      
      for (othcomb in c(1,2,3,4,5)){
        
        if (othcomb==1){
          j=0
          k=0
          l=0/3
          m=0/3
          n=0/3
          
        } else if (othcomb==2){
          j=0.333
          k=0.333
          l=0.333/3
          m=0.333/3
          n=0.333/3
          
        } else if (othcomb==3){
          j=0.5
          k=0.5
          l=0.5/3
          m=0.5/3
          n=0.5/3
          
        } else if (othcomb==4){
          j=1
          k=1
          l=1/3
          m=1/3
          n=1/3
          
        } else if (othcomb==5){
          j=-0.333
          k=-0.333
          l=-0.333/3
          m=-0.333/3
          n=-0.333/3
          
        }
        
        
        
        print(i)
        print(j)
        print(k)
        print(l)
        coeff_list<- c(i,j,k,l)
        
        #define models
        #Exposure
        a1=c+v+0.14*z1
        #Mediator
        a2=i*a1+e+0.18*z2+c
        #Outcome
        y=j*a1+k*a2+u+c+l*a1*a2+m*c*a1+n*c*a2
        
        a1a2=a1*a2
        z1z2=z1*z2
        z1z1=z1*z1
        z2z2=z2*z2  
        
        #observational
        fitols <- lm(y ~ a1 + a2 + a1a2)
        sum_ols=coef(summary(fitols))
        rob_fitols <- coeftest(fitols,vcov=vcovHC(fitols,type="HC0"))
        
        obs_vec=c(sum_ols[2,1],sum_ols[2,2],
                  sum_ols[3,1],sum_ols[3,2],
                  sum_ols[4,1],sum_ols[4,2])
        
        rob_obs_vec=c(rob_fitols[2,1],rob_fitols[2,2],
                      rob_fitols[3,1],rob_fitols[3,2],
                      rob_fitols[4,1],rob_fitols[4,2])
        
        
        
        #split the observations into 4 groups based on median splits of Instruments
        med_z1=median(z1)
        med_z2=median(z2)
        
        fact_group=c(1:nval)
        f0=c(1:nval)
        f1=c(1:nval)
        f2=c(1:nval)
        f3=c(1:nval)
        
        for (z in c(1:nval)){
          fact_group[z]=NA
          f0[z]=NA
          f1[z]=NA
          f2[z]=NA
          f3[z]=NA
        }
        
        for (z in c(1:nval)){
          if ((z1[z]<=med_z1) & (z2[z]<=med_z2)){
            fact_group[z]=0
            f0[z]=1
            f1[z]=0
            f2[z]=0
            f3[z]=0
          } else if ((z1[z]<=med_z1) & (z2[z]>med_z2)){
            fact_group[z]=1
            f0[z]=0
            f1[z]=1
            f2[z]=0
            f3[z]=0
          } else if ((z1[z]>med_z1) & (z2[z]<=med_z2)){
            fact_group[z]=2
            f0[z]=0
            f1[z]=0
            f2[z]=1
            f3[z]=0
          } else if ((z1[z]>med_z1) & (z2[z]>med_z2)){
            fact_group[z]=3
            f0[z]=0
            f1[z]=0
            f2[z]=0
            f3[z]=1
          }
        }
        
        
        
        #FACTORIAL MR (median split)
        fitfmr <- lm(y ~ f1 + f2 + f3)
        sum_fmr=coef(summary(fitfmr))
        rob_fitfmr <- coeftest(fitfmr,vcov=vcovHC(fitfmr,type="HC0"))
        
        #f0 as the baseline - see fig 2 reference paper
        
        #reality check on regression coefficients
        mean(y[f1==1]) - fitfmr$coefficients[2]
        mean(y[f2==1]) - fitfmr$coefficients[3]
        mean(y[f3==1]) - fitfmr$coefficients[4]
        
        mean(y[f0==1]) - fitfmr$coefficients[1] #intercept is mean outcome for baseline group
        
        #formal test for interaction using linear contrast of coeffs - use Wald test comparable to stata's 'test'
        lhyp='1*f1 + 1*f2 - 1*f3 = 0'
        lincon=linearHypothesis(fitfmr,lhyp,test='F')
        int_p=lincon$`Pr(>F)`[2]
        
        fmr_vec=c(sum_fmr[2,1],sum_fmr[2,2], #f1 vs f0
                  sum_fmr[3,1],sum_fmr[3,2], #f2 vs f0
                  sum_fmr[4,1],sum_fmr[4,2], #f3 vs f0
                  int_p)
        
        rob_fmr_vec=c(rob_fitfmr[2,1],rob_fitfmr[2,2], #f1 vs f0
                      rob_fitfmr[3,1],rob_fitfmr[3,2], #f2 vs f0
                      rob_fitfmr[4,1],rob_fitfmr[4,2], #f3 vs f0
                      "NA")
        
        
        #instrumental variable 2sls: inst=z1+z2+z1z2 [when no mediation expected]
        eqy<- y ~ a1 + a2 + a1a2 | z1 + z2 + z1z2 
        fit2sls <- ivreg(eqy)
        sum_2sls=coef(summary(fit2sls))
        
        #hetroskedastic robust SEs
        rob_fit2sls<-robust.se(fit2sls)
        
        tsls_vec=c(sum_2sls[2,1],sum_2sls[2,2],
                   sum_2sls[3,1],sum_2sls[3,2],
                   sum_2sls[4,1],sum_2sls[4,2])
        
        rob_tsls_vec=c(rob_fit2sls[2,1],rob_fit2sls[2,2],
                       rob_fit2sls[3,1],rob_fit2sls[3,2],
                       rob_fit2sls[4,1],rob_fit2sls[4,2])
        
        freg = ivreg(a1a2 ~ a1 + a2 | z1 + z2 + z1z2)
        fregres = freg$residuals
        F_a1a2 = summary(lm(fregres ~ z1 + z2 + z1z2))$fstatistic[1]
        F_adj_a1a2=F_a1a2*(3/(3-3+1))
        
        
        #instrumental variable 2sls: inst=z1+z2+z1z2+z1z1 [when mediation expected]          
        eq2y<- y ~ a1 + a2 + a1a2 | z1 + z2 + z1z2 + z1z1 
        fitz12sls <- ivreg(eq2y)
        sum_z12sls=coef(summary(fitz12sls))
        
        #hetroskedastic robust SEs
        rob_fitz12sls<-robust.se(fitz12sls)
        
        tslsz1_vec=c(sum_z12sls[2,1],sum_z12sls[2,2],
                     sum_z12sls[3,1],sum_z12sls[3,2],
                     sum_z12sls[4,1],sum_z12sls[4,2])
        
        rob_tslsz1_vec=c(rob_fitz12sls[2,1],rob_fitz12sls[2,2],
                         rob_fitz12sls[3,1],rob_fitz12sls[3,2],
                         rob_fitz12sls[4,1],rob_fitz12sls[4,2])
        
        
        z1freg = ivreg(a1a2 ~ a1 + a2 | z1 + z2 + z1z2 + z1z1)
        z1fregres = z1freg$residuals
        Fz1_a1a2 = summary(lm(z1fregres ~ z1 + z2 + z1z2 + z1z1))$fstatistic[1]
        Fz1_adj_a1a2=Fz1_a1a2*(4/(4-3+1))
        
        
        
        #all results
        res_vec=c(coeff_list,obs_vec,fmr_vec,tsls_vec,tslsz1_vec)
        rob_res_vec=c(coeff_list,rob_obs_vec,rob_fmr_vec,rob_tsls_vec,rob_tslsz1_vec)
        f_res_vec=c(coeff_list,F_adj_a1a2,Fz1_adj_a1a2)
        
        write("", file=paste(seedval,'_rep',rep,'_samp',nval,"_FMRres.txt",sep=""),append=TRUE, sep = "\n")
        write.table(res_vec, file=paste(seedval,'_rep',rep,'_samp',nval,"_FMRres.txt",sep=""),append=TRUE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
        
        write("", file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_ROBUST_res.txt",sep=""),append=TRUE, sep = "\n")
        write.table(rob_res_vec, file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_ROBUST_res.txt",sep=""),append=TRUE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
        
        write("", file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_F_res.txt",sep=""),append=TRUE, sep = "\n")
        write.table(f_res_vec, file=paste(seedval,'_rep',rep,'_samp',nval,"_FMR_F_res.txt",sep=""),append=TRUE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
        
        
        
        
      }
      
    }
    
    sink()
  }
}



