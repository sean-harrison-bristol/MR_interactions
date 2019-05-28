#script_name: meta_sensitivity_v8.r
#project: 4-way decomp
#script author: Teri North
#script purpose: simulation to examine impact of pleiotropy on performance of 2sls & FMR
#date created: 10/04/2018
#last edited: 25/01/2019

#read in seed
seedval=commandArgs(6)

print(seedval)
seedval=as.numeric(seedval)

#load packages
library("AER")


sessionInfo()

#define headers
headers=c("x_coeff_m","x_coeff_y","m_coeff_y","xm_coeff_y",
          "fmr_interac_p","xm_2sls","xm_2sls_se")



set.seed(seedval)

#2 models are run with 1000 repeats for each.
#model 1 - as normal
#model 2 - allow for A2 instrument to also causally affect A1

for (rep in c(1:100)){
  
  #log session for each sim repeat
  sink(file=paste(seedval,'_rep',rep,"_meta_log.txt",sep=""))
  
  nval=500000
  
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

  
  for (model in c(1:2)){
    #write file headers 
    #begin with append=FALSE to overwrite current file
    write.table(headers, file=paste(seedval,'_rep',rep,'_model',model,"_pleio_res.txt",sep=""),append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
    
    
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
        if (model==1){
          a1=c+v+0.14*z1
        } else if (model==2){
          a1=c+v+0.1*z1+0.1*z2
        }
        #Mediator
        a2=i*a1+e+0.18*z2+c
        #Outcome
        y=j*a1+k*a2+u+c+l*a1*a2+m*c*a1+n*c*a2
        
        a1a2=a1*a2
        z1z2=z1*z2
        z1z1=z1*z1
        
        #make a dataframe of the required vars
        targ=data.frame(a1,a2,y,a1a2,z1,z2,z1z2,z1z1)


        #split the observations into 4 groups based on median splits of Instruments
        med_z1=median(targ$z1)
        med_z2=median(targ$z2)
          
          
        f0=c(1:length(targ[,1]))
        f1=c(1:length(targ[,1]))
        f2=c(1:length(targ[,1]))
        f3=c(1:length(targ[,1]))
          
        for (z in c(1:length(targ[,1]))){
            
          f0[z]=NA
          f1[z]=NA
          f2[z]=NA
          f3[z]=NA
        }
          
        for (z in c(1:length(targ[,1]))){
          if ((targ$z1[z]<=med_z1) & (targ$z2[z]<=med_z2)){
              
            f0[z]=1
            f1[z]=0
            f2[z]=0
            f3[z]=0
          } else if ((targ$z1[z]<=med_z1) & (targ$z2[z]>med_z2)){
              
            f0[z]=0
            f1[z]=1
            f2[z]=0
            f3[z]=0
          } else if ((targ$z1[z]>med_z1) & (targ$z2[z]<=med_z2)){
              
            f0[z]=0
            f1[z]=0
            f2[z]=1
            f3[z]=0
          } else if ((targ$z1[z]>med_z1) & (targ$z2[z]>med_z2)){
              
            f0[z]=0
            f1[z]=0
            f2[z]=0
            f3[z]=1
          }
        }
        
          
          
        #FACTORIAL MR (median split)
        fitfmr <- lm(targ$y ~ f1 + f2 + f3)
        sum_fmr=coef(summary(fitfmr))
          
          
        #f0 as the baseline - see fig 2 ference paper
          
        #reality check on regression coefficients
        mean(targ$y[f1==1]) - fitfmr$coefficients[2]
        mean(targ$y[f2==1]) - fitfmr$coefficients[3]
        mean(targ$y[f3==1]) - fitfmr$coefficients[4]
          
        mean(targ$y[f0==1]) - fitfmr$coefficients[1] #intercept is mean outcome for baseline group
          
        #formal test for interaction using linear contrast of coeffs - use Wald test comparable to stata's 'test'
        lhyp='1*f1 + 1*f2 - 1*f3 = 0'
        lincon=linearHypothesis(fitfmr,lhyp,test='F')
        int_p=lincon$`Pr(>F)`[2]
        
        fmr_vec=c(int_p)
          

        #instrumental variable 2sls: inst=z1+z2+z1z2+z1z1 [when mediation expected]   
        eq2y<- targ$y ~ targ$a1 + targ$a2 + targ$a1a2 | targ$z1 + targ$z2 + targ$z1z2 + targ$z1z1 
        fitz12sls <- ivreg(eq2y)
        sum_z12sls=coef(summary(fitz12sls))
        
          
          
        tslsz1_vec=c(sum_z12sls[2,1],sum_z12sls[2,2],
                     sum_z12sls[3,1],sum_z12sls[3,2],
                     sum_z12sls[4,1],sum_z12sls[4,2])  
          

        mres=c(coeff_list,fmr_vec,tslsz1_vec[5],tslsz1_vec[6])
        
        
        write("", file=paste(seedval,'_rep',rep,'_model',model,"_pleio_res.txt",sep=""),append=TRUE, sep = "\n")
        write.table(mres, file=paste(seedval,'_rep',rep,'_model',model,"_pleio_res.txt",sep=""),append=TRUE, row.names=FALSE, col.names=FALSE, sep = "\t",eol="\t")
        
        
        
      }
    }
    
    
    
    
    
    
  }      
  
  sink()
}





