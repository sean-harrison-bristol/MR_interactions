#script_name: forest_plot_PAP1_V4.r
#project: 4-way decomp: paper 1
#script author: Teri North
#script purpose: plot results for mr interactions paper
#date created: 22/01/2019
#last edited:02/02/2019
#notes:

writeplot='' #Folder 2
dataloc='' #Folder 1
setwd(dataloc)

#COEFFICIENT PLOTS
#PREP THE N=50K DATA
N_50=data.frame(read.table(file=paste(50000,"_EXTRA_final_res.txt",sep=''),sep='\t',header=TRUE))
N_50=N_50[order(N_50$mean_betas.data.xm_coeff_y,N_50$mean_betas.data.x_coeff_m),]
tsls_xm_lower_N_50=N_50$mean_betas.data.xm_2sls-1.96*N_50$se.newdata.xm_2sls
tsls_xm_upper_N_50=N_50$mean_betas.data.xm_2sls+1.96*N_50$se.newdata.xm_2sls
tsls_z1_xm_lower_N_50=N_50$mean_betas.data.xm_z1_2sls-1.96*N_50$se.newdata.xm_z1_2sls
tsls_z1_xm_upper_N_50=N_50$mean_betas.data.xm_z1_2sls+1.96*N_50$se.newdata.xm_z1_2sls
obs_xm_lower_N_50=N_50$mean_betas.data.xm_obs-1.96*N_50$se.newdata.xm_obs
obs_xm_upper_N_50=N_50$mean_betas.data.xm_obs+1.96*N_50$se.newdata.xm_obs
N_50=data.frame(N_50,tsls_xm_lower_N_50,tsls_xm_upper_N_50,tsls_z1_xm_lower_N_50,tsls_z1_xm_upper_N_50,
                obs_xm_lower_N_50,obs_xm_upper_N_50)
N_50_interac_m3=N_50[round(N_50$mean_betas.data.xm_coeff_y,3)==-0.111,]
N_50_interac_0=N_50[round(N_50$mean_betas.data.xm_coeff_y,3)==0,]
N_50_interac_3=N_50[round(N_50$mean_betas.data.xm_coeff_y,3)==0.111,]
N_50_interac_5=N_50[round(N_50$mean_betas.data.xm_coeff_y,3)==0.167,]
N_50_interac_1=N_50[round(N_50$mean_betas.data.xm_coeff_y,3)==0.333,]

#2SLS Z=(Z1,Z2,Z1Z2) N=50K
setwd(writeplot)
tiff('coef_plot_z1z2_50k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_50_interac_m3$mean_betas.data.xm_2sls,N_50_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-2,2),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='2sls estimate of theta coefficient, Z=(Z1,Z2,Z1Z2)',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_50_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-2,-1.5,-1,-0.5,-0.111,0,0.5,1,1.5,2),labels=c('-2','-1.5','-1','-0.5','-0.111','0','0.5','1','1.5','2'),cex.axis=0.25)
lines(c(N_50_interac_m3$tsls_xm_lower_N_50[1],N_50_interac_m3$tsls_xm_upper_N_50[1]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[1],N_50_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_m3$tsls_xm_lower_N_50[2],N_50_interac_m3$tsls_xm_upper_N_50[2]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[2],N_50_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_m3$tsls_xm_lower_N_50[3],N_50_interac_m3$tsls_xm_upper_N_50[3]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[3],N_50_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_m3$tsls_xm_lower_N_50[4],N_50_interac_m3$tsls_xm_upper_N_50[4]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[4],N_50_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_m3$tsls_xm_lower_N_50[5],N_50_interac_m3$tsls_xm_upper_N_50[5]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[5],N_50_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_0$mean_betas.data.xm_2sls,N_50_interac_0$mean_betas.data.x_coeff_m,xlim=c(-2,2),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_50_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2),labels=c('-2','-1.5','-1','-0.5','0','0.5','1','1.5','2'),cex.axis=0.25)
lines(c(N_50_interac_0$tsls_xm_lower_N_50[1],N_50_interac_0$tsls_xm_upper_N_50[1]),c(N_50_interac_0$mean_betas.data.x_coeff_m[1],N_50_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_0$tsls_xm_lower_N_50[2],N_50_interac_0$tsls_xm_upper_N_50[2]),c(N_50_interac_0$mean_betas.data.x_coeff_m[2],N_50_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_0$tsls_xm_lower_N_50[3],N_50_interac_0$tsls_xm_upper_N_50[3]),c(N_50_interac_0$mean_betas.data.x_coeff_m[3],N_50_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_0$tsls_xm_lower_N_50[4],N_50_interac_0$tsls_xm_upper_N_50[4]),c(N_50_interac_0$mean_betas.data.x_coeff_m[4],N_50_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_0$tsls_xm_lower_N_50[5],N_50_interac_0$tsls_xm_upper_N_50[5]),c(N_50_interac_0$mean_betas.data.x_coeff_m[5],N_50_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_3$mean_betas.data.xm_2sls,N_50_interac_3$mean_betas.data.x_coeff_m,xlim=c(-2,2),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_50_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-2,-1.5,-1,-0.5,0,0.111,0.5,1,1.5,2),labels=c('-2','-1.5','-1','-0.5','0','0.111','0.5','1','1.5','2'),cex.axis=0.25)
lines(c(N_50_interac_3$tsls_xm_lower_N_50[1],N_50_interac_3$tsls_xm_upper_N_50[1]),c(N_50_interac_3$mean_betas.data.x_coeff_m[1],N_50_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_3$tsls_xm_lower_N_50[2],N_50_interac_3$tsls_xm_upper_N_50[2]),c(N_50_interac_3$mean_betas.data.x_coeff_m[2],N_50_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_3$tsls_xm_lower_N_50[3],N_50_interac_3$tsls_xm_upper_N_50[3]),c(N_50_interac_3$mean_betas.data.x_coeff_m[3],N_50_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_3$tsls_xm_lower_N_50[4],N_50_interac_3$tsls_xm_upper_N_50[4]),c(N_50_interac_3$mean_betas.data.x_coeff_m[4],N_50_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_3$tsls_xm_lower_N_50[5],N_50_interac_3$tsls_xm_upper_N_50[5]),c(N_50_interac_3$mean_betas.data.x_coeff_m[5],N_50_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_5$mean_betas.data.xm_2sls,N_50_interac_5$mean_betas.data.x_coeff_m,xlim=c(-2,2),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_50_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-2,-1.5,-1,-0.5,0,0.167,0.5,1,1.5,2),labels=c('-2','-1.5','-1','-0.5','0','0.167','0.5','1','1.5','2'),cex.axis=0.25)
lines(c(N_50_interac_5$tsls_xm_lower_N_50[1],N_50_interac_5$tsls_xm_upper_N_50[1]),c(N_50_interac_5$mean_betas.data.x_coeff_m[1],N_50_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_5$tsls_xm_lower_N_50[2],N_50_interac_5$tsls_xm_upper_N_50[2]),c(N_50_interac_5$mean_betas.data.x_coeff_m[2],N_50_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_5$tsls_xm_lower_N_50[3],N_50_interac_5$tsls_xm_upper_N_50[3]),c(N_50_interac_5$mean_betas.data.x_coeff_m[3],N_50_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_5$tsls_xm_lower_N_50[4],N_50_interac_5$tsls_xm_upper_N_50[4]),c(N_50_interac_5$mean_betas.data.x_coeff_m[4],N_50_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_5$tsls_xm_lower_N_50[5],N_50_interac_5$tsls_xm_upper_N_50[5]),c(N_50_interac_5$mean_betas.data.x_coeff_m[5],N_50_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_1$mean_betas.data.xm_2sls,N_50_interac_1$mean_betas.data.x_coeff_m,xlim=c(-2,2),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_50_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-2,-1.5,-1,-0.5,0,0.333,0.5,1,1.5,2),labels=c('-2','-1.5','-1','-0.5','0','0.333','0.5','1','1.5','2'),cex.axis=0.25)
lines(c(N_50_interac_1$tsls_xm_lower_N_50[1],N_50_interac_1$tsls_xm_upper_N_50[1]),c(N_50_interac_1$mean_betas.data.x_coeff_m[1],N_50_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_1$tsls_xm_lower_N_50[2],N_50_interac_1$tsls_xm_upper_N_50[2]),c(N_50_interac_1$mean_betas.data.x_coeff_m[2],N_50_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_1$tsls_xm_lower_N_50[3],N_50_interac_1$tsls_xm_upper_N_50[3]),c(N_50_interac_1$mean_betas.data.x_coeff_m[3],N_50_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_1$tsls_xm_lower_N_50[4],N_50_interac_1$tsls_xm_upper_N_50[4]),c(N_50_interac_1$mean_betas.data.x_coeff_m[4],N_50_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_1$tsls_xm_lower_N_50[5],N_50_interac_1$tsls_xm_upper_N_50[5]),c(N_50_interac_1$mean_betas.data.x_coeff_m[5],N_50_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()

#2SLS Z=(Z1,Z2,Z1Z2,Z1Z1) N=50K
tiff('coef_plot_z1z1_50k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_50_interac_m3$mean_betas.data.xm_z1_2sls,N_50_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.25,0.45),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='2sls estimate of theta coefficient, Z=(Z1,Z2,Z1Z2,Z1Z1)',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_50_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.25,-0.2,-0.15,-0.111,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),labels=c('-0.25','-0.2','-0.15','-0.111','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45'),cex.axis=0.25)
lines(c(N_50_interac_m3$tsls_z1_xm_lower_N_50[1],N_50_interac_m3$tsls_z1_xm_upper_N_50[1]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[1],N_50_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_m3$tsls_z1_xm_lower_N_50[2],N_50_interac_m3$tsls_z1_xm_upper_N_50[2]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[2],N_50_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_m3$tsls_z1_xm_lower_N_50[3],N_50_interac_m3$tsls_z1_xm_upper_N_50[3]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[3],N_50_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_m3$tsls_z1_xm_lower_N_50[4],N_50_interac_m3$tsls_z1_xm_upper_N_50[4]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[4],N_50_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_m3$tsls_z1_xm_lower_N_50[5],N_50_interac_m3$tsls_z1_xm_upper_N_50[5]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[5],N_50_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_0$mean_betas.data.xm_z1_2sls,N_50_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.25,0.45),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_50_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45),labels=c('-0.25','-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3','0.35','0.4','0.45'),cex.axis=0.25)
lines(c(N_50_interac_0$tsls_z1_xm_lower_N_50[1],N_50_interac_0$tsls_z1_xm_upper_N_50[1]),c(N_50_interac_0$mean_betas.data.x_coeff_m[1],N_50_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_0$tsls_z1_xm_lower_N_50[2],N_50_interac_0$tsls_z1_xm_upper_N_50[2]),c(N_50_interac_0$mean_betas.data.x_coeff_m[2],N_50_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_0$tsls_z1_xm_lower_N_50[3],N_50_interac_0$tsls_z1_xm_upper_N_50[3]),c(N_50_interac_0$mean_betas.data.x_coeff_m[3],N_50_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_0$tsls_z1_xm_lower_N_50[4],N_50_interac_0$tsls_z1_xm_upper_N_50[4]),c(N_50_interac_0$mean_betas.data.x_coeff_m[4],N_50_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_0$tsls_z1_xm_lower_N_50[5],N_50_interac_0$tsls_z1_xm_upper_N_50[5]),c(N_50_interac_0$mean_betas.data.x_coeff_m[5],N_50_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_3$mean_betas.data.xm_z1_2sls,N_50_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.25,0.45),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_50_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.111,0.15,0.2,0.25,0.3,0.35,0.4,0.45),labels=c('-0.25','-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.111','0.15','0.2','0.25','0.3','0.35','0.4','0.45'),cex.axis=0.25)
lines(c(N_50_interac_3$tsls_z1_xm_lower_N_50[1],N_50_interac_3$tsls_z1_xm_upper_N_50[1]),c(N_50_interac_3$mean_betas.data.x_coeff_m[1],N_50_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_3$tsls_z1_xm_lower_N_50[2],N_50_interac_3$tsls_z1_xm_upper_N_50[2]),c(N_50_interac_3$mean_betas.data.x_coeff_m[2],N_50_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_3$tsls_z1_xm_lower_N_50[3],N_50_interac_3$tsls_z1_xm_upper_N_50[3]),c(N_50_interac_3$mean_betas.data.x_coeff_m[3],N_50_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_3$tsls_z1_xm_lower_N_50[4],N_50_interac_3$tsls_z1_xm_upper_N_50[4]),c(N_50_interac_3$mean_betas.data.x_coeff_m[4],N_50_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_3$tsls_z1_xm_lower_N_50[5],N_50_interac_3$tsls_z1_xm_upper_N_50[5]),c(N_50_interac_3$mean_betas.data.x_coeff_m[5],N_50_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_5$mean_betas.data.xm_z1_2sls,N_50_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.25,0.45),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_50_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.167,0.2,0.25,0.3,0.35,0.4,0.45),labels=c('-0.25','-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.167','0.2','0.25','0.3','0.35','0.4','0.45'),cex.axis=0.25)
lines(c(N_50_interac_5$tsls_z1_xm_lower_N_50[1],N_50_interac_5$tsls_z1_xm_upper_N_50[1]),c(N_50_interac_5$mean_betas.data.x_coeff_m[1],N_50_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_5$tsls_z1_xm_lower_N_50[2],N_50_interac_5$tsls_z1_xm_upper_N_50[2]),c(N_50_interac_5$mean_betas.data.x_coeff_m[2],N_50_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_5$tsls_z1_xm_lower_N_50[3],N_50_interac_5$tsls_z1_xm_upper_N_50[3]),c(N_50_interac_5$mean_betas.data.x_coeff_m[3],N_50_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_5$tsls_z1_xm_lower_N_50[4],N_50_interac_5$tsls_z1_xm_upper_N_50[4]),c(N_50_interac_5$mean_betas.data.x_coeff_m[4],N_50_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_5$tsls_z1_xm_lower_N_50[5],N_50_interac_5$tsls_z1_xm_upper_N_50[5]),c(N_50_interac_5$mean_betas.data.x_coeff_m[5],N_50_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_1$mean_betas.data.xm_z1_2sls,N_50_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.25,0.45),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_50_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.333,0.35,0.4,0.45),labels=c('-0.25','-0.2','-0.15','-0.1','-0.05','0','0.05','0.1','0.15','0.2','0.25','0.3','0.333','0.35','0.4','0.45'),cex.axis=0.25)
lines(c(N_50_interac_1$tsls_z1_xm_lower_N_50[1],N_50_interac_1$tsls_z1_xm_upper_N_50[1]),c(N_50_interac_1$mean_betas.data.x_coeff_m[1],N_50_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_1$tsls_z1_xm_lower_N_50[2],N_50_interac_1$tsls_z1_xm_upper_N_50[2]),c(N_50_interac_1$mean_betas.data.x_coeff_m[2],N_50_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_1$tsls_z1_xm_lower_N_50[3],N_50_interac_1$tsls_z1_xm_upper_N_50[3]),c(N_50_interac_1$mean_betas.data.x_coeff_m[3],N_50_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_1$tsls_z1_xm_lower_N_50[4],N_50_interac_1$tsls_z1_xm_upper_N_50[4]),c(N_50_interac_1$mean_betas.data.x_coeff_m[4],N_50_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_1$tsls_z1_xm_lower_N_50[5],N_50_interac_1$tsls_z1_xm_upper_N_50[5]),c(N_50_interac_1$mean_betas.data.x_coeff_m[5],N_50_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()

#OBSERVATIONAL N=50K
tiff('coef_plot_obs_50k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_50_interac_m3$mean_betas.data.xm_obs,N_50_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='Ordinary least sqaures estimate of theta coefficient',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_50_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.111,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.111','-0.1','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_50_interac_m3$obs_xm_lower_N_50[1],N_50_interac_m3$obs_xm_upper_N_50[1]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[1],N_50_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_m3$obs_xm_lower_N_50[2],N_50_interac_m3$obs_xm_upper_N_50[2]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[2],N_50_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_m3$obs_xm_lower_N_50[3],N_50_interac_m3$obs_xm_upper_N_50[3]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[3],N_50_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_m3$obs_xm_lower_N_50[4],N_50_interac_m3$obs_xm_upper_N_50[4]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[4],N_50_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_m3$obs_xm_lower_N_50[5],N_50_interac_m3$obs_xm_upper_N_50[5]),c(N_50_interac_m3$mean_betas.data.x_coeff_m[5],N_50_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_0$mean_betas.data.xm_obs,N_50_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_50_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_50_interac_0$obs_xm_lower_N_50[1],N_50_interac_0$obs_xm_upper_N_50[1]),c(N_50_interac_0$mean_betas.data.x_coeff_m[1],N_50_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_0$obs_xm_lower_N_50[2],N_50_interac_0$obs_xm_upper_N_50[2]),c(N_50_interac_0$mean_betas.data.x_coeff_m[2],N_50_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_0$obs_xm_lower_N_50[3],N_50_interac_0$obs_xm_upper_N_50[3]),c(N_50_interac_0$mean_betas.data.x_coeff_m[3],N_50_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_0$obs_xm_lower_N_50[4],N_50_interac_0$obs_xm_upper_N_50[4]),c(N_50_interac_0$mean_betas.data.x_coeff_m[4],N_50_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_0$obs_xm_lower_N_50[5],N_50_interac_0$obs_xm_upper_N_50[5]),c(N_50_interac_0$mean_betas.data.x_coeff_m[5],N_50_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_3$mean_betas.data.xm_obs,N_50_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_50_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.111,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.111','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_50_interac_3$obs_xm_lower_N_50[1],N_50_interac_3$obs_xm_upper_N_50[1]),c(N_50_interac_3$mean_betas.data.x_coeff_m[1],N_50_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_3$obs_xm_lower_N_50[2],N_50_interac_3$obs_xm_upper_N_50[2]),c(N_50_interac_3$mean_betas.data.x_coeff_m[2],N_50_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_3$obs_xm_lower_N_50[3],N_50_interac_3$obs_xm_upper_N_50[3]),c(N_50_interac_3$mean_betas.data.x_coeff_m[3],N_50_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_3$obs_xm_lower_N_50[4],N_50_interac_3$obs_xm_upper_N_50[4]),c(N_50_interac_3$mean_betas.data.x_coeff_m[4],N_50_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_3$obs_xm_lower_N_50[5],N_50_interac_3$obs_xm_upper_N_50[5]),c(N_50_interac_3$mean_betas.data.x_coeff_m[5],N_50_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_5$mean_betas.data.xm_obs,N_50_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_50_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.167,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.167','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_50_interac_5$obs_xm_lower_N_50[1],N_50_interac_5$obs_xm_upper_N_50[1]),c(N_50_interac_5$mean_betas.data.x_coeff_m[1],N_50_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_5$obs_xm_lower_N_50[2],N_50_interac_5$obs_xm_upper_N_50[2]),c(N_50_interac_5$mean_betas.data.x_coeff_m[2],N_50_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_5$obs_xm_lower_N_50[3],N_50_interac_5$obs_xm_upper_N_50[3]),c(N_50_interac_5$mean_betas.data.x_coeff_m[3],N_50_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_5$obs_xm_lower_N_50[4],N_50_interac_5$obs_xm_upper_N_50[4]),c(N_50_interac_5$mean_betas.data.x_coeff_m[4],N_50_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_5$obs_xm_lower_N_50[5],N_50_interac_5$obs_xm_upper_N_50[5]),c(N_50_interac_5$mean_betas.data.x_coeff_m[5],N_50_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_50_interac_1$mean_betas.data.xm_obs,N_50_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_50_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.333,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.333','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_50_interac_1$obs_xm_lower_N_50[1],N_50_interac_1$obs_xm_upper_N_50[1]),c(N_50_interac_1$mean_betas.data.x_coeff_m[1],N_50_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_50_interac_1$obs_xm_lower_N_50[2],N_50_interac_1$obs_xm_upper_N_50[2]),c(N_50_interac_1$mean_betas.data.x_coeff_m[2],N_50_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_50_interac_1$obs_xm_lower_N_50[3],N_50_interac_1$obs_xm_upper_N_50[3]),c(N_50_interac_1$mean_betas.data.x_coeff_m[3],N_50_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_50_interac_1$obs_xm_lower_N_50[4],N_50_interac_1$obs_xm_upper_N_50[4]),c(N_50_interac_1$mean_betas.data.x_coeff_m[4],N_50_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_50_interac_1$obs_xm_lower_N_50[5],N_50_interac_1$obs_xm_upper_N_50[5]),c(N_50_interac_1$mean_betas.data.x_coeff_m[5],N_50_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()

#PREP THE N=100K DATA
setwd(dataloc)
N_100=data.frame(read.table(file=paste(100000,"_EXTRA_final_res.txt",sep=''),sep='\t',header=TRUE))
N_100=N_100[order(N_100$mean_betas.data.xm_coeff_y,N_100$mean_betas.data.x_coeff_m),]
tsls_xm_lower_N_100=N_100$mean_betas.data.xm_2sls-1.96*N_100$se.newdata.xm_2sls
tsls_xm_upper_N_100=N_100$mean_betas.data.xm_2sls+1.96*N_100$se.newdata.xm_2sls
tsls_z1_xm_lower_N_100=N_100$mean_betas.data.xm_z1_2sls-1.96*N_100$se.newdata.xm_z1_2sls
tsls_z1_xm_upper_N_100=N_100$mean_betas.data.xm_z1_2sls+1.96*N_100$se.newdata.xm_z1_2sls
obs_xm_lower_N_100=N_100$mean_betas.data.xm_obs-1.96*N_100$se.newdata.xm_obs
obs_xm_upper_N_100=N_100$mean_betas.data.xm_obs+1.96*N_100$se.newdata.xm_obs
N_100=data.frame(N_100,tsls_xm_lower_N_100,tsls_xm_upper_N_100,tsls_z1_xm_lower_N_100,tsls_z1_xm_upper_N_100,
                obs_xm_lower_N_100,obs_xm_upper_N_100)
N_100_interac_m3=N_100[round(N_100$mean_betas.data.xm_coeff_y,3)==-0.111,]
N_100_interac_0=N_100[round(N_100$mean_betas.data.xm_coeff_y,3)==0,]
N_100_interac_3=N_100[round(N_100$mean_betas.data.xm_coeff_y,3)==0.111,]
N_100_interac_5=N_100[round(N_100$mean_betas.data.xm_coeff_y,3)==0.167,]
N_100_interac_1=N_100[round(N_100$mean_betas.data.xm_coeff_y,3)==0.333,]

#2SLS Z=(Z1,Z2,Z1Z2) N=100K
setwd(writeplot)
tiff('coef_plot_z1z2_100k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_100_interac_m3$mean_betas.data.xm_2sls,N_100_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.6,0.5),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='2sls estimate of theta coefficient, Z=(Z1,Z2,Z1Z2)',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_100_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.111,-0.1,0,0.1,0.2,0.3,0.4,0.5),labels=c('-0.6','-0.5','-0.4','-0.3','-0.2','-0.111','-0.1','0','0.1','0.2','0.3','0.4','0.5'),cex.axis=0.25)
lines(c(N_100_interac_m3$tsls_xm_lower_N_100[1],N_100_interac_m3$tsls_xm_upper_N_100[1]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[1],N_100_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_m3$tsls_xm_lower_N_100[2],N_100_interac_m3$tsls_xm_upper_N_100[2]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[2],N_100_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_m3$tsls_xm_lower_N_100[3],N_100_interac_m3$tsls_xm_upper_N_100[3]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[3],N_100_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_m3$tsls_xm_lower_N_100[4],N_100_interac_m3$tsls_xm_upper_N_100[4]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[4],N_100_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_m3$tsls_xm_lower_N_100[5],N_100_interac_m3$tsls_xm_upper_N_100[5]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[5],N_100_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_0$mean_betas.data.xm_2sls,N_100_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.6,0.5),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_100_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5),labels=c('-0.6','-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5'),cex.axis=0.25)
lines(c(N_100_interac_0$tsls_xm_lower_N_100[1],N_100_interac_0$tsls_xm_upper_N_100[1]),c(N_100_interac_0$mean_betas.data.x_coeff_m[1],N_100_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_0$tsls_xm_lower_N_100[2],N_100_interac_0$tsls_xm_upper_N_100[2]),c(N_100_interac_0$mean_betas.data.x_coeff_m[2],N_100_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_0$tsls_xm_lower_N_100[3],N_100_interac_0$tsls_xm_upper_N_100[3]),c(N_100_interac_0$mean_betas.data.x_coeff_m[3],N_100_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_0$tsls_xm_lower_N_100[4],N_100_interac_0$tsls_xm_upper_N_100[4]),c(N_100_interac_0$mean_betas.data.x_coeff_m[4],N_100_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_0$tsls_xm_lower_N_100[5],N_100_interac_0$tsls_xm_upper_N_100[5]),c(N_100_interac_0$mean_betas.data.x_coeff_m[5],N_100_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_3$mean_betas.data.xm_2sls,N_100_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.6,0.5),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_100_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.111,0.2,0.3,0.4,0.5),labels=c('-0.6','-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.111','0.2','0.3','0.4','0.5'),cex.axis=0.25)
lines(c(N_100_interac_3$tsls_xm_lower_N_100[1],N_100_interac_3$tsls_xm_upper_N_100[1]),c(N_100_interac_3$mean_betas.data.x_coeff_m[1],N_100_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_3$tsls_xm_lower_N_100[2],N_100_interac_3$tsls_xm_upper_N_100[2]),c(N_100_interac_3$mean_betas.data.x_coeff_m[2],N_100_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_3$tsls_xm_lower_N_100[3],N_100_interac_3$tsls_xm_upper_N_100[3]),c(N_100_interac_3$mean_betas.data.x_coeff_m[3],N_100_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_3$tsls_xm_lower_N_100[4],N_100_interac_3$tsls_xm_upper_N_100[4]),c(N_100_interac_3$mean_betas.data.x_coeff_m[4],N_100_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_3$tsls_xm_lower_N_100[5],N_100_interac_3$tsls_xm_upper_N_100[5]),c(N_100_interac_3$mean_betas.data.x_coeff_m[5],N_100_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_5$mean_betas.data.xm_2sls,N_100_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.6,0.5),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_100_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.167,0.2,0.3,0.4,0.5),labels=c('-0.6','-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.167','0.2','0.3','0.4','0.5'),cex.axis=0.25)
lines(c(N_100_interac_5$tsls_xm_lower_N_100[1],N_100_interac_5$tsls_xm_upper_N_100[1]),c(N_100_interac_5$mean_betas.data.x_coeff_m[1],N_100_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_5$tsls_xm_lower_N_100[2],N_100_interac_5$tsls_xm_upper_N_100[2]),c(N_100_interac_5$mean_betas.data.x_coeff_m[2],N_100_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_5$tsls_xm_lower_N_100[3],N_100_interac_5$tsls_xm_upper_N_100[3]),c(N_100_interac_5$mean_betas.data.x_coeff_m[3],N_100_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_5$tsls_xm_lower_N_100[4],N_100_interac_5$tsls_xm_upper_N_100[4]),c(N_100_interac_5$mean_betas.data.x_coeff_m[4],N_100_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_5$tsls_xm_lower_N_100[5],N_100_interac_5$tsls_xm_upper_N_100[5]),c(N_100_interac_5$mean_betas.data.x_coeff_m[5],N_100_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_1$mean_betas.data.xm_2sls,N_100_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.6,0.5),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_100_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.333,0.4,0.5),labels=c('-0.6','-0.5','-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.333','0.4','0.5'),cex.axis=0.25)
lines(c(N_100_interac_1$tsls_xm_lower_N_100[1],N_100_interac_1$tsls_xm_upper_N_100[1]),c(N_100_interac_1$mean_betas.data.x_coeff_m[1],N_100_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_1$tsls_xm_lower_N_100[2],N_100_interac_1$tsls_xm_upper_N_100[2]),c(N_100_interac_1$mean_betas.data.x_coeff_m[2],N_100_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_1$tsls_xm_lower_N_100[3],N_100_interac_1$tsls_xm_upper_N_100[3]),c(N_100_interac_1$mean_betas.data.x_coeff_m[3],N_100_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_1$tsls_xm_lower_N_100[4],N_100_interac_1$tsls_xm_upper_N_100[4]),c(N_100_interac_1$mean_betas.data.x_coeff_m[4],N_100_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_1$tsls_xm_lower_N_100[5],N_100_interac_1$tsls_xm_upper_N_100[5]),c(N_100_interac_1$mean_betas.data.x_coeff_m[5],N_100_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()

#2SLS Z=(Z1,Z2,Z1Z2,Z1Z1) N=100K
tiff('coef_plot_z1z1_100k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_100_interac_m3$mean_betas.data.xm_z1_2sls,N_100_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.34,0.36),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='2sls estimate of theta coefficient, Z=(Z1,Z2,Z1Z2,Z1Z1)',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_100_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.4,-0.3,-0.2,-0.111,-0.1,0,0.1,0.2,0.3,0.4),labels=c('-0.4','-0.3','-0.2','-0.111','-0.1','0','0.1','0.2','0.3','0.4'),cex.axis=0.25)
lines(c(N_100_interac_m3$tsls_z1_xm_lower_N_100[1],N_100_interac_m3$tsls_z1_xm_upper_N_100[1]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[1],N_100_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_m3$tsls_z1_xm_lower_N_100[2],N_100_interac_m3$tsls_z1_xm_upper_N_100[2]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[2],N_100_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_m3$tsls_z1_xm_lower_N_100[3],N_100_interac_m3$tsls_z1_xm_upper_N_100[3]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[3],N_100_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_m3$tsls_z1_xm_lower_N_100[4],N_100_interac_m3$tsls_z1_xm_upper_N_100[4]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[4],N_100_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_m3$tsls_z1_xm_lower_N_100[5],N_100_interac_m3$tsls_z1_xm_upper_N_100[5]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[5],N_100_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_0$mean_betas.data.xm_z1_2sls,N_100_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.34,0.36),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_100_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4),labels=c('-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4'),cex.axis=0.25)
lines(c(N_100_interac_0$tsls_z1_xm_lower_N_100[1],N_100_interac_0$tsls_z1_xm_upper_N_100[1]),c(N_100_interac_0$mean_betas.data.x_coeff_m[1],N_100_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_0$tsls_z1_xm_lower_N_100[2],N_100_interac_0$tsls_z1_xm_upper_N_100[2]),c(N_100_interac_0$mean_betas.data.x_coeff_m[2],N_100_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_0$tsls_z1_xm_lower_N_100[3],N_100_interac_0$tsls_z1_xm_upper_N_100[3]),c(N_100_interac_0$mean_betas.data.x_coeff_m[3],N_100_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_0$tsls_z1_xm_lower_N_100[4],N_100_interac_0$tsls_z1_xm_upper_N_100[4]),c(N_100_interac_0$mean_betas.data.x_coeff_m[4],N_100_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_0$tsls_z1_xm_lower_N_100[5],N_100_interac_0$tsls_z1_xm_upper_N_100[5]),c(N_100_interac_0$mean_betas.data.x_coeff_m[5],N_100_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_3$mean_betas.data.xm_z1_2sls,N_100_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.34,0.36),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_100_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.111,0.2,0.3,0.4),labels=c('-0.4','-0.3','-0.2','-0.1','0','0.1','0.111','0.2','0.3','0.4'),cex.axis=0.25)
lines(c(N_100_interac_3$tsls_z1_xm_lower_N_100[1],N_100_interac_3$tsls_z1_xm_upper_N_100[1]),c(N_100_interac_3$mean_betas.data.x_coeff_m[1],N_100_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_3$tsls_z1_xm_lower_N_100[2],N_100_interac_3$tsls_z1_xm_upper_N_100[2]),c(N_100_interac_3$mean_betas.data.x_coeff_m[2],N_100_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_3$tsls_z1_xm_lower_N_100[3],N_100_interac_3$tsls_z1_xm_upper_N_100[3]),c(N_100_interac_3$mean_betas.data.x_coeff_m[3],N_100_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_3$tsls_z1_xm_lower_N_100[4],N_100_interac_3$tsls_z1_xm_upper_N_100[4]),c(N_100_interac_3$mean_betas.data.x_coeff_m[4],N_100_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_3$tsls_z1_xm_lower_N_100[5],N_100_interac_3$tsls_z1_xm_upper_N_100[5]),c(N_100_interac_3$mean_betas.data.x_coeff_m[5],N_100_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_5$mean_betas.data.xm_z1_2sls,N_100_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.34,0.36),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_100_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.167,0.2,0.3,0.4),labels=c('-0.4','-0.3','-0.2','-0.1','0','0.1','0.167','0.2','0.3','0.4'),cex.axis=0.25)
lines(c(N_100_interac_5$tsls_z1_xm_lower_N_100[1],N_100_interac_5$tsls_z1_xm_upper_N_100[1]),c(N_100_interac_5$mean_betas.data.x_coeff_m[1],N_100_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_5$tsls_z1_xm_lower_N_100[2],N_100_interac_5$tsls_z1_xm_upper_N_100[2]),c(N_100_interac_5$mean_betas.data.x_coeff_m[2],N_100_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_5$tsls_z1_xm_lower_N_100[3],N_100_interac_5$tsls_z1_xm_upper_N_100[3]),c(N_100_interac_5$mean_betas.data.x_coeff_m[3],N_100_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_5$tsls_z1_xm_lower_N_100[4],N_100_interac_5$tsls_z1_xm_upper_N_100[4]),c(N_100_interac_5$mean_betas.data.x_coeff_m[4],N_100_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_5$tsls_z1_xm_lower_N_100[5],N_100_interac_5$tsls_z1_xm_upper_N_100[5]),c(N_100_interac_5$mean_betas.data.x_coeff_m[5],N_100_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_1$mean_betas.data.xm_z1_2sls,N_100_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.34,0.36),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_100_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.333,0.4),labels=c('-0.4','-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.333','0.4'),cex.axis=0.25)
lines(c(N_100_interac_1$tsls_z1_xm_lower_N_100[1],N_100_interac_1$tsls_z1_xm_upper_N_100[1]),c(N_100_interac_1$mean_betas.data.x_coeff_m[1],N_100_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_1$tsls_z1_xm_lower_N_100[2],N_100_interac_1$tsls_z1_xm_upper_N_100[2]),c(N_100_interac_1$mean_betas.data.x_coeff_m[2],N_100_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_1$tsls_z1_xm_lower_N_100[3],N_100_interac_1$tsls_z1_xm_upper_N_100[3]),c(N_100_interac_1$mean_betas.data.x_coeff_m[3],N_100_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_1$tsls_z1_xm_lower_N_100[4],N_100_interac_1$tsls_z1_xm_upper_N_100[4]),c(N_100_interac_1$mean_betas.data.x_coeff_m[4],N_100_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_1$tsls_z1_xm_lower_N_100[5],N_100_interac_1$tsls_z1_xm_upper_N_100[5]),c(N_100_interac_1$mean_betas.data.x_coeff_m[5],N_100_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()

#OBSERVATIONAL N=100K

tiff('coef_plot_obs_100k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_100_interac_m3$mean_betas.data.xm_obs,N_100_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='Ordinary least sqaures estimate of theta coefficient',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_100_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.111,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.111','-0.1','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_100_interac_m3$obs_xm_lower_N_100[1],N_100_interac_m3$obs_xm_upper_N_100[1]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[1],N_100_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_m3$obs_xm_lower_N_100[2],N_100_interac_m3$obs_xm_upper_N_100[2]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[2],N_100_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_m3$obs_xm_lower_N_100[3],N_100_interac_m3$obs_xm_upper_N_100[3]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[3],N_100_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_m3$obs_xm_lower_N_100[4],N_100_interac_m3$obs_xm_upper_N_100[4]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[4],N_100_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_m3$obs_xm_lower_N_100[5],N_100_interac_m3$obs_xm_upper_N_100[5]),c(N_100_interac_m3$mean_betas.data.x_coeff_m[5],N_100_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_0$mean_betas.data.xm_obs,N_100_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_100_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_100_interac_0$obs_xm_lower_N_100[1],N_100_interac_0$obs_xm_upper_N_100[1]),c(N_100_interac_0$mean_betas.data.x_coeff_m[1],N_100_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_0$obs_xm_lower_N_100[2],N_100_interac_0$obs_xm_upper_N_100[2]),c(N_100_interac_0$mean_betas.data.x_coeff_m[2],N_100_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_0$obs_xm_lower_N_100[3],N_100_interac_0$obs_xm_upper_N_100[3]),c(N_100_interac_0$mean_betas.data.x_coeff_m[3],N_100_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_0$obs_xm_lower_N_100[4],N_100_interac_0$obs_xm_upper_N_100[4]),c(N_100_interac_0$mean_betas.data.x_coeff_m[4],N_100_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_0$obs_xm_lower_N_100[5],N_100_interac_0$obs_xm_upper_N_100[5]),c(N_100_interac_0$mean_betas.data.x_coeff_m[5],N_100_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_3$mean_betas.data.xm_obs,N_100_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_100_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.111,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.111','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_100_interac_3$obs_xm_lower_N_100[1],N_100_interac_3$obs_xm_upper_N_100[1]),c(N_100_interac_3$mean_betas.data.x_coeff_m[1],N_100_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_3$obs_xm_lower_N_100[2],N_100_interac_3$obs_xm_upper_N_100[2]),c(N_100_interac_3$mean_betas.data.x_coeff_m[2],N_100_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_3$obs_xm_lower_N_100[3],N_100_interac_3$obs_xm_upper_N_100[3]),c(N_100_interac_3$mean_betas.data.x_coeff_m[3],N_100_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_3$obs_xm_lower_N_100[4],N_100_interac_3$obs_xm_upper_N_100[4]),c(N_100_interac_3$mean_betas.data.x_coeff_m[4],N_100_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_3$obs_xm_lower_N_100[5],N_100_interac_3$obs_xm_upper_N_100[5]),c(N_100_interac_3$mean_betas.data.x_coeff_m[5],N_100_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_5$mean_betas.data.xm_obs,N_100_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_100_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.167,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.167','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_100_interac_5$obs_xm_lower_N_100[1],N_100_interac_5$obs_xm_upper_N_100[1]),c(N_100_interac_5$mean_betas.data.x_coeff_m[1],N_100_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_5$obs_xm_lower_N_100[2],N_100_interac_5$obs_xm_upper_N_100[2]),c(N_100_interac_5$mean_betas.data.x_coeff_m[2],N_100_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_5$obs_xm_lower_N_100[3],N_100_interac_5$obs_xm_upper_N_100[3]),c(N_100_interac_5$mean_betas.data.x_coeff_m[3],N_100_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_5$obs_xm_lower_N_100[4],N_100_interac_5$obs_xm_upper_N_100[4]),c(N_100_interac_5$mean_betas.data.x_coeff_m[4],N_100_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_5$obs_xm_lower_N_100[5],N_100_interac_5$obs_xm_upper_N_100[5]),c(N_100_interac_5$mean_betas.data.x_coeff_m[5],N_100_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_100_interac_1$mean_betas.data.xm_obs,N_100_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_100_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.333,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.333','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_100_interac_1$obs_xm_lower_N_100[1],N_100_interac_1$obs_xm_upper_N_100[1]),c(N_100_interac_1$mean_betas.data.x_coeff_m[1],N_100_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_100_interac_1$obs_xm_lower_N_100[2],N_100_interac_1$obs_xm_upper_N_100[2]),c(N_100_interac_1$mean_betas.data.x_coeff_m[2],N_100_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_100_interac_1$obs_xm_lower_N_100[3],N_100_interac_1$obs_xm_upper_N_100[3]),c(N_100_interac_1$mean_betas.data.x_coeff_m[3],N_100_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_100_interac_1$obs_xm_lower_N_100[4],N_100_interac_1$obs_xm_upper_N_100[4]),c(N_100_interac_1$mean_betas.data.x_coeff_m[4],N_100_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_100_interac_1$obs_xm_lower_N_100[5],N_100_interac_1$obs_xm_upper_N_100[5]),c(N_100_interac_1$mean_betas.data.x_coeff_m[5],N_100_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()


#PREP THE N=500K DATA
setwd(dataloc)
N_500=data.frame(read.table(file=paste(500000,"_EXTRA_final_res.txt",sep=''),sep='\t',header=TRUE))
N_500=N_500[order(N_500$mean_betas.data.xm_coeff_y,N_500$mean_betas.data.x_coeff_m),]
tsls_xm_lower_N_500=N_500$mean_betas.data.xm_2sls-1.96*N_500$se.newdata.xm_2sls
tsls_xm_upper_N_500=N_500$mean_betas.data.xm_2sls+1.96*N_500$se.newdata.xm_2sls
tsls_z1_xm_lower_N_500=N_500$mean_betas.data.xm_z1_2sls-1.96*N_500$se.newdata.xm_z1_2sls
tsls_z1_xm_upper_N_500=N_500$mean_betas.data.xm_z1_2sls+1.96*N_500$se.newdata.xm_z1_2sls
obs_xm_lower_N_500=N_500$mean_betas.data.xm_obs-1.96*N_500$se.newdata.xm_obs
obs_xm_upper_N_500=N_500$mean_betas.data.xm_obs+1.96*N_500$se.newdata.xm_obs
N_500=data.frame(N_500,tsls_xm_lower_N_500,tsls_xm_upper_N_500,tsls_z1_xm_lower_N_500,tsls_z1_xm_upper_N_500,
                 obs_xm_lower_N_500,obs_xm_upper_N_500)
N_500_interac_m3=N_500[round(N_500$mean_betas.data.xm_coeff_y,3)==-0.111,]
N_500_interac_0=N_500[round(N_500$mean_betas.data.xm_coeff_y,3)==0,]
N_500_interac_3=N_500[round(N_500$mean_betas.data.xm_coeff_y,3)==0.111,]
N_500_interac_5=N_500[round(N_500$mean_betas.data.xm_coeff_y,3)==0.167,]
N_500_interac_1=N_500[round(N_500$mean_betas.data.xm_coeff_y,3)==0.333,]

#2SLS Z=(Z1,Z2,Z1Z2) N=500K
setwd(writeplot)
tiff('coef_plot_z1z2_500k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_500_interac_m3$mean_betas.data.xm_2sls,N_500_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='2sls estimate of theta coefficient, Z=(Z1,Z2,Z1Z2)',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_500_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.111,-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.111','-0.10','-0.05','0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_m3$tsls_xm_lower_N_500[1],N_500_interac_m3$tsls_xm_upper_N_500[1]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[1],N_500_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_m3$tsls_xm_lower_N_500[2],N_500_interac_m3$tsls_xm_upper_N_500[2]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[2],N_500_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_m3$tsls_xm_lower_N_500[3],N_500_interac_m3$tsls_xm_upper_N_500[3]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[3],N_500_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_m3$tsls_xm_lower_N_500[4],N_500_interac_m3$tsls_xm_upper_N_500[4]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[4],N_500_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_m3$tsls_xm_lower_N_500[5],N_500_interac_m3$tsls_xm_upper_N_500[5]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[5],N_500_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_0$mean_betas.data.xm_2sls,N_500_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_500_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_0$tsls_xm_lower_N_500[1],N_500_interac_0$tsls_xm_upper_N_500[1]),c(N_500_interac_0$mean_betas.data.x_coeff_m[1],N_500_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_0$tsls_xm_lower_N_500[2],N_500_interac_0$tsls_xm_upper_N_500[2]),c(N_500_interac_0$mean_betas.data.x_coeff_m[2],N_500_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_0$tsls_xm_lower_N_500[3],N_500_interac_0$tsls_xm_upper_N_500[3]),c(N_500_interac_0$mean_betas.data.x_coeff_m[3],N_500_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_0$tsls_xm_lower_N_500[4],N_500_interac_0$tsls_xm_upper_N_500[4]),c(N_500_interac_0$mean_betas.data.x_coeff_m[4],N_500_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_0$tsls_xm_lower_N_500[5],N_500_interac_0$tsls_xm_upper_N_500[5]),c(N_500_interac_0$mean_betas.data.x_coeff_m[5],N_500_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_3$mean_betas.data.xm_2sls,N_500_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_500_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.111,0.15,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.111','0.15','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_3$tsls_xm_lower_N_500[1],N_500_interac_3$tsls_xm_upper_N_500[1]),c(N_500_interac_3$mean_betas.data.x_coeff_m[1],N_500_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_3$tsls_xm_lower_N_500[2],N_500_interac_3$tsls_xm_upper_N_500[2]),c(N_500_interac_3$mean_betas.data.x_coeff_m[2],N_500_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_3$tsls_xm_lower_N_500[3],N_500_interac_3$tsls_xm_upper_N_500[3]),c(N_500_interac_3$mean_betas.data.x_coeff_m[3],N_500_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_3$tsls_xm_lower_N_500[4],N_500_interac_3$tsls_xm_upper_N_500[4]),c(N_500_interac_3$mean_betas.data.x_coeff_m[4],N_500_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_3$tsls_xm_lower_N_500[5],N_500_interac_3$tsls_xm_upper_N_500[5]),c(N_500_interac_3$mean_betas.data.x_coeff_m[5],N_500_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_5$mean_betas.data.xm_2sls,N_500_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_500_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.15,0.167,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.167','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_5$tsls_xm_lower_N_500[1],N_500_interac_5$tsls_xm_upper_N_500[1]),c(N_500_interac_5$mean_betas.data.x_coeff_m[1],N_500_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_5$tsls_xm_lower_N_500[2],N_500_interac_5$tsls_xm_upper_N_500[2]),c(N_500_interac_5$mean_betas.data.x_coeff_m[2],N_500_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_5$tsls_xm_lower_N_500[3],N_500_interac_5$tsls_xm_upper_N_500[3]),c(N_500_interac_5$mean_betas.data.x_coeff_m[3],N_500_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_5$tsls_xm_lower_N_500[4],N_500_interac_5$tsls_xm_upper_N_500[4]),c(N_500_interac_5$mean_betas.data.x_coeff_m[4],N_500_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_5$tsls_xm_lower_N_500[5],N_500_interac_5$tsls_xm_upper_N_500[5]),c(N_500_interac_5$mean_betas.data.x_coeff_m[5],N_500_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_1$mean_betas.data.xm_2sls,N_500_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_500_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.333,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.20','0.25','0.30','0.333','0.35'),cex.axis=0.25)
lines(c(N_500_interac_1$tsls_xm_lower_N_500[1],N_500_interac_1$tsls_xm_upper_N_500[1]),c(N_500_interac_1$mean_betas.data.x_coeff_m[1],N_500_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_1$tsls_xm_lower_N_500[2],N_500_interac_1$tsls_xm_upper_N_500[2]),c(N_500_interac_1$mean_betas.data.x_coeff_m[2],N_500_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_1$tsls_xm_lower_N_500[3],N_500_interac_1$tsls_xm_upper_N_500[3]),c(N_500_interac_1$mean_betas.data.x_coeff_m[3],N_500_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_1$tsls_xm_lower_N_500[4],N_500_interac_1$tsls_xm_upper_N_500[4]),c(N_500_interac_1$mean_betas.data.x_coeff_m[4],N_500_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_1$tsls_xm_lower_N_500[5],N_500_interac_1$tsls_xm_upper_N_500[5]),c(N_500_interac_1$mean_betas.data.x_coeff_m[5],N_500_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()

#2SLS Z=(Z1,Z2,Z1Z2,Z1Z1) N=500K
tiff('coef_plot_z1z1_500k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_500_interac_m3$mean_betas.data.xm_z1_2sls,N_500_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='2sls estimate of theta coefficient, Z=(Z1,Z2,Z1Z2,Z1Z1)',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_500_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.111,-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.111','-0.10','-0.05','0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_m3$tsls_z1_xm_lower_N_500[1],N_500_interac_m3$tsls_z1_xm_upper_N_500[1]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[1],N_500_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_m3$tsls_z1_xm_lower_N_500[2],N_500_interac_m3$tsls_z1_xm_upper_N_500[2]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[2],N_500_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_m3$tsls_z1_xm_lower_N_500[3],N_500_interac_m3$tsls_z1_xm_upper_N_500[3]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[3],N_500_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_m3$tsls_z1_xm_lower_N_500[4],N_500_interac_m3$tsls_z1_xm_upper_N_500[4]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[4],N_500_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_m3$tsls_z1_xm_lower_N_500[5],N_500_interac_m3$tsls_z1_xm_upper_N_500[5]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[5],N_500_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_0$mean_betas.data.xm_z1_2sls,N_500_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_500_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_0$tsls_z1_xm_lower_N_500[1],N_500_interac_0$tsls_z1_xm_upper_N_500[1]),c(N_500_interac_0$mean_betas.data.x_coeff_m[1],N_500_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_0$tsls_z1_xm_lower_N_500[2],N_500_interac_0$tsls_z1_xm_upper_N_500[2]),c(N_500_interac_0$mean_betas.data.x_coeff_m[2],N_500_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_0$tsls_z1_xm_lower_N_500[3],N_500_interac_0$tsls_z1_xm_upper_N_500[3]),c(N_500_interac_0$mean_betas.data.x_coeff_m[3],N_500_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_0$tsls_z1_xm_lower_N_500[4],N_500_interac_0$tsls_z1_xm_upper_N_500[4]),c(N_500_interac_0$mean_betas.data.x_coeff_m[4],N_500_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_0$tsls_z1_xm_lower_N_500[5],N_500_interac_0$tsls_z1_xm_upper_N_500[5]),c(N_500_interac_0$mean_betas.data.x_coeff_m[5],N_500_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_3$mean_betas.data.xm_z1_2sls,N_500_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_500_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.111,0.15,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.111','0.15','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_3$tsls_z1_xm_lower_N_500[1],N_500_interac_3$tsls_z1_xm_upper_N_500[1]),c(N_500_interac_3$mean_betas.data.x_coeff_m[1],N_500_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_3$tsls_z1_xm_lower_N_500[2],N_500_interac_3$tsls_z1_xm_upper_N_500[2]),c(N_500_interac_3$mean_betas.data.x_coeff_m[2],N_500_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_3$tsls_z1_xm_lower_N_500[3],N_500_interac_3$tsls_z1_xm_upper_N_500[3]),c(N_500_interac_3$mean_betas.data.x_coeff_m[3],N_500_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_3$tsls_z1_xm_lower_N_500[4],N_500_interac_3$tsls_z1_xm_upper_N_500[4]),c(N_500_interac_3$mean_betas.data.x_coeff_m[4],N_500_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_3$tsls_z1_xm_lower_N_500[5],N_500_interac_3$tsls_z1_xm_upper_N_500[5]),c(N_500_interac_3$mean_betas.data.x_coeff_m[5],N_500_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_5$mean_betas.data.xm_z1_2sls,N_500_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_500_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.15,0.167,0.20,0.25,0.30,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.167','0.20','0.25','0.30','0.35'),cex.axis=0.25)
lines(c(N_500_interac_5$tsls_z1_xm_lower_N_500[1],N_500_interac_5$tsls_z1_xm_upper_N_500[1]),c(N_500_interac_5$mean_betas.data.x_coeff_m[1],N_500_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_5$tsls_z1_xm_lower_N_500[2],N_500_interac_5$tsls_z1_xm_upper_N_500[2]),c(N_500_interac_5$mean_betas.data.x_coeff_m[2],N_500_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_5$tsls_z1_xm_lower_N_500[3],N_500_interac_5$tsls_z1_xm_upper_N_500[3]),c(N_500_interac_5$mean_betas.data.x_coeff_m[3],N_500_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_5$tsls_z1_xm_lower_N_500[4],N_500_interac_5$tsls_z1_xm_upper_N_500[4]),c(N_500_interac_5$mean_betas.data.x_coeff_m[4],N_500_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_5$tsls_z1_xm_lower_N_500[5],N_500_interac_5$tsls_z1_xm_upper_N_500[5]),c(N_500_interac_5$mean_betas.data.x_coeff_m[5],N_500_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_1$mean_betas.data.xm_z1_2sls,N_500_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.15,0.35),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_500_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.15,-0.10,-0.05,0,0.05,0.10,0.15,0.20,0.25,0.30,0.333,0.35),labels=c('-0.15','-0.10','-0.05','0','0.05','0.10','0.15','0.20','0.25','0.30','0.333','0.35'),cex.axis=0.25)
lines(c(N_500_interac_1$tsls_z1_xm_lower_N_500[1],N_500_interac_1$tsls_z1_xm_upper_N_500[1]),c(N_500_interac_1$mean_betas.data.x_coeff_m[1],N_500_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_1$tsls_z1_xm_lower_N_500[2],N_500_interac_1$tsls_z1_xm_upper_N_500[2]),c(N_500_interac_1$mean_betas.data.x_coeff_m[2],N_500_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_1$tsls_z1_xm_lower_N_500[3],N_500_interac_1$tsls_z1_xm_upper_N_500[3]),c(N_500_interac_1$mean_betas.data.x_coeff_m[3],N_500_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_1$tsls_z1_xm_lower_N_500[4],N_500_interac_1$tsls_z1_xm_upper_N_500[4]),c(N_500_interac_1$mean_betas.data.x_coeff_m[4],N_500_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_1$tsls_z1_xm_lower_N_500[5],N_500_interac_1$tsls_z1_xm_upper_N_500[5]),c(N_500_interac_1$mean_betas.data.x_coeff_m[5],N_500_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()


#OBSERVATIONAL N=500K
tiff('coef_plot_obs_500k.tif',width=3.5,height=10,units='in',res=400)
par(mfrow=c(5,1),mar=c(4,4.1,2,2.1))
plot(N_500_interac_m3$mean_betas.data.xm_obs,N_500_interac_m3$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=-0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n',
     main='Ordinary least sqaures estimate of theta coefficient',cex.main=0.7)
abline(v=-0.111)
axis(side=2,at=N_500_interac_m3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.111,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.111','-0.1','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_500_interac_m3$obs_xm_lower_N_500[1],N_500_interac_m3$obs_xm_upper_N_500[1]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[1],N_500_interac_m3$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_m3$obs_xm_lower_N_500[2],N_500_interac_m3$obs_xm_upper_N_500[2]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[2],N_500_interac_m3$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_m3$obs_xm_lower_N_500[3],N_500_interac_m3$obs_xm_upper_N_500[3]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[3],N_500_interac_m3$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_m3$obs_xm_lower_N_500[4],N_500_interac_m3$obs_xm_upper_N_500[4]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[4],N_500_interac_m3$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_m3$obs_xm_lower_N_500[5],N_500_interac_m3$obs_xm_upper_N_500[5]),c(N_500_interac_m3$mean_betas.data.x_coeff_m[5],N_500_interac_m3$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_0$mean_betas.data.xm_obs,N_500_interac_0$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0)
axis(side=2,at=N_500_interac_0$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_500_interac_0$obs_xm_lower_N_500[1],N_500_interac_0$obs_xm_upper_N_500[1]),c(N_500_interac_0$mean_betas.data.x_coeff_m[1],N_500_interac_0$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_0$obs_xm_lower_N_500[2],N_500_interac_0$obs_xm_upper_N_500[2]),c(N_500_interac_0$mean_betas.data.x_coeff_m[2],N_500_interac_0$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_0$obs_xm_lower_N_500[3],N_500_interac_0$obs_xm_upper_N_500[3]),c(N_500_interac_0$mean_betas.data.x_coeff_m[3],N_500_interac_0$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_0$obs_xm_lower_N_500[4],N_500_interac_0$obs_xm_upper_N_500[4]),c(N_500_interac_0$mean_betas.data.x_coeff_m[4],N_500_interac_0$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_0$obs_xm_lower_N_500[5],N_500_interac_0$obs_xm_upper_N_500[5]),c(N_500_interac_0$mean_betas.data.x_coeff_m[5],N_500_interac_0$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_3$mean_betas.data.xm_obs,N_500_interac_3$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.111',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.111)
axis(side=2,at=N_500_interac_3$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.111,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.111','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_500_interac_3$obs_xm_lower_N_500[1],N_500_interac_3$obs_xm_upper_N_500[1]),c(N_500_interac_3$mean_betas.data.x_coeff_m[1],N_500_interac_3$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_3$obs_xm_lower_N_500[2],N_500_interac_3$obs_xm_upper_N_500[2]),c(N_500_interac_3$mean_betas.data.x_coeff_m[2],N_500_interac_3$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_3$obs_xm_lower_N_500[3],N_500_interac_3$obs_xm_upper_N_500[3]),c(N_500_interac_3$mean_betas.data.x_coeff_m[3],N_500_interac_3$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_3$obs_xm_lower_N_500[4],N_500_interac_3$obs_xm_upper_N_500[4]),c(N_500_interac_3$mean_betas.data.x_coeff_m[4],N_500_interac_3$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_3$obs_xm_lower_N_500[5],N_500_interac_3$obs_xm_upper_N_500[5]),c(N_500_interac_3$mean_betas.data.x_coeff_m[5],N_500_interac_3$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_5$mean_betas.data.xm_obs,N_500_interac_5$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.167',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.167)
axis(side=2,at=N_500_interac_5$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.167,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.167','0.2','0.3','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_500_interac_5$obs_xm_lower_N_500[1],N_500_interac_5$obs_xm_upper_N_500[1]),c(N_500_interac_5$mean_betas.data.x_coeff_m[1],N_500_interac_5$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_5$obs_xm_lower_N_500[2],N_500_interac_5$obs_xm_upper_N_500[2]),c(N_500_interac_5$mean_betas.data.x_coeff_m[2],N_500_interac_5$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_5$obs_xm_lower_N_500[3],N_500_interac_5$obs_xm_upper_N_500[3]),c(N_500_interac_5$mean_betas.data.x_coeff_m[3],N_500_interac_5$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_5$obs_xm_lower_N_500[4],N_500_interac_5$obs_xm_upper_N_500[4]),c(N_500_interac_5$mean_betas.data.x_coeff_m[4],N_500_interac_5$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_5$obs_xm_lower_N_500[5],N_500_interac_5$obs_xm_upper_N_500[5]),c(N_500_interac_5$mean_betas.data.x_coeff_m[5],N_500_interac_5$mean_betas.data.x_coeff_m[5]))
plot(N_500_interac_1$mean_betas.data.xm_obs,N_500_interac_1$mean_betas.data.x_coeff_m,xlim=c(-0.3,0.8),cex=0.5,
     xlab='theta=0.333',ylab='alpha coefficient',cex.lab=0.6,yaxt='n',xaxt='n')
abline(v=0.333)
axis(side=2,at=N_500_interac_1$mean_betas.data.x_coeff_m,labels=c('-0.333','0','0.333','0.5','1'),cex.axis=0.4)
axis(side=1,at=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.333,0.4,0.5,0.6,0.7,0.8),labels=c('-0.3','-0.2','-0.1','0','0.1','0.2','0.3','0.333','0.4','0.5','0.6','0.7','0.8'),cex.axis=0.25)
lines(c(N_500_interac_1$obs_xm_lower_N_500[1],N_500_interac_1$obs_xm_upper_N_500[1]),c(N_500_interac_1$mean_betas.data.x_coeff_m[1],N_500_interac_1$mean_betas.data.x_coeff_m[1]))
lines(c(N_500_interac_1$obs_xm_lower_N_500[2],N_500_interac_1$obs_xm_upper_N_500[2]),c(N_500_interac_1$mean_betas.data.x_coeff_m[2],N_500_interac_1$mean_betas.data.x_coeff_m[2]))
lines(c(N_500_interac_1$obs_xm_lower_N_500[3],N_500_interac_1$obs_xm_upper_N_500[3]),c(N_500_interac_1$mean_betas.data.x_coeff_m[3],N_500_interac_1$mean_betas.data.x_coeff_m[3]))
lines(c(N_500_interac_1$obs_xm_lower_N_500[4],N_500_interac_1$obs_xm_upper_N_500[4]),c(N_500_interac_1$mean_betas.data.x_coeff_m[4],N_500_interac_1$mean_betas.data.x_coeff_m[4]))
lines(c(N_500_interac_1$obs_xm_lower_N_500[5],N_500_interac_1$obs_xm_upper_N_500[5]),c(N_500_interac_1$mean_betas.data.x_coeff_m[5],N_500_interac_1$mean_betas.data.x_coeff_m[5]))
dev.off()

#POWER PLOTS
setwd(writeplot)

tiff('50_forestplot_power.tif',width=12.7,height=9,units='cm',res=400)
#layout(matrix(c(1,2)),widths=c(lcm(4.5),2),heights=c(lcm(12.7),1),byrow=FALSE)
layout(matrix(c(1,2),1,2,byrow=FALSE))
#vals=par('usr')
#rem_0interac=N_50[2:26,]
rem_0interac=N_50
rem_0interac=rem_0interac[rem_0interac$mean_betas.data.xm_coeff_y!=0,]
powervals=rbind(rem_0interac$xm_z1_2sls_detec.xm_z1_2sls_detec,rem_0interac$z1_coverage.z1_coverage)
powervals=as.matrix(powervals)
powervals=powervals/10
acolnames=cbind(round(rem_0interac$mean_betas.data.x_coeff_m,3),round(rem_0interac$mean_betas.data.xm_coeff_y,3))
colnamesconcat_1=paste('alpha:',acolnames[,1],sep='')
colnamesconcat_2=paste('theta:',acolnames[,2],sep='')
colnamesconcat=paste(colnamesconcat_1,colnamesconcat_2,sep=', ')
par(las=1)
spacevec1=c(1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
spacevec2=c(1.5,
           0.5,0.5,0.5,0.5,
           1.5,
           0.5,0.5,0.5,0.5,
           1.5,
           0.5,0.5,0.5,0.5,
           1.5,
           0.5,0.5,0.5,0.5)
colnames(powervals)=colnamesconcat
fmrpowervals=rem_0interac$fmr_y_int.fmr_y_int
fmrpowervals=fmrpowervals/10
#widthspec=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1)
barplot(powervals,horiz=TRUE,cex.names=0.3,cex.axis=0.4,xlab='Power and Coverage of 2sls Estimator of theta (%)',cex.lab=0.4,
        xlim=c(0,100),main="2sls with Z=(Z1,Z2,Z1Z2,Z1Z1)", cex.main="0.5",
        beside=TRUE,col=c('blue','grey'),space=spacevec1)
        #width=widthspec
        #space=spacevec
        #names.arg=colnamesconcat
#text(x=91,y=66,labels="95%",cex=0.3)
mtext("95%",side=1,at=95,cex=0.3)
abline(v=80,lty=2,col='blue')
abline(v=95,lty=2,col='black')
barplot(fmrpowervals,horiz=TRUE,names.arg=colnamesconcat,cex.names=0.3,cex.axis=0.4,xlab='Power of FMR to detect non-zero theta (%)',cex.lab=0.4,col='blue',
        space=spacevec2,xlim=c(0,100),main="Factorial MR",cex.main="0.5")
abline(v=80,lty=2,col='blue')
legend(x='bottomleft',inset=c(-0.6,-0.55),legend=c('power','coverage'),col=c('blue','grey'),pch=19,xpd=TRUE,cex=0.4)
dev.off()



tiff('100_forestplot_power.tif',width=12.7,height=9,units='cm',res=400)
#layout(matrix(c(1,2)),widths=c(lcm(4.5),2),heights=c(lcm(12.7),1),byrow=FALSE)
layout(matrix(c(1,2),1,2,byrow=FALSE))
#vals=par('usr')
#crem_0interac=N_100[2:26,]
crem_0interac=N_100
crem_0interac=crem_0interac[crem_0interac$mean_betas.data.xm_coeff_y!=0,]
cpowervals=rbind(crem_0interac$xm_z1_2sls_detec.xm_z1_2sls_detec,crem_0interac$z1_coverage.z1_coverage)
cpowervals=as.matrix(cpowervals)
cpowervals=cpowervals/10
ccolnames=cbind(round(crem_0interac$mean_betas.data.x_coeff_m,3),round(crem_0interac$mean_betas.data.xm_coeff_y,3))
ccolnamesconcat_1=paste('alpha:',ccolnames[,1],sep='')
ccolnamesconcat_2=paste('theta:',ccolnames[,2],sep='')
ccolnamesconcat=paste(ccolnamesconcat_1,ccolnamesconcat_2,sep=', ')
par(las=1)
cspacevec1=c(1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
cspacevec2=c(1.5,
            0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5)
colnames(cpowervals)=ccolnamesconcat
cfmrpowervals=crem_0interac$fmr_y_int.fmr_y_int
cfmrpowervals=cfmrpowervals/10
#widthspec=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1)
barplot(cpowervals,horiz=TRUE,cex.names=0.3,cex.axis=0.4,xlab='Power and Coverage of 2sls Estimator of theta (%)',cex.lab=0.4,
        xlim=c(0,100),main="2sls with Z=(Z1,Z2,Z1Z2,Z1Z1)", cex.main="0.5",
        beside=TRUE,col=c('blue','grey'),space=cspacevec1)
#width=widthspec
#space=spacevec
#names.arg=colnamesconcat
#text(x=91,y=66,labels="95%",cex=0.3)
mtext("95%",side=1,at=95,cex=0.3)
abline(v=80,lty=2,col='blue')
abline(v=95,lty=2,col='black')
barplot(cfmrpowervals,horiz=TRUE,names.arg=ccolnamesconcat,cex.names=0.3,cex.axis=0.4,xlab='Power of FMR to detect non-zero theta (%)',cex.lab=0.4,col='blue',
        space=cspacevec2,xlim=c(0,100),main="Factorial MR",cex.main="0.5")
abline(v=80,lty=2,col='blue')
legend(x='bottomleft',inset=c(-0.6,-0.55),legend=c('power','coverage'),col=c('blue','grey'),pch=19,xpd=TRUE,cex=0.4)
dev.off()





tiff('500_forestplot_power.tif',width=12.7,height=9,units='cm',res=400)
layout(matrix(c(1,2),1,2,byrow=FALSE))
#vals=par('usr')
#brem_0interac=N_500[2:26,]
brem_0interac=N_500
brem_0interac=brem_0interac[brem_0interac$mean_betas.data.xm_coeff_y!=0,]
bpowervals=rbind(brem_0interac$xm_z1_2sls_detec.xm_z1_2sls_detec,brem_0interac$z1_coverage.z1_coverage)
bpowervals=as.matrix(bpowervals)
bpowervals=bpowervals/10
bcolnames=cbind(round(brem_0interac$mean_betas.data.x_coeff_m,3),round(brem_0interac$mean_betas.data.xm_coeff_y,3))
bcolnamesconcat_1=paste('alpha:',bcolnames[,1],sep='')
bcolnamesconcat_2=paste('theta:',bcolnames[,2],sep='')
bcolnamesconcat=paste(bcolnamesconcat_1,bcolnamesconcat_2,sep=', ')
par(las=1)
bspacevec1=c(1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
            1.5,
            0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
bspacevec2=c(1.5,
           0.5,0.5,0.5,0.5,
           1.5,
           0.5,0.5,0.5,0.5,
           1.5,
           0.5,0.5,0.5,0.5,
           1.5,
           0.5,0.5,0.5,0.5)
colnames(bpowervals)=bcolnamesconcat
bfmrpowervals=brem_0interac$fmr_y_int.fmr_y_int
bfmrpowervals=bfmrpowervals/10
barplot(bpowervals,horiz=TRUE,cex.names=0.3,cex.axis=0.4,xlab='Power and Coverage of 2sls Estimator of theta (%)',cex.lab=0.4,
        xlim=c(0,100),main="2sls with Z=(Z1,Z2,Z1Z2,Z1Z1)", cex.main="0.5",
        beside=TRUE,col=c('blue','grey'),space=bspacevec1)
#text(x=91,y=65,labels="95%",cex=0.3)
#text(x=91,y=-1,labels="95%",cex=0.3)
mtext("95%",side=1,at=95,cex=0.3)
abline(v=80,lty=2,col='blue')
abline(v=95,lty=2,col='black')
barplot(bfmrpowervals,horiz=TRUE,names.arg=bcolnamesconcat,cex.names=0.3,cex.axis=0.4,xlab='Power of FMR to detect non-zero theta (%)',cex.lab=0.4,col='blue',
        space=bspacevec2,xlim=c(0,100),main="Factorial MR",cex.main="0.5")
abline(v=80,lty=2,col='blue')
legend(x='bottomleft',inset=c(-0.6,-0.55),legend=c('power','coverage'),col=c('blue','grey'),pch=19,xpd=TRUE,cex=0.4)
dev.off()



#TYPE I ERROR PLOTS

tiff('50_forestplot_typei.tif',width=12.7,height=9,units='cm',res=400)
#vals=par('usr')
#keep_0interac=N_50[2:26,]
keep_0interac=N_50
keep_0interac=keep_0interac[keep_0interac$mean_betas.data.xm_coeff_y==0,]
typeivals=rbind(keep_0interac$xm_z1_2sls_detec.xm_z1_2sls_detec,keep_0interac$fmr_y_int.fmr_y_int)
typeivals=as.matrix(typeivals)
typeivals=typeivals/10
ticolnames=cbind(round(keep_0interac$mean_betas.data.x_coeff_m,3),round(keep_0interac$mean_betas.data.xm_coeff_y,3))
ticolnamesconcat_1=paste('alpha:',ticolnames[,1],sep='')
ticolnamesconcat_2=paste('theta:',ticolnames[,2],sep='')
ticolnamesconcat=paste(ticolnamesconcat_1,ticolnamesconcat_2,sep=', ')
par(las=1)
colnames(typeivals)=ticolnamesconcat
#widthspec=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1)
barplot(typeivals,horiz=TRUE,cex.names=0.3,cex.axis=0.4,xlab='Type I Error Rate',cex.lab=0.4,
        xlim=c(0,6),main='2sls vs FMR: N=50,000', cex.main="0.5",
        beside=TRUE,col=c('blue','red'))
#space=spacevec1
#width=widthspec
#space=spacevec
#names.arg=colnamesconcat
abline(v=5,lty=2,col='black')
legend(x='bottomright',legend=c('2sls','FMR'),col=c('blue','red'),pch=19,xpd=TRUE,cex=0.4)
dev.off()



tiff('100_forestplot_typei.tif',width=12.7,height=9,units='cm',res=400)
#vals=par('usr')
#ckeep_0interac=N_100[2:26,]
ckeep_0interac=N_100
ckeep_0interac=ckeep_0interac[ckeep_0interac$mean_betas.data.xm_coeff_y==0,]
ctypeivals=rbind(ckeep_0interac$xm_z1_2sls_detec.xm_z1_2sls_detec,ckeep_0interac$fmr_y_int.fmr_y_int)
ctypeivals=as.matrix(ctypeivals)
ctypeivals=ctypeivals/10
cticolnames=cbind(round(ckeep_0interac$mean_betas.data.x_coeff_m,3),round(ckeep_0interac$mean_betas.data.xm_coeff_y,3))
cticolnamesconcat_1=paste('alpha:',cticolnames[,1],sep='')
cticolnamesconcat_2=paste('theta:',cticolnames[,2],sep='')
cticolnamesconcat=paste(cticolnamesconcat_1,cticolnamesconcat_2,sep=', ')
par(las=1)
colnames(ctypeivals)=cticolnamesconcat
#widthspec=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1)
barplot(ctypeivals,horiz=TRUE,cex.names=0.3,cex.axis=0.4,xlab='Type I Error Rate',cex.lab=0.4,
        xlim=c(0,6),main='2sls vs FMR: N=100,000', cex.main="0.5",
        beside=TRUE,col=c('blue','red'))
#space=spacevec1
#width=widthspec
#space=spacevec
#names.arg=colnamesconcat
abline(v=5,lty=2,col='black')
legend(x='bottomright',legend=c('2sls','FMR'),col=c('blue','red'),pch=19,xpd=TRUE,cex=0.4)
dev.off()





tiff('500_forestplot_typei.tif',width=12.7,height=9,units='cm',res=400)
#vals=par('usr')
#bkeep_0interac=N_500[2:26,]
bkeep_0interac=N_500
bkeep_0interac=bkeep_0interac[bkeep_0interac$mean_betas.data.xm_coeff_y==0,]
btypeivals=rbind(bkeep_0interac$xm_z1_2sls_detec.xm_z1_2sls_detec,bkeep_0interac$fmr_y_int.fmr_y_int)
btypeivals=as.matrix(btypeivals)
btypeivals=btypeivals/10
bticolnames=cbind(round(bkeep_0interac$mean_betas.data.x_coeff_m,3),round(bkeep_0interac$mean_betas.data.xm_coeff_y,3))
bticolnamesconcat_1=paste('alpha:',bticolnames[,1],sep='')
bticolnamesconcat_2=paste('theta:',bticolnames[,2],sep='')
bticolnamesconcat=paste(bticolnamesconcat_1,bticolnamesconcat_2,sep=', ')
par(las=1)
colnames(btypeivals)=bticolnamesconcat
#widthspec=c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,1)
barplot(btypeivals,horiz=TRUE,cex.names=0.3,cex.axis=0.4,xlab='Type I Error Rate',cex.lab=0.4,
        xlim=c(0,6),main='2sls vs FMR: N=500,000', cex.main="0.5",
        beside=TRUE,col=c('blue','red'))
#space=spacevec1
#width=widthspec
#space=spacevec
#names.arg=colnamesconcat
abline(v=5,lty=2,col='black')
legend(x='bottomright',legend=c('2sls','FMR'),col=c('blue','red'),pch=19,xpd=TRUE,cex=0.4)
dev.off()


