library(mvtnorm)

#Example 1
EX1 <- read.csv("~/PhD/Infering_idividual_causal/SIM.csv")


Fig<-c()
Fig1<-c()
Fig2<-c()

col<-0

for(i in c(1:1000)){
    col<-col+1
    data<-EX1[EX1$ID==i,]
    #Simulation parameters
    beta0<-120
    beta1<- -10
    beta2<- -5
    betaC<-5
    
    sigma0<-1
    tau0<-5
    tau1<-10
    tau2<-5
    
    ind<-0
    for(s in c(3,10,100)){
        print(c(i,s))
        ind<-ind+1
        #<-100
        #s<-length(data[,1])
        
        #Simulated data
        Y<-data$Y_lag_0[1:s]
        C<-data$C_lag_1[1:s]
        A1<-data$A_lag_1[1:s]
        A2<-data$A_lag_2[1:s]
        
        k<-3
        Y.k<-Y[k]
        A.1<-A1[k]
        A.2<-A2[k]
        
        mu.u<-rep(0,3)
        mu.y<-beta0+C*betaC+A1*beta1+A2*beta2
        
        
        sigma.11<-diag(c(tau0^2,tau1^2,tau2^2))
        sigma.12<-as.matrix(rbind(rep(tau0^2,length(A1)),A1^2*tau1^2,A2^2*tau2^2))
        sigma.22<- diag(rep(sigma0,length(A1)))  +tau0^2 +tau1^2*A1%*%t(A1)+tau2^2*A2%*%t(A2)
        
        
        mu<-sigma.12%*%solve(sigma.22)%*%as.matrix(Y-mu.y)
        
        # print(cbind(mu,unique(c(data$u0,data$u1,data$u2))))
        sigma<-sigma.11 - sigma.12%*%solve(sigma.22)%*%t(sigma.12)
        
        sigma[2,3]<-sigma[3,2]
        
    
        
        
        #Real ICE    
        a.1<-1
        a.2<-1
        
        ICE<-a.1*(beta1+unique(data$u1))+a.2*(beta2+unique(data$u2))
        Y.ak<-Y.k+(a.1-A.1)*(beta1+unique(data$u1))+(a.2-A.2)*(beta2+unique(data$u2))
        Y.0<-Y.k+(0-A.1)*(beta1+unique(data$u1))+(0-A.2)*(beta2+unique(data$u2))
      
        
        mu.ICE<-a.1*(beta1)+a.2*(beta2)+matrix(c(0,(a.1),(a.2)),nrow=1)%*%mu
        sigma.ICE<-matrix(c(0,(a.1),(a.2)),nrow=1)%*%sigma%*%t(matrix(c(0,(a.1),(a.2)),nrow=1))
        
        
        Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,sum(A1),mu.ICE,sigma.ICE))
        
    }
}

Fig1<-as.data.frame(Fig1)
colnames(Fig1)<-c("ID","Ya","Y0","ICE","Repeats","TRTs","mu.ICE","sigma.ICE")

#Fig 8a
library(tikzDevice)

a<- -30
b<-5

tikz(file = "D:/Documents/PhD/Infering_idividual_causal/TikZFigures/GaussianCWC.tex", width = 7, height = 6)


plot(seq(Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),7]-3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),8]),Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),7]+3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),8]),0.01),
     sapply(seq(Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),7]-3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),8]),Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),7]+3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),8]),0.01),
       function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),7],sd=sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==1),8]))}),
     lwd=2,xlab="CWCE",ylab="pdf",type="l",xlim=c(a,b),ylim=c(0,0.75),col="red")

lines(seq(Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),7]-3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),8]),Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),7]+3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),8]),0.01),
     sapply(seq(Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),7]-3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),8]),Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),7]+3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),8]),0.01),
            function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),7],sd=sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==1),8]))}),
     lwd=2,col="red",lty=2)

lines(seq(Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),7]-3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),8]),Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),7]+3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),8]),0.01),
      sapply(seq(Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),7]-3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),8]),Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),7]+3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),8]),0.01),
             function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),7],sd=sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),8]))}),
      lwd=2,col="red",lty=3)

abline(v=Fig1[which(Fig1$Repeats==100 & Fig1$ID==1),4],col="red",lwd=2)

lines(seq(Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),7]-3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),8]),Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),7]+3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),8]),0.01),
      sapply(seq(Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),7]-3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),8]),Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),7]+3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),8]),0.01),
             function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),7],sd=sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==2),8]))}),
      lwd=2,col="blue",lty=3)

lines(seq(Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),7]-3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),8]),Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),7]+3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),8]),0.01),
      sapply(seq(Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),7]-3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),8]),Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),7]+3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),8]),0.01),
             function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),7],sd=sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==2),8]))}),
      lwd=2,col="blue",lty=2)

lines(seq(Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),7]-3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),8]),Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),7]+3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),8]),0.01),
      sapply(seq(Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),7]-3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),8]),Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),7]+3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),8]),0.01),
             function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),7],sd=sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),8]))}),
      lwd=2,col="blue",lty=1)

abline(v=Fig1[which(Fig1$Repeats==100 & Fig1$ID==2),4],col="blue",lwd=2)

lines(seq(Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),7]-3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),8]),Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),7]+3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),8]),0.01),
      sapply(seq(Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),7]-3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),8]),Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),7]+3*sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),8]),0.01),
             function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),7],sd=sqrt(Fig1[which(Fig1$Repeats==3 & Fig1$ID==15),8]))}),
      lwd=2,lty=3)

lines(seq(Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),7]-3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),8]),Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),7]+3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),8]),0.01),
      sapply(seq(Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),7]-3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),8]),Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),7]+3*sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),8]),0.01),
             function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),7],sd=sqrt(Fig1[which(Fig1$Repeats==10 & Fig1$ID==15),8]))}),
      lwd=2,lty=2)

lines(seq(Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),7]-3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),8]),Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),7]+3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),8]),0.01),
      sapply(seq(Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),7]-3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),8]),Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),7]+3*sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),8]),0.01),
             function(x){dnorm(x,mean=Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),7],sd=sqrt(Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),8]))}),
      lwd=2,lty=1)

abline(v=Fig1[which(Fig1$Repeats==100 & Fig1$ID==15),4],lwd=2)

legend(x=-12.5,y=0.7,legend=c(3,10,100),lty=c(3,2,1),col="grey",xpd=T,lwd=2,title="Repeats")
dev.off()

############################################
##########    Inference ####################
############################################

library(mvtnorm)

#Example 1
EX1 <- read.csv("~/PhD/Infering_idividual_causal/SIM.csv")

random <- read.csv("~/PhD/Infering_idividual_causal/EX1_cov.csv")
fixed <- read.csv("~/PhD/Infering_idividual_causal/EX1_fixed.csv")

design<-data.frame(sc=as.factor(c("A","B","C","D","E","F","G","H","I")),s=rep(c(3,10,100),3),n=rep(c(100,500,1000),each=3))
design$beta0<-sapply(design$sc,function(x){fixed$Estimate[which(fixed$scenario==x)[1]]})
design$beta1<-sapply(design$sc,function(x){fixed$Estimate[which(fixed$scenario==x)[3]]})
design$beta2<-sapply(design$sc,function(x){fixed$Estimate[which(fixed$scenario==x)[2]]})
design$betaC<-sapply(design$sc,function(x){fixed$Estimate[which(fixed$scenario==x)[4]]})

design$sigma0<-sapply(design$sc,function(x){random$Estimate[which(fixed$scenario==x)[4]]})
design$tau0<-sapply(design$sc,function(x){random$Estimate[which(fixed$scenario==x)[1]]})
design$tau1<-sapply(design$sc,function(x){random$Estimate[which(fixed$scenario==x)[3]]})
design$tau2<-sapply(design$sc,function(x){random$Estimate[which(fixed$scenario==x)[2]]})

Fig2<-NULL

for (sc in rev(c("A","B","C","D","E","F","G","H","I"))){
    
    ind<-which(design$sc==sc)
    
    n<-design$n[ind]
    s<-design$s[ind]
    
    for(i in c(1:n)){
        print(c(sc,i))
        data<-EX1[EX1$ID==i,]
        #Simulation parameters
        
        
        beta0<-design$beta0[ind]
        beta1<- design$beta1[ind]
        beta2<- design$beta2[ind]
        betaC<- design$betaC[ind]
        
        sigma0<- design$sigma0[ind]
        tau0<- design$tau0[ind]
        tau1<- design$tau1[ind]
        tau2<- design$tau2[ind]
        
            Y<-data$Y_lag_0[1:s]
            C<-data$C_lag_1[1:s]
            A1<-data$A_lag_1[1:s]
            A2<-data$A_lag_2[1:s]
            
            k<-3
            Y.k<-Y[k]
            A.1<-A1[k]
            A.2<-A2[k]
            
            mu.u<-rep(0,3)
            mu.y<-beta0+C*betaC+A1*beta1+A2*beta2
            
            
            sigma.11<-diag(c(tau0^2,tau1^2,tau2^2))
            sigma.12<-as.matrix(rbind(rep(tau0^2,length(A1)),A1^2*tau1^2,A2^2*tau2^2))
            sigma.22<- diag(rep(sigma0,length(A1)))  +tau0^2 +tau1^2*A1%*%t(A1)+tau2^2*A2%*%t(A2)
            
            
            mu<-sigma.12%*%solve(sigma.22)%*%as.matrix(Y-mu.y)
            
            # print(cbind(mu,unique(c(data$u0,data$u1,data$u2))))
            sigma<-sigma.11 - sigma.12%*%solve(sigma.22)%*%t(sigma.12)
            
            sigma[2,3]<-sigma[3,2]
            
            
            
            a.1<-1
            a.2<-1
            
            
            mu.ICE<-a.1*(beta1)+a.2*(beta2)+matrix(c(0,(a.1),(a.2)),nrow=1)%*%mu
            sigma.ICE<-matrix(c(0,(a.1),(a.2)),nrow=1)%*%sigma%*%t(matrix(c(0,(a.1),(a.2)),nrow=1))
            
            
            Fig2<-rbind(Fig2,c(sc,i,mu.ICE,sigma.ICE))
            
        }
}

Fig2<-as.data.frame(Fig2)
Fig2[,2]<-as.numeric(Fig2[,2])
Fig2[,3]<-as.numeric(Fig2[,3])
Fig2[,4]<-as.numeric(Fig2[,4])

colnames(Fig2)<-c("sc","ID","mu.ICE","sigma.ICE")

True.ICE<-unique(Fig1[,c(1,4)])

Fig3<-merge(Fig2,True.ICE)


library(RColorBrewer)
par(mfrow=c(3,3),mar=c(2,2,2,2))


#Plot ICE vs estimate E[ICE|H] 
col<-0
for (sc in c("A","B","C","D","E","F","G","H","I")){
    col<-col+1
    plot((Fig3[Fig3$sc==sc,c(5,3)]),lwd=2,col=brewer.pal(11, "Spectral")[c(1,3,4,5,7,8,9,10,11)][col],xlab="ICE",ylab="ICE estimated",main=sc)
    abline(a=0,b=1,col="grey")
}


#Plot distribution ICE vs estimated distribution ICE
col<-0
for (sc in c("A","B","C","D","E","F","G","H","I")){
    col<-col+1
    mu.true<- -10 -5
    var.true<-10^2+5^2
    
    x<-seq(mu.true-3*sqrt(var.true),mu.true+3*sqrt(var.true),0.1)
    plot(x,sapply(x,function(y){dnorm(y,mean=mu.true,sd=sqrt(var.true))}),lwd=2,xlab="ICE",ylab="pdf",type="l",main=sc)
    
    z<-NULL
    for (ID in 1:length(Fig3[Fig3$sc==sc,]$mu.ICE)){
        if(ID%%100==0){print(c(sc,ID))}
        z<-rbind(z,unlist(sapply(x,function(y){dnorm(y,mean=Fig3[which(Fig3$sc==sc & Fig3$ID==ID),]$mu.ICE,sd=sqrt(Fig3[which(Fig3$sc==sc & Fig3$ID==ID),]$sigma.ICE))})))
    }
    
    lines(x, colMeans(z),col=brewer.pal(11, "Spectral")[c(1,3,4,5,7,8,9,10,11)][col])
}

#plot distribution of modes

col<-0
for (sc in c("A","B","C","D","E","F","G","H","I")){
    
    col<-col+1
    mu.true<- -10 -5
    var.true<-10^2+5^2
    
    x<-seq(mu.true-3*sqrt(var.true),mu.true+3*sqrt(var.true),0.1)
    plot(x,sapply(x,function(y){dnorm(y,mean=mu.true,sd=sqrt(var.true))}),lwd=2,xlab="ICE",ylab="pdf",type="l",main=sc)
    
    test<-approx(density(Fig3[Fig3$sc==sc,c(5)]), xout = x)
    lines(x,test$y,col=brewer.pal(11, "Spectral")[c(1,3,4,5,7,8,9,10,11)][col])
    
    
    
}





