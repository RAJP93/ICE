library(mvtnorm)

#Example 1
EX1 <- read.csv("~/PhD/Infering_idividual_causal/SIM2.csv")


Fig<-NULL
Fig1<-NULL

col<-0

for(i in c(1:1000)){
    col<-col+1
    data<-EX1[EX1$ID==i,]
    #Simulation parameters
    beta0<-0
    beta1<- -0.2
    beta2<- -0.1
    betaC<- 4
    
    sigma0<-0.25
    tau0<-0.25
    tau1<-0.5
    tau2<-0.25
    
    ind<-0
    
    for(s in c(3,10,100)){
        print(c(i,s))
        ind<-ind+1
        #<-100
        #s<-length(data[,1])
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
        Y.ak<-Y.k+(a.1-A.1)*(beta1+unique(data$u1))+(a.2-A.2)*(beta2+unique(data$u2))
        Y.0<-Y.k+(0-A.1)*(beta1+unique(data$u1))+(0-A.2)*(beta2+unique(data$u2))
        ICE<-exp(Y.ak)-exp(Y.0)
        
        library(cubature)
        
        
        int.CWC<-function(y){
            if(A.2==a.2){
                sum(na.omit(sapply(seq(mu[2]-3*sqrt(sigma[2,2]),mu[2]+3*sqrt(sigma[2,2]),0.01), function(x){(0.01)*dmvnorm(c(x,-log(max(0,(exp(Y.k+(a.1-A.1)*(beta1+x))-y)/(exp(Y.k-A.1*(beta1+x)))))-beta2),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
            }
            else{
                if(A.1==a.1){
                    sum(na.omit(sapply(seq(mu[3]-3*sqrt(sigma[3,3]),mu[3]+3*sqrt(sigma[3,3]),0.01), function(x){(0.01)*dmvnorm(c(-log(max(0,(exp(Y.k+(a.2-A.2)*(beta2+x))-y)/(exp(Y.k-A.2*(beta2+x)))))-beta1,x),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
                }
                else{
                    # sum(na.omit(sapply(seq(mu[2]-3*sqrt(sigma[2,2]),mu[2]+3*sqrt(sigma[2,2]),0.01), function(x){
                    #     (0.01)*dmvnorm(c(x,log(max((exp(Y.k)+y)/(exp(Y.k+(beta1+x))),0))-beta2),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
                    sum(na.omit(sapply(seq(mu[2]-3*sqrt(sigma[2,2]),mu[2]+3*sqrt(sigma[2,2]),0.01), function(x){
                        (0.01)*dmvnorm(c(x,log(max((exp(Y.k)+y),0))-Y.k-(beta1+x)-beta2),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
                    
                }
            }
        }
        
        a<-1
        b<-1
        x<- 0
        
        dist<-c(round(ICE,digits=1)+x,int.CWC(round(ICE,digits=1)+x))
        
        while(!(a==0 && b==0)){
            x<- x +0.1
            if(a==1){
                a2<-int.CWC(round(ICE,digits=1)+x)
                dist<-rbind(dist,c(round(ICE,digits=1)+x,a2))
                a=(a2>10^-4)
                
            }
            if(b==1){
                b2<-int.CWC(round(ICE,digits=1)-x)
                dist<-rbind(dist,c(round(ICE,digits=1)-x,b2))
                b=(b2>10^-4)
            }
            
        }
        
        dist<-dist[order(dist[,1]),]
        
        
        Fig1<-rbind(Fig1,c(i,ICE,s,dist[which(dist[,2]==max(dist[,2])),]))
        Fig<-rbind(Fig,cbind(rep(i,length(dist[,1])),rep(ICE,length(dist[,1])),rep(s,length(dist[,1])),dist[,1],dist[,2]))
        
        
    }
}

Fig1<-as.data.frame(Fig1)
colnames(Fig1)<-c("ID","ICE","Repeats","ICE.hat","pdf")

Fig<-as.data.frame(Fig)
colnames(Fig)<-c("ID","ICE","Repeats","x","pdf")

#colnames(Fig1)<-c("ID","Ya","Y0","ICE","Repeats","TRTs","mu.ICE","sigma.ICE")

#Fig 8b

par(mfrow=c(1,1),mar=c(5,5,5,5))

a<- -20
b<- 30

#tikz(file = "D:/Documents/PhD/Infering_idividual_causal/TikZFigures/Fig9LN.tex", width = 7, height = 6)

plot((Fig[which(Fig$Repeats==100 & Fig$ID==100),c(4,5)]),lwd=2,xlab="CWC",ylab="pdf",type="l",xlim=c(a,b),col=rgb(231/255,41/255,138/255))
lines((Fig[which(Fig$Repeats==10 & Fig$ID==100),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",lty=2,col=rgb(231/255,41/255,138/255))
lines((Fig[which(Fig$Repeats==3 & Fig$ID==100),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",lty=3,col=rgb(231/255,41/255,138/255))

abline(v=Fig[which(Fig$Repeats==3 & Fig$ID==100),2][1],col=rgb(231/255,41/255,138/255),lwd=2)

lines((Fig[which(Fig$Repeats==100 & Fig$ID==263),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",lty=1,xlim=c(a,b),col=rgb(102/255,166/255,30/255))
lines((Fig[which(Fig$Repeats==10 & Fig$ID==263),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",lty=2,xlim=c(a,b),col=rgb(102/255,166/255,30/255))
lines((Fig[which(Fig$Repeats==3 & Fig$ID==263),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",lty=3,xlim=c(a,b),col=rgb(102/255,166/255,30/255))

abline(v=Fig[which(Fig$Repeats==3 & Fig$ID==263),2][1],xlim=c(a,b),col=rgb(102/255,166/255,30/255),lwd=2)

lines((Fig[which(Fig$Repeats==100 & Fig$ID==974),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",xlim=c(a,b),col=rgb(230/255,171/255,2/255))
lines((Fig[which(Fig$Repeats==10 & Fig$ID==974),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",lty=2,col=rgb(230/255,171/255,2/255))
lines((Fig[which(Fig$Repeats==3 & Fig$ID==974),c(4,5)]),lwd=2,xlab=NA,ylab=NA,type="l",lty=3,col=rgb(230/255,171/255,2/255))

abline(v=Fig[which(Fig$Repeats==3 & Fig$ID==974),2][1],lwd=2,col=rgb(230/255,171/255,2/255))

legend("topright",inset=c(-0.2,0),legend=c(3,10,100),lty=c(3,2,1),,xpd=T,lwd=2,title="Repeats")
#dev.off()


############################################
##########    Inference ####################
############################################

library(mvtnorm)

#Example 1
EX1 <- read.csv("~/PhD/Infering_idividual_causal/SIM2.csv")

random <- read.csv("~/PhD/Infering_idividual_causal/EX2_cov.csv")
fixed <- read.csv("~/PhD/Infering_idividual_causal/EX2_fixed.csv")

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
Fig3<-NULL

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
        
        #s<-length(data[,1])
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
        Y.ak<-Y.k+(a.1-A.1)*(beta1+unique(data$u1))+(a.2-A.2)*(beta2+unique(data$u2))
        Y.0<-Y.k+(0-A.1)*(beta1+unique(data$u1))+(0-A.2)*(beta2+unique(data$u2))
        ICE<-exp(Y.ak)-exp(Y.0)
        
        library(cubature)
        
        
        int.CWC<-function(y){
            if(A.2==a.2){
                sum(na.omit(sapply(seq(mu[2]-5*sqrt(sigma[2,2]),mu[2]+5*sqrt(sigma[2,2]),0.01), function(x){(0.01)*dmvnorm(c(x,-log(max(0,(exp(Y.k+(a.1-A.1)*(beta1+x))-y)/(exp(Y.k-A.1*(beta1+x)))))-beta2),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
            }
            else{
                if(A.1==a.1){
                    sum(na.omit(sapply(seq(mu[3]-5*sqrt(sigma[3,3]),mu[3]+5*sqrt(sigma[3,3]),0.01), function(x){(0.01)*dmvnorm(c(-log(max(0,(exp(Y.k+(a.2-A.2)*(beta2+x))-y)/(exp(Y.k-A.2*(beta2+x)))))-beta1,x),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
                }
                else{
                    # sum(na.omit(sapply(seq(mu[2]-3*sqrt(sigma[2,2]),mu[2]+3*sqrt(sigma[2,2]),0.01), function(x){
                    #     (0.01)*dmvnorm(c(x,log(max((exp(Y.k)+y)/(exp(Y.k+(beta1+x))),0))-beta2),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
                    sum(na.omit(sapply(seq(mu[2]-5*sqrt(sigma[2,2]),mu[2]+5*sqrt(sigma[2,2]),0.01), function(x){
                        (0.01)*dmvnorm(c(x,log(max((exp(Y.k)+y),0))-Y.k-(beta1+x)-beta2),mu[c(2:3)],sigma[c(2:3),c(2:3)])})))
                    
                }
            }
        }
        
        a<-1
        b<-1
        x<- 0
        
        dist<-c(round(ICE,digits=1)+x,int.CWC(round(ICE,digits=1)+x))
        
        while(!(a==0 && b==0)){
            x<- x +0.1
            if(a==1){
                a2<-int.CWC(round(ICE,digits=1)+x)
                dist<-rbind(dist,c(round(ICE,digits=1)+x,a2))
                a=(a2>10^-4)
                
            }
            if(b==1){
                b2<-int.CWC(round(ICE,digits=1)-x)
                dist<-rbind(dist,c(round(ICE,digits=1)-x,b2))
                b=(b2>10^-4)
            }
            
        }
        
        dist<-dist[order(dist[,1]),]
        
        
        Fig2<-rbind(Fig2,c(sc,i,dist[which(dist[,2]==max(dist[,2])),]))
        Fig3<-rbind(Fig3,cbind(rep(sc,length(dist[,1])),rep(i,length(dist[,1])),dist[,1],dist[,2]))
        
        
    }
}

Fig2<-as.data.frame(Fig2)
Fig2[,2]<-as.numeric(Fig2[,2])
Fig2[,3]<-as.numeric(Fig2[,3])
Fig2[,4]<-as.numeric(Fig2[,4])

colnames(Fig2)<-c("sc","ID","ICE.hat","pdf")

True.ICE<-unique(Fig1[,c(1,2)])

Fig2b<-merge(Fig2,True.ICE)




library(RColorBrewer)
par(mfrow=c(3,3),mar=c(2,2,2,2))


#Plot ICE vs estimate E[ICE|H] 
col<-0
for (sc in c("A","B","C","D","E","F","G","H","I")){
    col<-col+1
    plot((Fig2b[Fig2b$sc==sc,c(5,3)]),lwd=2,col=brewer.pal(11, "Spectral")[c(1,3,4,5,7,8,9,10,11)][col],xlab="ICE",ylab="ICE estimated",main=sc)
    abline(a=0,b=1,col="grey")
}

write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNA.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='A'),c(5)],P=Fig2b [which(Fig2b [,2]=='A'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNB.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='B'),c(5)],P=Fig2b [which(Fig2b [,2]=='B'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNC.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='C'),c(5)],P=Fig2b [which(Fig2b [,2]=='C'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLND.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='D'),c(5)],P=Fig2b [which(Fig2b [,2]=='D'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNE.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='E'),c(5)],P=Fig2b [which(Fig2b [,2]=='E'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNF.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='F'),c(5)],P=Fig2b [which(Fig2b [,2]=='F'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNG.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='G'),c(5)],P=Fig2b [which(Fig2b [,2]=='G'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNH.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='H'),c(5)],P=Fig2b [which(Fig2b [,2]=='H'),c(3)]))
write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FigLNI.csv",data.frame(PT=Fig2b [which(Fig2b [,2]=='I'),c(5)],P=Fig2b [which(Fig2b [,2]=='I'),c(3)]))



#Plot distribution ICE vs estimated distribution ICE

Fig3<-as.data.frame(Fig3)
Fig3[,2]<-as.numeric(Fig3[,2])
Fig3[,3]<-as.numeric(Fig3[,3])
Fig3[,4]<-as.numeric(Fig3[,4])

colnames(Fig3)<-c("sc","ID","x","y")

col<-0
for (sc in (c("A","B","C","D","E","F","G","H","I"))){
    col<-col+1
    
    test.x<-seq(min(Fig2b[Fig2b$sc==sc,5])*0.9,max(Fig2b[Fig2b$sc==sc,5])*1.1,0.1)
    
    test0<-approx(density(Fig2b[Fig2b$sc==sc,5]), xout = test.x)
    
    plot(test.x,test0$y,lwd=2,xlab="ICE",ylab="pdf",type="l",main=sc)
    
    test.y<-NULL
    
    ind<-which(design$sc==sc)
    
    n<-design$n[ind]
    s<-design$s[ind]
    
    for(ID in c(1:n)){
    test<-approx(Fig3[which(Fig3$sc==sc & Fig3$ID==ID),c(3,4)], xout = test.x)
    test$y[which(is.na(test$y))]<-0
    test$y<-test$y/sum(test$y*0.1)
    test.y<-rbind(test.y,test$y)

    }
    
    lines(test.x,colMeans(test.y),col=brewer.pal(11, "Spectral")[c(1,3,4,5,7,8,9,10,11)][col])
    
}


#plot distribution of modes

col<-0
for (sc in c("A","B","C","D","E","F","G","H","I")){
    
    col<-col+1
    
    x<-seq(min(Fig2b[Fig2b$sc==sc,5])*0.9,max(Fig2b[Fig2b$sc==sc,5])*1.1,0.1)
    
    test0<-approx(density(Fig2b[Fig2b$sc==sc,5]), xout = x)

    plot(x,test0$y,lwd=2,xlab="ICE",ylab="pdf",type="l",main=sc)
    
    test<-approx(density(Fig2b[Fig2b$sc==sc,3]), xout = x)
    lines(x,test$y,col=brewer.pal(11, "Spectral")[c(1,3,4,5,7,8,9,10,11)][col])
}





