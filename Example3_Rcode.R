library(mvtnorm)

#Example 1
EX1 <- read.csv("~/PhD/Infering_idividual_causal/SIM.csv")

Fig1<-NULL

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
        Y.ak<-Y.k+(a.1-A.1)*(beta1+unique(data$u1))+(a.2-A.2)*(beta2+unique(data$u2))
        Y.0<-Y.k+(0-A.1)*(beta1+unique(data$u1))+(0-A.2)*(beta2+unique(data$u2))
        ICE<-as.numeric(Y.ak>120)-as.numeric(Y.0>120)
        
        mu.Y<-c(Y.k+(a.1-A.1)*(beta1)+(a.2-A.2)*(beta2),Y.k-A.1*beta1-A.2*beta2)+as.matrix(rbind(c(0,(a.1-A.1),(a.2-A.2)),c(0,-A.1,-A.2)))%*%mu
        
        sigma.Y<-as.matrix(rbind(c(0,(a.1-A.1),(a.2-A.2)),c(0,-A.1,-A.2)))%*%sigma%*%t(as.matrix(rbind(c(0,(a.1-A.1),(a.2-A.2)),c(0,-A.1,-A.2))))
        
        if(sigma.Y[1,1]==0){
            
            P00=as.numeric(mu.Y[1]<120)*pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2]))
            P11=as.numeric(mu.Y[1]>120)*(1-pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2])))
            P01=as.numeric(mu.Y[1]<120)*(1-pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2])))
            P10=as.numeric(mu.Y[1]>120)*pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2]))
        }
        else{
            if(sigma.Y[2,2]==0){
                P00=as.numeric(mu.Y[2]<120)*pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1]))
                P11=as.numeric(mu.Y[2]>120)*(1-pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1])))
                P01=as.numeric(mu.Y[2]>120)*pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1]))
                P10=as.numeric(mu.Y[2]<120)*(1-pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1])))    
            }
            else{
                library(mvtnorm)
                P00=pmvnorm(lower=c(min(0,mu.Y[1]-3*sqrt(sigma.Y[1,1])),min(0,mu.Y[2]-3*sqrt(sigma.Y[2,2]))),upper=c(120,120),as.vector(t(mu.Y)),sigma=sigma.Y)
                P11=pmvnorm(lower=c(120,120),upper=c(max(120*1.1,mu.Y[1]+3*sqrt(sigma.Y[1,1])),max(mu.Y[2]+3*sqrt(sigma.Y[2,2]),120*1.1)),as.vector(t(mu.Y)),sigma=sigma.Y)
                P01=pmvnorm(lower=c(min(0,mu.Y[1]-3*sqrt(sigma.Y[1,1])),120),upper=c(120,max(1.1*120,mu.Y[2]+3*sqrt(sigma.Y[2,2]))),as.vector(t(mu.Y)),sigma=sigma.Y)
                P10=pmvnorm(lower=c(120,min(0,mu.Y[2]-3*sqrt(sigma.Y[2,2]))),upper=c(max(1.1*120,mu.Y[1]+3*sqrt(sigma.Y[1,1])),120),as.vector(t(mu.Y)),sigma=sigma.Y)
            }
        }
        
        Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,P00,P11,P01,P10))
        
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"00",P00))
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"11",P11))
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"01",P01))
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"10",P10)) 
        
        
    }
    
    
}
Fig1<-as.data.frame(Fig1)
colnames(Fig1)<-c("ID","Y.ak","Y.0","ICE","Repeats","P00","P11","P01","P10")

mean((Fig1$P00+Fig1$P11)*0+Fig1$P01*(-1)+Fig1$P10*(1))

#Update EX1
mean(((Fig1$P00+Fig1$P11)*0+Fig1$P01*(-1)+Fig1$P10*(1))[EX1[which(EX1$PER==3&EX1$C_lag_1==-0.3),]$ID])

mean(Fig1[EX1[which(EX1$PER==3&EX1$C_lag_1==0.7),]$ID,]$Y.ak>120)-mean(Fig1[EX1[which(EX1$PER==3&EX1$C_lag_1==0.7),]$ID,]$Y.0>120)
mean(Fig1[EX1[which(EX1$PER==3&EX1$C_lag_1==-0.3),]$ID,]$Y.ak>120)-mean(Fig1[EX1[which(EX1$PER==3&EX1$C_lag_1==-0.3),]$ID,]$Y.0>120)


#Fig 9
par(mfrow=c(1,1))
for(i in c(23,290,250)){
    
    
    P <- as.matrix(Fig1[c((1+(i-1)*3):(i*3)),c(6:9)])
    colnames(P) <- c("00","11","01","10")
    rownames(P) <- c(3,10,100)
    P <- prop.table(as.table(P),1)
    print(Fig1[(1+(i-1)*3),c(1:3)])
   # tikz(paste("D:/Documents/PhD/Infering_idividual_causal/TikzFigures/Dist3",i,"b.tex",sep=""),width=3.5,height=3.5)
        barplot(t(P),beside=TRUE,legend=colnames(P),args.legend=list(x="top",title="title",cex=0.75,horiz=T),ylim=c(0,1),xlab="Repeats")
   # dev.off()
        }

FIG10a<-Fig1[Fig1$Repeats=="3",c(1,8,4)]
FIG10b<-Fig1[Fig1$Repeats=="10",c(1,8,4)]
FIG10c<-Fig1[Fig1$Repeats=="100",c(1,8,4)]

colnames(FIG10a)<-c("x","y","label")
colnames(FIG10b)<-c("x","y","label")
colnames(FIG10c)<-c("x","y","label")

#write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FIG10a.csv",FIG10a,row.names=F)
#write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FIG10b.csv",FIG10b,row.names=F)
#write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FIG10c.csv",FIG10c,row.names=F)

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
    print(sc)
    
    ind<-which(design$sc==sc)
    
    n<-design$n[ind]
    s<-design$s[ind]
    
    for(i in c(1:n)){
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
            
        mu.Y<-c(Y.k+(a.1-A.1)*(beta1)+(a.2-A.2)*(beta2),Y.k-A.1*beta1-A.2*beta2)+as.matrix(rbind(c(0,(a.1-A.1),(a.2-A.2)),c(0,-A.1,-A.2)))%*%mu
        
        sigma.Y<-as.matrix(rbind(c(0,(a.1-A.1),(a.2-A.2)),c(0,-A.1,-A.2)))%*%sigma%*%t(as.matrix(rbind(c(0,(a.1-A.1),(a.2-A.2)),c(0,-A.1,-A.2))))
        
        if(sigma.Y[1,1]==0){
            
            P00=as.numeric(mu.Y[1]<120)*pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2]))
            P11=as.numeric(mu.Y[1]>120)*(1-pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2])))
            P01=as.numeric(mu.Y[1]<120)*(1-pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2])))
            P10=as.numeric(mu.Y[1]>120)*pnorm(120,mu.Y[2],sqrt(sigma.Y[2,2]))
        }
        else{
            if(sigma.Y[2,2]==0){
                P00=as.numeric(mu.Y[2]<120)*pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1]))
                P11=as.numeric(mu.Y[2]>120)*(1-pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1])))
                P01=as.numeric(mu.Y[2]>120)*pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1]))
                P10=as.numeric(mu.Y[2]<120)*(1-pnorm(120,mu.Y[1],sqrt(sigma.Y[1,1])))    
            }
            else{
                library(mvtnorm)
                P00=pmvnorm(lower=c(min(0,mu.Y[1]-3*sqrt(sigma.Y[1,1])),min(0,mu.Y[2]-3*sqrt(sigma.Y[2,2]))),upper=c(120,120),as.vector(t(mu.Y)),sigma=sigma.Y)
                P11=pmvnorm(lower=c(120,120),upper=c(max(120*1.1,mu.Y[1]+3*sqrt(sigma.Y[1,1])),max(mu.Y[2]+3*sqrt(sigma.Y[2,2]),120*1.1)),as.vector(t(mu.Y)),sigma=sigma.Y)
                P01=pmvnorm(lower=c(min(0,mu.Y[1]-3*sqrt(sigma.Y[1,1])),120),upper=c(120,max(1.1*120,mu.Y[2]+3*sqrt(sigma.Y[2,2]))),as.vector(t(mu.Y)),sigma=sigma.Y)
                P10=pmvnorm(lower=c(120,min(0,mu.Y[2]-3*sqrt(sigma.Y[2,2]))),upper=c(max(1.1*120,mu.Y[1]+3*sqrt(sigma.Y[1,1])),120),as.vector(t(mu.Y)),sigma=sigma.Y)
            }
        }
        
        Fig2<-rbind(Fig2,c(sc,i,as.numeric(P10),as.numeric(P00+P11),as.numeric(P01)))
        
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"00",P00))
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"11",P11))
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"01",P01))
        #Fig1<-rbind(Fig1,c(i,Y.ak,Y.0,ICE,s,"10",P10)) 
        
        
    }
    
    
}

Fig2<-as.data.frame(Fig2)
Fig2[,c(2)]<-as.numeric(Fig2[,c(2)])
Fig2[,c(3)]<-as.numeric(Fig2[,c(3)])
Fig2[,c(4)]<-as.numeric(Fig2[,c(4)])
Fig2[,c(5)]<-as.numeric(Fig2[,c(5)])

colnames(Fig2)<-c("sc","ID","1","0","-1")

Fig2$MPN<-apply(Fig2[,c(3:5)],1,function(x){c(1,0,-1)[which(x==max(x))]})

True.ICE<-unique(Fig1[,c(1,4)])

Fig3<-merge(Fig2,True.ICE)


#Classification tables
classification<-list()
i<-0
for (sc in c("A","B","C","D","E","F","G","H","I")){
    i<-i+1
    subset<-Fig3[Fig3$sc==sc,]
   classification[[i]]<-prop.table(table(subset$ICE,subset$MPN))
}

#ID vs P(-1) plot
for (sc in c("A","B","C","D","E","F","G","H","I")){
write.csv(file=paste("~/PhD/Infering_idividual_causal/TikZFigures/FigHyp",sc,".csv",sep=""),Fig3[Fig3$sc==sc,c(1,5,7)],row.names=F)
    #write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FIG10b.csv",FIG10b,row.names=F)
    #write.csv(file="~/PhD/Infering_idividual_causal/TikZFigures/FIG10c.csv",FIG10c,row.names=F)
    
}

#Fig 9
par(mfrow=c(3,3),mar=c(2,2,2,2))

for (sc in c("A","B","C","D","E","F","G","H","I")){
    P<-colMeans(Fig3[Fig3$sc==sc,c(3:5)])
    print(P)
    barplot(P,col=c(1,2,3),ylim=c(0,1),main=sc)
    }



