NTS5_BF=function(alpha0, P, Y, Ystar, Ystar2, Nb=500, tiny=0.0001, Zcv=1.96){

# This program calculates point-estimates for two error rates Alpha, Beta, 
# edge density, 2-star edge count, triangle count, and clustering coefficient.
# It also calculates confidence intervals for edge density delta, 2-star
# edge count, triangle count, and clustering coefficient based on asymptotic
# normality

        # P: number of nodes
	# Y, Ystar, Ystar2: 3 observed PxP adjacent matrices
        # Nb: No. bootstrap replications
	# Zcv: critical value for confidence interval

P012=P*(P-1)*(P-2); PP=P*P; P01=P*(P-1); P01h=P01/2; rP01h=sqrt(P01h)
H=matrix(1:6, nrow=2)
G=matrix(1:6, nrow=2)
DD=matrix(1:4, nrow=2)
S2=matrix(1:9, ncol=3)
Sv=matrix(nrow=Nb, ncol=2)

# Part 3: Point estimate for Alpha, Beta and Delta
#     NO Part 1, Part 2!!!
u1=sum(Y)/P01
u2=sum(abs(Y-Ystar))/(2*P01)
u3=length(Y[((Ystar2-2*Ystar+Y)==1)|((Ystar2-2*Ystar+Y)==-2)])/(3*P01)
 alpha=alpha0; alpha0=alpha+10*tiny
while(abs(alpha-alpha0)>tiny) { alpha0=alpha
   beta=(u2-alpha0+alpha0*u1)/(u1-alpha0)
   delta=((u1-alpha0)^2)/(u1-u2-2*u1*alpha0+alpha0^2)
   alpha=(u3-delta*(1-beta)*beta^2)/((1-delta)*(1-alpha0)^2)
}
eS=rep(0,24)
eS[1]=max(0,alpha); eS[2]=max(0,beta); eS[3]=delta

# Part 4: confidence interval for Delta
k3=1-eS[1]-eS[2]  # if k3 <=0, DISASTER!!!
tmp1=max(eS[1], eS[2]); tmp2=min(eS[1], eS[2])
if(((k3<=0)|(k3==1))|(tmp1*(1-tmp2)>0.25)) {
    #cat("No decent estimation, exit!")
    #quit()
    slist=list(Alpha=eS[1], Beta=eS[2], Delta=eS[3], Nedge=eS[6], N2s=eS[7], ScaleR0=eS[8],Mdegree=eS[15],Mdegree2=eS[18],
    cfDelta=c(eS[4], eS[5]),
    cfNedge=c(eS[9], eS[10]),
    cfN2s=c(eS[11], eS[12]),
    cfScaleR0=c(eS[13], eS[14]),
    cfMdegree=c(eS[16], eS[17]),
    cfMdegree2=c(eS[19], eS[20]),
    cfAlpha=c(eS[21], eS[22]),
    cfBeta=c(eS[23], eS[24])
    )
    return(slist)
}
k1=eS[1]*(1-eS[1]); k2=eS[2]*(1-eS[2]); k4=eS[2]-eS[1]
# S2 is 3x3 matrix Sigma defined in Appendix
S2[1,1]=eS[3]*k2+(1-eS[3])*k1
S2[1,2]=eS[3]*k2*(eS[2]-.5)+(1-eS[3])*k1*(.5-eS[1]); S2[2,1]=S2[1,2]
S2[1,3]=eS[3]*k2*(eS[2]^2-2*k2)/3+(1-eS[3])*k1*((1-eS[1])^2-2*k1)/3
S2[3,1]=S2[1,3]
S2[2,2]=eS[3]*k2*(.5-k2)+(1-eS[3])*k1*(.5-k1);
S2[2,3]=eS[3]*eS[2]*k2*(1/3-k2)+(1-eS[3])*(1-eS[1])*k1*(1/3-k1)
S2[3,2]=S2[2,3]
S2[3,3]=eS[3]*eS[2]*k2*(1/3-eS[2]*k2)+(1-eS[3])*k1*(1-eS[1])*(1/3-k1*(1-eS[1]))
U3=c(3*k3+6*eS[1]*eS[2]-2, 3*k3+6*eS[2]-2, -2)/(k3^3) # 3rd row of matrix W in Appendix
sigmaE=t(U3)%*%S2%*%U3; sigmaE=sqrt(sigmaE/P01h)
eS[4]=eS[3]-Zcv*sigmaE # lower bound for Delta
eS[5]=eS[3]+Zcv*sigmaE # upper bound for Delta

U1=c((1-2*eS[2])*eS[1]+eS[2]^2, eS[1]-2*eS[2], 1)/(k3^2*(1-eS[3]) ) # 1st row of matrix W in Appendix
sigmaE=t(U1)%*%S2%*%U1; sigmaE=sqrt(sigmaE/P01h)
eS[21]=eS[1]-Zcv*sigmaE # lower bound for alpha
eS[22]=eS[1]+Zcv*sigmaE # upper bound for alpha
U2=c(-(1-2*eS[1])*eS[2]-eS[1]^2, eS[2]-2*eS[1]+1, -1)/(k3^2*eS[3]) # 2nd row of matrix W in Appendix
sigmaE=t(U2)%*%S2%*%U2; sigmaE=sqrt(sigmaE/P01h)
eS[23]=eS[2]-Zcv*sigmaE # lower bound for beta
eS[24]=eS[2]+Zcv*sigmaE # upper bound for beta

# Part 5: Point estimate Nos. of 2-stars, kappa
Y0=Y-eS[1]; diag(Y0)=0; YY=Y0%*%Y0
hatNedge=sum(Y0)/(2*k3)      
diag(YY)=0; hatN2s=sum(YY)/(2*k3*k3) 
eS[7]=hatN2s
eS[6]=hatNedge  
#eS[6]=eS[3]*P01/2 
# density of 2* edges = 2*eS[6]/P012
# density of triangles = 6*eS[7]/P012
eS[8]=eS[7]/eS[6] # point estimate for kappa - 1
eS[15]=2*eS[6]/P # point estimate for average degree
eS[18]=2*(eS[7]+eS[6])/P # point estimate for average degree square


# Part 6: Compute 2x2 asymptotic covariance matrix of eS[6]=hatNedge and eS[7]=hatN2s
                                       # in the form V1n + V2n + V3n
# Part 6.1: compute V2n, V3n
meanY0=sum(Y0)/P01
tt=c(6*k4, 3*(k4^2-k1-k2),2*(k4*(-6*eS[1]*eS[2]+3*k3^2-4*k3)+(1-eS[1])*(eS[2]-2*eS[1])))
H[2,]=tt*4*eS[7]/(3*P012)+c(6*k1, 3*k1*(1-2*eS[1]), 2*k1*(1-eS[1])*(1-3*eS[1]))*2/(3*k3^2)*meanY0
H[1,]=tt*eS[6]*2/P01/3+c(6*k1, 3*k1*(1-2*eS[1]), 2*k1*(1-eS[1])*(1-3*eS[1]))/(3*k3)
G[1,]=c((1-2*eS[2])*eS[1]+eS[2]^2, eS[1]-2*eS[2], 1)/((1-eS[3])*k3^2)
G[2,]=c(-(1-2*eS[1])*eS[2]-eS[1]^2, eS[2]-2*eS[1]+1, -1)/(eS[3]*k3^2)
DD[2,]=c(eS[7]*2/P012-meanY0/k3, eS[7]*2/P012)*2/k3
DD[1,]=c(eS[6]*2/P01-1, eS[6]*2/P01)/k3
tt=DD%*%G
V3n=tt%*%t(H); V3n=(V3n+t(V3n))*0.5
V2n=tt%*%S2%*%t(tt)


# Part 6.2: bootstrap computing V1n
# Find eta1 eta2
if(abs(eS[1]-eS[2])<tiny) { eta2=eS[1]; eta1=1-2*eta2}
if((eS[2]-eS[1])>tiny) {  tmp1=sqrt(1-4*eS[1]*(1-eS[2])); tmp2=sqrt(1-4*eS[2]*(1-eS[1]))  
	eta2=0.5*(1-tmp1)
	if((tmp1+tmp2)<0.5) eta1=0.5*(tmp1+tmp2) 
	else eta1=0.5*(tmp1-tmp2)
}
if((eS[1]-eS[2])>tiny) { tmp1=sqrt(1-4*eS[1]*(1-eS[2])); tmp2=sqrt(1-4*eS[2]*(1-eS[1])) 
	eta2=0.5*(1+tmp1); eta1=0.5*(-tmp1+tmp2) 
}

for(nb in 1:Nb) {
E1=matrix(rep(0, PP), nrow=P)
E1[upper.tri(E1, diag=F)]=sample(c(-1,0,1), P01h, replace=T, prob=c(1-eta1-eta2, eta1, eta2))
E2=t(E1); E1[lower.tri(E1, diag=F)]=E2[lower.tri(E2, diag=F)] # make E1 symmetric
Yb=Y*ifelse(E1==0, 1, 0)+ifelse(E1==1, 1, 0)
Yb=Yb-Y*eta1-eta2; diag(Yb)=0
YbY0=Yb%*%Y0; Y0Yb=Y0%*%Yb; 
Sv[nb,1]=sum(Yb)*(rP01h/(P01*k3))
diag(YbY0)=0; diag(Y0Yb)=0
Sv[nb,2]=sum(YbY0+Y0Yb)*(rP01h/(P012*k3*k3))
}  # bootstrap replications end here
Vn=(var(Sv)+V2n+V3n)/P01h
tmp1=sqrt(Vn[1,1])*P01/2; tmp2=sqrt(Vn[2,2])*P012/2
eS[9]=eS[6]-Zcv*tmp1; eS[10]=eS[6]+Zcv*tmp1 # confidence bounds for No. edge
eS[16]=eS[15]-Zcv*tmp1*2/P; eS[17]=eS[15]+Zcv*tmp1*2/P # confidence bounds for average degree
eS[11]=eS[7]-Zcv*tmp2; eS[12]=eS[7]+Zcv*tmp2 # confidence bounds for No. 2-star edge
tmp=sqrt(eS[7]^2*tmp1^2/eS[6]^4-eS[7]/eS[6]^3*Vn[1,2]*P012*P01/2+tmp2^2/eS[6]^2)  
# var of estimated kappa - 1
eS[13]=eS[8]-Zcv*tmp; eS[14]=eS[8]+Zcv*tmp  # confidence bounds for kappa - 1

tmp=2/P*sqrt(tmp1^2+tmp2^2+Vn[1,2]*P012*P01/2 )
# var of estimated  average degree square
eS[19]=eS[18]-Zcv*tmp; eS[20]=eS[18]+Zcv*tmp  # confidence bounds for average degree square

# Summary statistics from the outputs:
slist=list(Alpha=eS[1], Beta=eS[2], Delta=eS[3], Nedge=eS[6], N2s=eS[7], ScaleR0=eS[8],Mdegree=eS[15],Mdegree2=eS[18],
cfDelta=c(eS[4], eS[5]),
cfNedge=c(eS[9], eS[10]),
cfN2s=c(eS[11], eS[12]),
cfScaleR0=c(eS[13], eS[14]),
cfMdegree=c(eS[16], eS[17]),
cfMdegree2=c(eS[19], eS[20]),
cfAlpha=c(eS[21], eS[22]),
cfBeta=c(eS[23], eS[24])
)
return(slist)
} # The end 
