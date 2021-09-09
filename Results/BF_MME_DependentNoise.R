library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)

################ functions ##################
# create noisy networks
getObs<-function(data.true,Nv,alpha,beta,frac)
{
  data.true <- as.matrix(data.true)
  ind0 <- which(data.true==0,arr.ind = TRUE)
  ind0 <- ind0[ind0[,1]<ind0[,2],] 
  ind1 <- which(data.true==1,arr.ind = TRUE)
  ind1 <- ind1[ind1[,1]<ind1[,2],] 
  num_obs0 <- min(nrow(ind0),max(0,floor(rnorm(1,nrow(ind0)*alpha,sqrt(frac*nrow(ind0)*alpha*(1-alpha))))))
  num_obs1 <- min(nrow(ind1),max(0,floor(rnorm(1,nrow(ind1)*beta,sqrt(frac*nrow(ind1)*beta*(1-beta))))))
  ind0_sample <- ind0[sample(1:nrow(ind0),num_obs0),]
  ind1_sample <- ind1[sample(1:nrow(ind1),nrow(ind1)-num_obs1),]
  data.obs <- sparseMatrix(i=c(ind0_sample[,1],ind1_sample[,1]),j=c(ind0_sample[,2],ind1_sample[,2]),x = 1, symmetric = TRUE,dims = c(Nv,Nv))
  return(data.obs)
}
# get the value of branching factor - 1
getscaleR0 <- function(degree)
{
  mean(degree^2)/mean(degree)-1
}
# get MME
getOne <- function(data.true,alpha,beta,frac,alpha0=0.09, Nb=500, tiny=0.0001, Zcv=1.96)
{
  scaleR0_true <- getscaleR0(colSums( data.true ))
  P <- ncol(data.true)
  Y <- as.matrix(getObs(data.true,P,alpha,beta,frac) )
  Ystar <- as.matrix(getObs(data.true,P,alpha,beta,frac) )
  Ystar2 <-   as.matrix(getObs(data.true,P,alpha,beta,frac) )
  scaleR0_obs <- getscaleR0(colSums( Y ))
  res_est <- MME_BF(alpha0, P, Y, Ystar, Ystar2,Nb,tiny,Zcv)
  list_est<- c(res_est$Alpha,res_est$Beta,res_est$ScaleR0,res_est$cfScaleR0[1]<scaleR0_true&&res_est$cfScaleR0[2]>scaleR0_true,res_est$cfScaleR0[2]-res_est$cfScaleR0[1])
  return(c(alpha,beta,frac,scaleR0_true,scaleR0_obs,list_est))
}

source("../R/MME_BF.R")
################ load data ##################
indir <-  "../Data/RandomNetworks/"
load( paste0(indir,"randomNetwork.RData") )

################ simulation ##################
# number of trials
m <- 500
rep_num <- 120
# set alpha and beta 
Beta <- rep(c(0.15,0.2),each=4)
Alpha <- Beta/10
Frac <- rep(c(0.5,1,1.5,2),2)
num_case <- length(Beta)

# get MME
data.true1 <- sparseMatrix(i=type5[,1],j=type5[,2],x = 1, symmetric = T)
data.true2 <- sparseMatrix(i=type6[,1],j=type6[,2],x = 1, symmetric = T)
data.true3 <- sparseMatrix(i=type7[,1],j=type7[,2],x = 1, symmetric = T)
data.true4 <- sparseMatrix(i=type8[,1],j=type8[,2],x = 1, symmetric = T)
res_MME1 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true1,Alpha[i],Beta[i],Frac[i]))))
res_MME2 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true2,Alpha[i],Beta[i],Frac[i]))))
res_MME3 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true3,Alpha[i],Beta[i],Frac[i]))))
res_MME4 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true4,Alpha[i],Beta[i],Frac[i]))))
getRes <- function(case)
{
  res <- cbind(case[,1:3],case[,4]+1,abs(case[,5]- case[,4]),abs(case[,6]- case[,1]),abs(case[,7]- case[,2]),abs(case[,8]- case[,4])/case[,4],case[,9:10])
  return(colMeans(res))
}
table_all1 <- t(sapply(1:num_case, function(i) getRes(res_MME1[[i]])))
table_all2 <- t(sapply(1:num_case, function(i) getRes(res_MME2[[i]])))
table_all3 <- t(sapply(1:num_case, function(i) getRes(res_MME3[[i]])))
table_all4 <- t(sapply(1:num_case, function(i) getRes(res_MME4[[i]])))

res1 <- table_all1[,c(4,1:3,8:10)]
res2 <- table_all2[,c(4,1:3,8:10)]
res3 <- table_all3[,c(4,1:3,8:10)]
res4 <- table_all4[,c(4,1:3,8:10)]

res <- rbind(res1,res2,res3,res4)
res[,1] <- c(rep(1,num_case),rep(2,num_case),rep(3,num_case),rep(4,num_case))
colnames(res) <- c("kappa","psi","phi","omega","MAE","RF","Length")
res <- data.frame(res) %>% mutate(Data=case_when(kappa==1 ~ "Homogeneous, average degree 50",kappa==2 ~ "Homogeneous, average degree 100",
kappa==3 ~ "Inhomogeneous, average degree 50",kappa==4 ~ "Inhomogeneous, average degree 100")  ) %>%
  select(-kappa)
################ Figure 6 ##################
res.long <- melt(res, id=c("Data","psi","phi","omega"),measure=c("MAE","RF","Length"))%>%
  mutate(psi = case_when(psi==0.015 ~ "psi=0.015, phi=0.15",psi==0.02 ~ "psi=0.02, phi=0.2"))

 
Fig6 <- res.long %>% mutate(omega = as.numeric(as.character(omega)) + as.numeric(factor(Data))/40- .075 ) %>%
  ggplot(aes(x = omega, y = value,color = factor(Data)))+
  geom_point(alpha = 0.7)   +
  facet_grid(variable~psi, scales = "free") + theme_bw(base_size = 14) +
  scale_x_continuous(breaks = c(0.5,1,1.5,2), name=expression(omega))+  ylab(" ")+
  geom_blank(aes(y = 0))  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),panel.background = element_blank(),legend.position="bottom") +
  labs(color  = "Network")+guides(color=guide_legend(nrow=2,byrow=TRUE))



