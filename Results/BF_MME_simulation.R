library(igraph)
library(Matrix)
library(ggplot2)
library(dplyr)

################ functions ##################
# create noisy networks
getObs<-function(data.true,Nv,alpha,beta)
{
  data.obs <- rsparsematrix(Nv, Nv, nnz = rbinom(1,Nv^2,alpha), rand.x = function(n) return(1))
  ind1 <- which(data.true==1)
  data.obs[ind1] <- rbinom(length(ind1),1,1-beta)
  diag(data.obs)<- 0
  data.obs <-  Matrix::forceSymmetric(data.obs,uplo="L")
  return(data.obs)
}
# get the value of branching factor - 1
getscaleR0 <- function(degree)
{
  mean(degree^2)/mean(degree)-1
}
 
# get MME
getOne <- function(data.network,alpha,beta,alpha0=0.09, Nb=500, tiny=0.0001, Zcv=1.96)
{
  data.true <- get.adjacency(data.network)
  scaleR0_true <- getscaleR0(colSums( data.true ))
  P <- ncol(data.true)
  Y <- as.matrix(getObs(data.true,P,alpha,beta) )
  Ystar <- as.matrix(getObs(data.true,P,alpha,beta) )
  Ystar2 <-   as.matrix(getObs(data.true,P,alpha,beta) )
  scaleR0_obs <- getscaleR0(colSums( Y ))
  res_est <- MME_BF(alpha0, P, Y, Ystar, Ystar2,Nb,tiny,Zcv)
  list_est<- c(res_est$Alpha,res_est$Beta,res_est$ScaleR0,res_est$cfScaleR0[1]<scaleR0_true&&res_est$cfScaleR0[2]>scaleR0_true,res_est$cfScaleR0[2]-res_est$cfScaleR0[1])
  return(c(alpha,beta,scaleR0_true,scaleR0_obs,list_est))
}

source("../../R/MME_BF.R")
################ load data ##################
indir <-  "../../Data/HospitalNetworks/"
load( paste0(indir,"dataTrue.RData") )
indir <-  "../../Data/SchoolNetworks/"
load( paste0(indir,"dataTrue.RData") )

################ simulation ##################
# number of trials
m <- 500

# set alpha and beta 
beta <- rep(c(0.1,0.15,0.2),2)
alpha <- rep(c(0.005,0.01),each=3)

# hospital
res_MME <- lapply(1:length(beta),function (i) t(replicate(m ,getOne(g_true,alpha[i],beta[i]))) ) 
getRes <- function(case)
{
  res <- cbind(case[,1:2],case[,3]+1,abs(case[,4]- case[,3]),abs(case[,5]- case[,1]),abs(case[,6]- case[,2]),abs(case[,7]- case[,3])/case[,3],case[,8:9])
  return(colMeans(res))
}
table_all <- t(sapply(1:length(beta), function(i) getRes(res_MME[[i]])))
res <- table_all[,c(3,1:2,7:9)]
res[,1] <- rep(1,6)
colnames(res) <- c("kappa","alpha","beta","MAE","RF","Length")
res_hosp <- data.frame(res) %>% mutate(Data="Hospital") %>% select(-kappa)  %>%  mutate(contact=1)

# school
res_MME1 <- lapply(1:length(beta),function (i) t(replicate(m ,getOne(school1_undiCol_true,alpha[i],beta[i]))) ) 
res_MME2 <- lapply(1:length(beta),function (i) t(replicate(m ,getOne(school2_undiCol_true,alpha[i],beta[i]))) ) 
res_MME3 <- lapply(1:length(beta),function (i) t(replicate(m ,getOne(school3_undiCol_true,alpha[i],beta[i]))) ) 
res_MME4 <- lapply(1:length(beta),function (i) t(replicate(m ,getOne(school4_undiCol_true,alpha[i],beta[i]))) ) 

table_all1 <- t(sapply(1:length(beta), function(i) getRes(res_MME1[[i]])))
table_all2 <- t(sapply(1:length(beta), function(i) getRes(res_MME2[[i]])))
table_all3 <- t(sapply(1:length(beta), function(i) getRes(res_MME3[[i]])))
table_all4 <- t(sapply(1:length(beta), function(i) getRes(res_MME4[[i]])))

res1 <- table_all1[,c(3,1:2,7:9)]
res2 <- table_all2[,c(3,1:2,7:9)]
res3 <- table_all3[,c(3,1:2,7:9)]
res4 <- table_all4[,c(3,1:2,7:9)]

res <- rbind(res1,res2,res3,res4)
res[,1] <- c(rep(1,6),rep(2,6),rep(3,6),rep(4,6))
colnames(res) <- c("kappa","alpha","beta","MAE","RF","Length")
res_school <- data.frame(res) %>% mutate(Data=case_when(kappa==1 ~ "School 1",kappa==2 ~ "School 2",
                                                        kappa==3 ~ "School 3",kappa==4 ~ "School 4") ) %>%
  select(-kappa)%>%  mutate(contact=2)


################ Figure 5.1 ##################
res <- rbind(res_hosp,res_school)
res.long <- melt(res, id=c("Data","alpha","beta","contact"),measure=c("MAE","RF","Length"))%>% 
  mutate(alpha = case_when(alpha==0.005 ~ "alpha: 0.005",alpha==0.01 ~ "alpha: 0.01" ))
 
Fig5.1 <-res.long %>% mutate(beta = as.numeric(as.character(beta)) + as.numeric(factor(Data))/400- .0075 ) %>%
  ggplot(aes(x = beta, y = value,color = factor(Data),shape=factor(contact)))+
  geom_point(alpha = 0.7)   +
  facet_grid(variable~alpha, scales = "free") + theme_bw() +
  scale_x_continuous(breaks = c(0.1, .15, .2), name=expression(beta))+  ylab(" ")+
  geom_blank(aes(y = 0))  +
  theme(plot.title = element_text(hjust = 0.5,size=16),axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),panel.background = element_blank()) +
  guides(colour = guide_legend(override.aes = list(shape = c(16,rep(17,4) ) )))+scale_shape(guide = FALSE)+
  labs(color  = "Network")



