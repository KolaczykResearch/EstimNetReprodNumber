library(igraph)
library(Matrix)
library(dplyr)
library(ggplot2)

################ functions ##################
# get the values of beta and alpha for heterogeneous noise
getAlphaBeta<-function(data.true,Nv,alpha_seq,beta_seq)
{
  degree <- colSums( data.true )
  cutoff <- quantile(degree,probs = c(.2,.8))
  degree_low_index <- which(degree <= cutoff[1])
  degree_high_index <- which(degree >= cutoff[2])
  degree_mid_index <- setdiff(1:Nv,c(degree_low_index,degree_high_index))
  mid_index <- c(degree_low_index,degree_mid_index)
  num_01 <- alpha_seq[2]*(Nv*(Nv-1)-sum(degree))
  num_10 <- beta_seq[2]*sum(degree)
  num_10_low <- beta_seq[1]*sum(data.true[degree_low_index,degree_low_index])
  num_01_low <- alpha_seq[1]*(length(degree_low_index)*(length(degree_low_index)-1) -sum(data.true[degree_low_index,degree_low_index]))
  num_10_mid <- beta_seq[2]*(sum(data.true[mid_index,mid_index])-sum(data.true[degree_low_index,degree_low_index]))
  num_01_mid <- alpha_seq[2]*( length(mid_index)*(length(mid_index)-1)-length(degree_low_index)*(length(degree_low_index)-1)-(sum(data.true[mid_index,mid_index])-sum(data.true[degree_low_index,degree_low_index])) )
  beta_high <- (num_10-num_10_low-num_10_mid)/(sum(degree)-sum(data.true[mid_index,mid_index]))
  alpha_high <- (num_01-num_01_low-num_01_mid)/(Nv*(Nv-1)-length(mid_index)*(length(mid_index)-1)-(sum(degree)-sum(data.true[mid_index,mid_index])) )
  return(list(alpha_seq=c(alpha_seq,alpha_high),beta_seq=c(beta_seq,beta_high)))
}
# get observed networks
getObs<-function(data.true,Nv,alpha_seq,beta_seq)
{
  degree <- colSums( data.true )
  cutoff <- quantile(degree,probs = c(.2,.8))
  degree_low_index <- which(degree <= cutoff[1])
  degree_high_index <- which(degree >= cutoff[2])
  degree_mid_index <- setdiff(1:Nv,c(degree_low_index,degree_high_index))
  data.obs <- rsparsematrix(Nv, Nv, nnz = rbinom(1,Nv^2,alpha_seq[3]), rand.x = function(n) return(1))
  mid_index <- c(degree_low_index,degree_mid_index)
  data.obs[mid_index,mid_index] <- rsparsematrix(length(mid_index),length(mid_index), nnz = rbinom(1,length(mid_index)^2,alpha_seq[2]),rand.x = function(n) return(1))
  data.obs[degree_low_index,degree_low_index] <- rsparsematrix(length(degree_low_index),length(degree_low_index), nnz = rbinom(1,length(degree_low_index)^2,alpha_seq[1]), rand.x = function(n) return(1))
  index_edge <- which(data.true==1,arr.ind = TRUE)
  data.obs[index_edge] <- rbinom(nrow(index_edge),1,1-beta_seq[3])
  index_edge_mid <- index_edge[which(index_edge[,1]%in%mid_index | index_edge[,2]%in%mid_index),]
  data.obs[index_edge_mid] <- rbinom(nrow(index_edge_mid),1,1-beta_seq[2])
  index_edge_low <- index_edge[which(index_edge[,1]%in%degree_low_index | index_edge[,2]%in%degree_low_index),]
  data.obs[index_edge_low] <- rbinom(nrow(index_edge_low),1,1-beta_seq[1])
  diag(data.obs) <- 0
  data.obs <- Matrix::forceSymmetric(data.obs)
  return(data.obs)
}
# get the value of branching factor - 1
getscaleR0 <- function(degree)
{
  mean(degree^2)/mean(degree)-1
}
# get MME
getOne <- function(data.true,alpha,beta,alpha0=0.09, Nb=500, tiny=0.0001, Zcv=1.96)
{
  scaleR0_true <- getscaleR0(colSums( data.true ))
  P <- ncol(data.true)
  AlphaBeta <- getAlphaBeta(data.true,P,alpha,beta)
  alpha_seq <- AlphaBeta$alpha_seq
  beta_seq <- AlphaBeta$beta_seq
  Y <- as.matrix(getObs(data.true,P,alpha_seq,beta_seq) )
  Ystar <- as.matrix(getObs(data.true,P,alpha_seq,beta_seq) )
  Ystar2 <-   as.matrix(getObs(data.true,P,alpha_seq,beta_seq) )
  scaleR0_obs <- getscaleR0(colSums( Y ))
  res_est <- MME_BF(alpha0, P, Y, Ystar, Ystar2,Nb,tiny,Zcv)
  list_est<- c(res_est$Alpha,res_est$Beta,res_est$ScaleR0,res_est$cfScaleR0[1]<scaleR0_true&&res_est$cfScaleR0[2]>scaleR0_true,res_est$cfScaleR0[2]-res_est$cfScaleR0[1])
  return(c(alpha_seq,beta_seq,scaleR0_true,scaleR0_obs,list_est))
}


source("../R/MME_BF.R")

################ load data ##################
indir <-  "../Data/RandomNetworks/"
load( paste0(indir,"randomNetwork.RData") )

################ simulation ##################
# number of trials
m <- 500
 
# set alpha and beta
beta1 <- rep(0.1,4)
Beta <- rbind(beta1*c(1,.9,.8,.7),beta1)
Alpha <- Beta/10
num_case <- length(beta1)

# get MME
data.true1 <- sparseMatrix(i=type5[,1],j=type5[,2],x = 1, symmetric = T)
data.true2 <- sparseMatrix(i=type6[,1],j=type6[,2],x = 1, symmetric = T)
data.true3 <- sparseMatrix(i=type7[,1],j=type7[,2],x = 1, symmetric = T)
data.true4 <- sparseMatrix(i=type8[,1],j=type8[,2],x = 1, symmetric = T)
res_MME1 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true1,Alpha[,i],Beta[,i]))))
res_MME2 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true2,Alpha[,i],Beta[,i]))))
res_MME3 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true3,Alpha[,i],Beta[,i]))))
res_MME4 <- lapply(1:num_case,function (i) t(replicate(m ,getOne(data.true4,Alpha[,i],Beta[,i]))))
getRes <- function(case)
{
  case <- case[,-c(3,2,6,5)]
  res <- cbind(case[,1:2],case[,3]+1,abs(case[,4]- case[,3]),abs(case[,5]- case[,1]),abs(case[,6]- case[,2]),abs(case[,7]- case[,3])/case[,3],case[,8:9])
  return(colMeans(res))
}
table_all1 <- t(sapply(1:num_case, function(i) getRes(res_MME1[[i]])))
table_all2 <- t(sapply(1:num_case, function(i) getRes(res_MME2[[i]])))
table_all3 <- t(sapply(1:num_case, function(i) getRes(res_MME3[[i]])))
table_all4 <- t(sapply(1:num_case, function(i) getRes(res_MME4[[i]])))

res1 <- table_all1[,c(3,1:2,7:9)]
res2 <- table_all2[,c(3,1:2,7:9)]
res3 <- table_all3[,c(3,1:2,7:9)]
res4 <- table_all4[,c(3,1:2,7:9)]

res <- rbind(res1,res2,res3,res4)
res[,1] <- c(rep(1,num_case),rep(2,num_case),rep(3,num_case),rep(4,num_case))
colnames(res) <- c("kappa","alpha","beta","MAE","RF","Length")
res <- data.frame(res) %>% mutate(Data=case_when(kappa==1 ~ "Homogeneous, average degree 50",kappa==2 ~ "Homogeneous, average degree 100",
                                                        kappa==3 ~ "Inhomogeneous, average degree 50",kappa==4 ~ "Inhomogeneous, average degree 100") ) %>%
  select(-kappa)  %>% mutate(Case=case_when(alpha==0.01 ~ "case 1",alpha==.009 ~ "case 2",alpha==.008 ~ "case 3",alpha==.007 ~ "case 4") ) %>%
  mutate(beta=0.1-beta)

################ Figure 5 ##################
res.long <- melt(res, id=c("Data","alpha","beta","Case"),measure=c("MAE","RF","Length"))


Fig5 <- res.long %>% mutate(beta = as.numeric(as.character(beta)) + as.numeric(factor(Data))/4000- .00075 ) %>%
  ggplot(aes(x =beta, y = value,color = factor(Data)))+
  geom_point(alpha = 0.7)   +
  facet_grid(variable~. , scales = "free") + theme_bw(base_size = 14) +
  scale_x_continuous(name="Edge noise settings",breaks = c(0,0.01,0.02,0.03),labels = c("case 1","case 2","case 3","case 4"))+  ylab(" ")+
  geom_blank(aes(y = 0))  +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),panel.background = element_blank()) +
  labs(color  = "Network")




 





