library(sand)
library(dplyr)
library(networkDynamic)
library(igraph)
data(hc)

# set a destination for output
outdir <- "../../Data/HospitalNetworks/"

vids <- sort(unique(c(hc$ID1, hc$ID2)))
g.week <- graph.data.frame(hc[, c("ID1", "ID2","Time")], vertices=data.frame(vids),directed=FALSE)
E(g.week)$Time <- E(g.week)$Time / (60 * 60)
status <- unique(rbind(data.frame(id=hc$ID1,status=hc$S1), data.frame(id=hc$ID2, status=hc$S2)))
V(g.week)$Status <-as.character(status[order(status[,1]),2])
E(g.week)$weight <- 1

# create 3 networks (Tues, Wed, Thur)
star_hour <- c(11,35,59)
g.sl3<- lapply(1:3, function(i) {
  g <- subgraph.edges(g.week,E(g.week)[Time > star_hour[i] &Time <= (star_hour[i]+24)])
  edge_data <- data.frame(ID1=get.edgelist(g)[,1],ID2=get.edgelist(g)[,2],Time=E(g)$Time) %>% 
    group_by(ID1,ID2) %>% summarise(counts=n()) %>% subset(counts>5*3)
  g_sparse <- graph_from_edgelist(as.matrix(edge_data[,-3]), directed = FALSE)
   igraph::simplify(g_sparse)
})


## create true network 
library(plyr)
edge_all <- sapply(g.sl3,get.edgelist)
df <- data.frame(a=c(edge_all[[1]][,1],edge_all[[2]][,1],edge_all[[3]][,1]),b=c(edge_all[[1]][,2],edge_all[[2]][,2],edge_all[[3]][,2]) )
edge_rep <- ddply(df,.(a,b),nrow)
edge_true <- as.matrix( edge_rep[edge_rep[,3]>1,-3])
g_true <- graph_from_edgelist(edge_true, directed = FALSE) #ecount(g_true) 337 vcount(g_true) 50
save(g_true,file =paste0(outdir,"dataTrue.RData"))

## create network replicates
node_com <- Reduce(intersect, sapply(1:3,function(i) V(g.sl3[[i]])$name))
adj_list <-  lapply(1:3, function(i) {
g_sub <- induced_subgraph(g.sl3[[i]],node_com)
get.adjacency(g_sub)
})
save(adj_list,file = paste0(outdir,"Replicates.RData"))



 
