library(igraph)

# Set a destination for output
outdir <- "../../Data/RandomNetworks/"

# type 1: ER networks, average degree = 50
set.seed(1)
Nv <- 1e4
p=50/Nv
g <- erdos.renyi.game(Nv, p)
g_edgelist <- get.edgelist(g)
type1 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree1 <- degree(g)

# type 2: ER networks, average degree = 100
set.seed(2)
p=100/Nv
g <- erdos.renyi.game(Nv, p)
g_edgelist <- get.edgelist(g)
type2 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree2 <- degree(g)

# type 3: scale-free networks, average degree = 50
set.seed(3)
g <- barabasi.game(Nv, power = 2, m = 25,directed=FALSE) 
g_edgelist <- get.edgelist(g)
type3 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree3 <- degree(g)

# type 4: scale-free networks, average degree = 100
set.seed(4)
g <- barabasi.game(Nv, power = 2, m = 50,directed=FALSE) 
g_edgelist <- get.edgelist(g)
type4 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree4 <- degree(g)
 
# degree matrix for 4 networks
degree_matrix <- cbind(degree1,degree2,degree3,degree4)
 
# save data

save(type1,type2,type3,type4,degree_matrix,file = paste0(outdir,"randomNetwork.RData"))


