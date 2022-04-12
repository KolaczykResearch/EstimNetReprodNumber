library(igraph)

# Set a destination for output
outdir <- "../Data/RandomNetworks/"

# type 1: ER networks, average degree = 50, number of nodes = 1e4
set.seed(1)
Nv <- 1e4
p=50/Nv
g <- erdos.renyi.game(Nv, p)
g_edgelist <- get.edgelist(g)
type1 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree1 <- degree(g)

# type 1: configuration model 
g_conf <- sample_degseq(degree1, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type1_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree1_conf <- degree(g_conf)

# type 2: ER networks, average degree = 100, number of nodes = 1e4
set.seed(2)
p=100/Nv
g <- erdos.renyi.game(Nv, p)
g_edgelist <- get.edgelist(g)
type2 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree2 <- degree(g)

# type 2: configuration model 
g_conf <- sample_degseq(degree2, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type2_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree2_conf <- degree(g_conf)

# type 3: scale-free networks, average degree = 50, number of nodes = 1e4
set.seed(3)
g <- barabasi.game(Nv, power = 2, m = 25,directed=FALSE) 
g_edgelist <- get.edgelist(g)
type3 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree3 <- degree(g)

# type 3: configuration model 
g_conf <- sample_degseq(degree3, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type3_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree3_conf <- degree(g_conf)

# type 4: scale-free networks, average degree = 100, number of nodes = 1e4
set.seed(4)
g <- barabasi.game(Nv, power = 2, m = 50,directed=FALSE) 
g_edgelist <- get.edgelist(g)
type4 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree4 <- degree(g)

# type 4: configuration model 
g_conf <- sample_degseq(degree4, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type4_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree4_conf <- degree(g_conf)
 
# degree matrix for 4 networks
degree_matrix <- cbind(degree1,degree2,degree3,degree4)
degree_matrix_conf <- cbind(degree1_conf,degree2_conf,degree3_conf,degree4_conf) 

# type 5: ER networks, average degree = 50, number of nodes = 500
set.seed(5)
Nv <- 5e2
p=50/Nv
g <- erdos.renyi.game(Nv, p)
g_edgelist <- get.edgelist(g)
type5 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree5 <- degree(g)

# type 5: configuration model 
g_conf <- sample_degseq(degree5, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type5_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree5_conf <- degree(g_conf)
 

# type 6: ER networks, average degree = 100, number of nodes = 500
set.seed(6)
p=100/Nv
g <- erdos.renyi.game(Nv, p)
g_edgelist <- get.edgelist(g)
type6 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree6 <- degree(g)

# type 6: configuration model 
g_conf <- sample_degseq(degree6, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type6_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree6_conf <- degree(g_conf)

# type 7: scale-free networks, average degree = 50, number of nodes = 500
set.seed(7)
g <- barabasi.game(Nv, power = 2, m = 25,directed=FALSE) 
g_edgelist <- get.edgelist(g)
type7 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree7 <- degree(g)

# type 7: configuration model 
g_conf <- sample_degseq(degree7, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type7_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree7_conf <- degree(g_conf)

# type 8: scale-free networks, average degree = 100, number of nodes = 500
set.seed(8)
g <- barabasi.game(Nv, power = 2, m = 50,directed=FALSE) 
g_edgelist <- get.edgelist(g)
type8 <- cbind(g_edgelist[,2],g_edgelist[,1])
degree8 <- degree(g)

# type 8: configuration model 
g_conf <- sample_degseq(degree8, method="vl")
g_conf_edgelist <- get.edgelist(g_conf)
type8_conf <- cbind(g_conf_edgelist[,2],g_conf_edgelist[,1])
degree8_conf <- degree(g_conf)

# save data
save(type1,type2,type3,type4,type5,type6,type7,type8,degree_matrix,file = paste0(outdir,"randomNetwork.RData"))

type1 <- type1_conf
type2 <- type2_conf
type3 <- type3_conf
type4 <- type4_conf
type5 <- type5_conf
type6 <- type6_conf
type7 <- type7_conf
type8 <- type8_conf
degree_matrix <- degree_matrix_conf
save(type1,type2,type3,type4,type5,type6,type7,type8,degree_matrix,file = paste0(outdir,"randomNetwork_configuration.RData"))

