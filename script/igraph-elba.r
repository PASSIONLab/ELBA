#!/usr/bin/env Rscript

# install.packages("igraph", repos = "http://cran.us.r-project.org")
# install.packages("Matrix", repos = "http://cran.us.r-project.org")

library(igraph)
library(Matrix)

if (length(args) == 0)
{
  	stop("One argument must be supplied root vertex", call.=FALSE)
} 	else if (length(args) == 1) 
{
	args <- commandArgs(trailingOnly = TRUE)
  	root <- strtoi(args[1])
}

adjm <- readMM("ecoli.hifi.result.mtx")
G <- graph_from_adjacency_matrix(adjm, weighted=TRUE)

# E(G)$weight
# str(adjm)

# print("Graph loaded")

# G2 <- graph_from_adjacency_matrix(adjm, weighted=TRUE, mode="max") # "max": an undirected graph will be created and max(A(i,j), A(j,i)) gives the edge weights.
# E(G2)$weight

## Callback function to stop seach past 10 dist
f <- function(graph, data, extra) {
  data['dist'] == 10
}

# bfs(G, root, unreachable = FALSE, restricted = NULL, order = TRUE,
#   rank = FALSE, father = FALSE, pred = FALSE, succ = FALSE,
#   dist = TRUE, callback = NULL, extra = NULL,
#   rho = parent.frame())


tmp <- bfs(G, root, unreachable = FALSE, restricted = NULL, order = TRUE,
  rank = FALSE, father = FALSE, pred = FALSE, succ = FALSE,
  dist = TRUE, callback = f, extra = NULL,
  rho = parent.frame())

# str(tmp)
# tmp

vid <- 1
vidvect <- vector()

for (hoop in tmp[['dist']])
{
	if(!is.nan(hoop))
	{
    	# print(hoop)
    	# print(vid)
    	vidvect <- append(vidvect, vid)
	}

	vid = vid + 1
}

# print(vidvect)

## Extract subgraph
S <- induced_subgraph(G, vidvect, impl = c("create_from_scratch")) # "create_from_scratch" constructs the result graph from scratch and then copies the attributes accordingly
E(S)$weight

write_graph(S, "testr1.dot", format = c("dot"))

