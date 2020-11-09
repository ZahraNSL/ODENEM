#' construct a network structure
#'
#' @param S_num number of Hidden nodes
#' @param P_num number of perturbation nodes
#' @param E_num number of Observations nodes per each hidden nodes
#' @param theta matrix of prior information for S_E connections
#' @param g_s hidden layer network structure
#'
#' @return an igraph of all nodes and connection as true network
#' @export
sample_network <- function (S_num,
                            P_num,
                            E_num,
                            theta,
                            g_s) {
  # set nodes' names
  P_names = paste("P", 1:P_num, sep = "")
  E_names = paste("E", 1:(E_num * S_num), sep = "")
  rownames(theta) = E_names

  # make an empty graph from S_nodes, E_nodes and P_nodes
  newgraph = g_s + P_names + E_names
  #graph(g,n=S_num+P_num+(E_num*S_num),directed = T)

  #set color for nodes
  V(newgraph)$color <- c(paste(rep("yellow", S_num), sep = ""),
                         paste(rep("red", P_num), sep = ""),
                         paste(rep("white", (E_num * S_num)), sep = ""))

  #connect P layer and S layer
  edgelist = (rbind(P_names, S_names))
  newgraph = igraph::add_edges(newgraph, edgelist)

  igraph::plot.igraph(newgraph)

  #connect S layer and E layer regarding theta matrix
  edgelist =   unlist(sapply(1:(nrow(theta)), function(j) {
    rbind(S_names[which(theta[j,] != 0)], rownames(theta)[j])
  }))

  newgraph = igraph::add_edges(newgraph, edgelist)

  return(newgraph)
}
