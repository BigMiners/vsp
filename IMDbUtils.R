###############################################################################
### VSP implementation for Knowledge Discovery in Graphs Through Vertex     ###
### Separation                                                              ###
###                                                                         ###  
### Copyright (C) 2017  Marc Sarfati, Marc Queudot,                         ###
###                     Catherine Mancel, Marie-Jean Meurs                  ###
###                                                                         ###
### Permission is hereby granted, free of charge, to any person obtaining a ###
### copy of this software and associated documentation files                ### 
### (the "Software"), to deal in the Software without restriction,          ###
### including without limitation the rights to use, copy, modify, merge,    ###
### publish, distribute, sublicense, and/or sell copies of the Software,    ###
### and to permit persons to whom the Software is furnished to do so,       ###
### subject to the following conditions:                                    ###
###                                                                         ###
### The above copyright notice and this permission notice shall be included ###
### in all copies or substantial portions of the Software.                  ###
###                                                                         ###
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS ###
### OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF              ###
### MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  ###
### IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY    ###
### CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,    ###
### TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE       ###
### SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                  ###
###############################################################################

getActorName <- function(artistid){
  return(as.character(artists$fullname[artists$artistid==artistid]))
}

getMovieName <- function(movieid){
  return(as.character(movies$imdbtitle[movies$movieid==movieid]))
}

getActorsMovies <- function(artistid){
  return(roles$movieid[roles$artistid==artistid])
}

printActorsMovies <- function(artistid){
  print(sapply(getActorsMovies(artistid), getMovieName))
}

getNodeActorNames <- function(node, actorsInNode=actorsInNode) {
  actors <- actorsInNode[[node]]
  return (sapply(actors, getActorName))
}

getNodeMovies<- function(node, actorsInNode=actorsInNode) {
  return(getActorsMovies(actorsInNode[[node]][1]))
}

printNodeMovies<- function(node, actorsInNode=actorsInNode) {
  return(sapply(getActorsMovies(actorsInNode[[node]][1]), getMovieName))
}

reduceInputSize <- function(nMovies = 20){
  return(rolesData[rolesData$movieid <= sort(unique(rolesData$movieid))[nMovies], ])
}

keepBiggestConnectedComponent <- function(g){
  cl <- clusters(g)
  group <- which(cl$csize == max(cl$csize))
  return(induced_subgraph(g, cl$membership == group))
  
}

readSolutionFromScip <- function(filename, n){
  fileConn <- file(filename)
  lines <- readLines(fileConn)
  
  ones <- sapply(3:length(lines), function(i) {
    varname <- strsplit(lines[i], " ")[[1]][1]
    as.numeric(substr(varname, 2, nchar(varname)))
  })
  
  close(fileConn)
  
  xy <- array(F, 2*n)
  xy[ones] <- T
  return(xy)
}

regroupActors <- function(g){
  visited <- array(F, length(V(g)))
  cluster <- 1
  mapping <- array(0, length(V(g)))
  
  
  while(prod(visited) == 0){
    currentNode <- min(which(visited==F))
    visited[currentNode] <- T
    mapping[currentNode] <- cluster
    ## compare list of movies
    nodeMovies <- getNodeMovies(currentNode, actorsInNode)
    for(n in neighbors(g, currentNode)){
      neighborMovies <- getNodeMovies(n, actorsInNode)
      if(setequal(nodeMovies, neighborMovies)) {
        visited[n] <- T
        mapping[n] <- cluster
      }
    }
    cluster <- cluster +1
  }
  return(mapping)
}
