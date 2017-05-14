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

source("config.R")
source("Solver.R")
source("IMDbUtils.R")


# Load the data if they are not loaded yet
if(!exists("movies")){

  library(gdata)
  library(igraph)
  yearregex <- "\\(([0-9]{4})"
  path <- IMDbdatasetLocation
  
  # Read the 3 tables : roles, artists, movies
  
  roles <- data.frame(read.table(paste(path, "Roles.csv", sep=""), sep="\t"))
  names(roles) <- c("artistid", "movieid", "isDirector")
  artists <- data.frame(read.table(paste(path, "Artists.csv", sep=""), sep="\t"))
  names(artists) <- c("artistid", "fullname")
  movies <- data.frame(read.table(paste(path, "AllMovies.csv", sep=""), sep="\t"))
  names(movies) <- c("movieid", "votes", "nvotes", "rating", "genre", "country", "imdbtitle")
  
  # Parse the movies year of production and define a new column
  movies$year <- 0
  movies$year[grepl(yearregex, movies$imdbtitle)] = as.numeric(substr(regmatches(movies$imdbtitle, regexpr(yearregex, movies$imdbtitle)), 2,5))
  
  # Add ratings and number of movies of each artist in the artist table, for statistics purposes
  # Compute a score for the artists equal to the sum of the ratings of the movies he/she contributed to
  rolesmodified <- merge(roles, movies[, c("movieid", "rating")])
  ratings <- aggregate(rating ~ artistid, rolesmodified, mean)
  artists <- merge(artists, ratings, by="artistid")
  
  counts <- aggregate(rating ~ artistid, rolesmodified, length)
  artists <- merge(artists, counts, by="artistid")
  names(artists) <- c("artistid", "fullname", "rating", "numMovies")
  
  artists$score <- artists$rating*artists$numMovies 
    
  rolesData <- roles
}


# Reduce the size of the graph
roles <- reduceInputSize(nMovies = nMovies)

# Actors in node is a list indexed by the vertices. actorsInNode[[u]] contains the ids of all the actors in u
actorsInNode <- list()
E <- rbind()

# Create the graph. For each film add a link between all the actors
for(movieid in unique(roles$movieid)){
  actors <- roles$artistid[roles$movieid==movieid]
  if(length(actors)>1){
    E <- rbind(E, t(combn(actors, 2)))
  } else { E <- rbind(E, c(actors[1], actors[1]))}
}

# Simplify the graph to remove double edges or loops
g <- simplify(as.undirected(graph.data.frame(E)))
totalActors <- length(V(g))
print(paste("Number of actors :", totalActors))
print(paste("Number of movies :", length(unique(unlist(sapply(as.numeric(V(g)$name), getActorsMovies))))))

# Keep the biggest connected component
g <- keepBiggestConnectedComponent(g)

totalActorsConnectedComponent <- length(V(g))

print(paste("Number of actors in connected component :", totalActorsConnectedComponent))

dict <- as.numeric(V(g)$name)
V(g)$name <- 1:length(V(g))

# Fill actorInNode
for(i in 1:length(dict)){ 
  actorsInNode[[length(actorsInNode)+1]] <- c(dict[i])
}

# Regroup all the actors who have the exact same set of movies
mapping <- regroupActors(g)
g <- simplify(contract.vertices(g, mapping,
                                vertex.attr.comb=list(weight="sum", "ignore")))

# Rename the vertices  1 to n.
# Update actorsInNode
V(g)$name <- 1:length(V(g))
updateMap <- list()
first <- T
for(i in 1:max(mapping)){ 
  first <- T
  for(k in which(mapping==i)){
    if(first){
      first <- F
      updateMap[[i]] <- unique(actorsInNode[[k]])
    } else {
      updateMap[[i]] <- unique(c(updateMap[[i]], actorsInNode[[k]]))
    }
  }
}

actorsInNode <- updateMap

# Set the weights of the nodes to the number of actors it contains
freq <- as.data.frame(table(mapping))$Freq
V(g)$size <- 3+10*log(1+5*freq/max(freq))

n <- length(V(g))

# Solve the Iterative version of the VSP
membership <- solve_vsp_iterative(g, lambda, freq, MAX_ITER, MIN_SIZE, verbose=verbose)

#fileConn <- file("~/Desktop/IMDb1700.graph.solution")
#membership <- readLines(fileConn)
#close(fileConn)

# Only keep the relevant components (not too small)

colors <- c("lightblue", "yellow", "orange", "pink", "green", "darkolivegreen1", "dodgerblue", "gold", "deeppink", "hotpink", "green1", "yellow1", "blue1")
uniques <- unique(membership)
firstsep <- c("A", "B", "S")
V(g)$color <- sapply(membership, function(group){
  return(colors[1 + which(uniques==group) %% length(colors)])
})

V(g)$color[grepl("S", membership)] <-"red"
# colors <- c("green", "yellow", "red")
# V(g)$color <- unlist(sapply(substr(membership, 1, 1), function(letter) colors[which(firstsep==letter)]))
set.seed(0)

png(paste(temporaryDirectory, "separated_graph.png", sep=""))
plot(g, vertex.label=NA)
dev.off()

# Remove separators from groups
uniques <- uniques[-grep("S",uniques)]
# Remove small connected sets
too_small <- sapply(uniques, function(group) {sum(membership==group) < MIN_SIZE})
uniques <- uniques[!too_small]
uniques <- sort(uniques)



nodes_in_groups <- lapply(uniques, function(group) {
  nodes_in_group <- which(membership==group)
})

# Print statistics about the movies in these groups
movies_in_groups <- list()
for(i in 1:length(nodes_in_groups)){
   movies_in_current_group <- unique(unlist(lapply(nodes_in_groups[[i]], function(node) { getNodeMovies(node, actorsInNode) })))
   movies_in_groups[[i]] <- movies[movies$movieid %in% movies_in_current_group,]
   print(summary(movies_in_groups[[i]][c("year", "country", "genre")]))
}

uniques

allCountries <- unique(movies$country)

## Plot statistics
years <- lapply(1:length(uniques), function(i) movies_in_groups[[i]]$year)
countries <- lapply(1:length(uniques), function(i) sapply(movies_in_groups[[i]]$country, function(country) which(country==allCountries)))

png(paste(temporaryDirectory, "years.png", sep=""))
plot(years[[1]], ylim=c(1900, 2020))
dev.off()

#for(i in 2:length(years)) points(years[[i]], col=colors[i])

#i <- 0
#i <- i+1
#plot(countries[[i]], years[[i]], ylim=c(1900, 2020), xlim=c(0, max(unlist(countries))), col=colors[i], pch=i)
#for(i in 2:length(years)) points(countries[[i]], years[[i]], col=colors[i], pch=i)

