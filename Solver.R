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
library(gurobi)

#############################################
###                                       ###
###   READ A SOLUTION FILE AND RETURN XY  ###
###                                       ###
#############################################

readSolutionFromFile <- function(filename, n, gurobi){
  
  fileConn <- file(filename)
  lines <- readLines(fileConn)
  if(gurobi){
    ones <- unlist(sapply(2:length(lines), function(i){
      varname <- strsplit(lines[i], " ")[[1]][1]
      value <- strsplit(lines[i], " ")[[1]][2]
      if(value == 1 && substr(varname,1,1)=="x")
        as.numeric(substr(varname, 2, nchar(varname)))
    }))
    xy <- array(F, 2*n)
    xy[ones] <- T
  } else{
    ones <- sapply(3:length(lines), function(i) {
      varname <- strsplit(lines[i], " ")[[1]][1]
      if(substr(varname,1,1)=="x")
        as.numeric(substr(varname, 2, nchar(varname)))
    })
    ones <- unlist(ones)
    xy <- array(F, 2*n)
    xy[ones] <- T
  }
  close(fileConn)
  
  return(xy)
}


#############################################
###                                       ###
###  CREATE A LINEAR PROGRAM FOR THE VSP  ###
###                                       ###
#############################################

createLinearProgram <- function(g, lambda, weights, forced_in_A, forced_in_B){
  
  n <- length(V(g))
  E <- length(E(g))
  
  # Maximize |A| + lambda * |B| subject to |A| <= |B| and that there is no edge between A and B
  # The variable is the concatenation of x and y. (size 2*n)
  
  # Set the objective function to |A| + lambda * |B|
  objective_function <- c(weights, lambda*weights)
  
  # Initialize the constraints
  number_of_constraints <- 2*E+n+4
  constraint_matrix <- matrix(0,number_of_constraints,2*n)
  direction <- array(0, number_of_constraints)
  right_hand_side <- array(0, number_of_constraints)
  
  # For each edge (u,v) add the constraints  (i) x_u + y_v <= 1  and (ii) x_v + y_u <= 1
  k <- 1
  for (i in 1:E){
    edge <- as.numeric(ends(g, i))
    
    constraint_matrix[k, c(edge[1], edge[2]+n)] <- 1
    direction[k] <- "<="
    right_hand_side[k] <- 1
    k <- k+1
    
    constraint_matrix[k, c(edge[1]+n, edge[2])] <- 1
    direction[k] <- "<="
    right_hand_side[k] <- 1
    k <- k+1
  }
  
  # For each vertex u add the constraint (iii) x_u + y_u <= 1
  for (i in 1:n){
    constraint_matrix[k, c(i, i+n)] <- 1
    direction[k] <- "<="
    right_hand_side[k] <- 1
    k <- k+1
  }
  
  # Specify that there must be at least one vertex in each pile
  constraint_matrix[k,] <- c(array(1, n), array(0, n))
  direction[k] <- ">="
  right_hand_side[k] <- 1
  k <- k+1
  
  constraint_matrix[k,] <- c(array(0, n), array(1, n))
  direction[k] <- ">="
  right_hand_side[k] <- 1
  k <- k+1
  
  # Specifiy that |A| <= |B|
  constraint_matrix[k,] <- c(weights, -weights)
  direction[k] <- "<="
  right_hand_side[k] <- 0
  k <- k+1
  
  # Set the lower bound for the separator size
  constraint_matrix[k,] <- array(1, 2*n)
  direction[k] <- "<="
  right_hand_side[k] <- n-graph.cohesion(g)
  
  
  # Force the vertices in forced_in_A to be in A and the vertices in forced_in_B to be in B
  number_of_force_constraints <- length(forced_in_A) + length(forced_in_B) 
  if(number_of_force_constraints > 0){
    #    force_constraints <- matrix(0, number_of_force_constraints, 2*n)
    #    p <- 1
    #    for(node_in_A in forced_in_A){
    #      force_constraints[p, node_in_A] <- 1
    #      p <- p + 1
    #    }
    #    for(node_in_B in forced_in_B){
    #      force_constraints[p, node_in_B + n] <- 1
    #      p <- p + 1
    #    }
    number_of_force_constraints <- 2*number_of_force_constraints -  1
    
    # The first inequality set that a1 is either in A or B (not in separator)
    # The two after set that b1 is in the opposite subset of a1
    force_constraints <- matrix(0, number_of_force_constraints, 2*n)
    a1 <- forced_in_A[1]
    b1 <- forced_in_B[1]    
    force_constraints[1,c(a1, a1+n)] <- 1
    force_constraints[2, c(a1, b1+n) ] <- c(1, -1)
    force_constraints[3, c(a1+n, b1)] <- c(1, -1)
    p <- 4
    
    # For all others elements, set the a_i in the same subset as a_1. And also set the b_i in the same subset as b1
    if(length(forced_in_A)>=2){
      for(i in 2:length(forced_in_A)){
        force_constraints[p ,c(a1, forced_in_A[i])] <- c(1,-1)
        p <- p+1
        force_constraints[p ,c(a1+n, forced_in_A[i]+n)] <- c(1,-1)
        p <- p+1
      }
    }
    if(length(forced_in_B)>=2){
      for(i in 2:length(forced_in_B)){
        force_constraints[p ,c(b1, forced_in_B[i])] <- c(1,-1)
        p <- p+1
        force_constraints[p ,c(b1+n, forced_in_B[i]+n)] <- c(1,-1)
        p <- p+1
      }
    }   
    
    constraint_matrix <- rbind(constraint_matrix, force_constraints)
    direction <- c(direction, array("=", number_of_force_constraints))
    right_hand_side <- c(right_hand_side, c(1, array(0, number_of_force_constraints-1)))
    
  }
  
  
  return(list("objective_function"=objective_function, "constraint_matrix"=constraint_matrix, "direction"=direction, "right_hand_side"=right_hand_side))
  
}

#############################################
###                                       ###
###  SOLVE VSP IN R USING GUROBI LIBRARY  ###
###                                       ###
#############################################


solve_vsp_R <- function(g, lambda, weights=array(1, length(V(g))), forced_in_A, forced_in_B, verbose=F) {
  # Create Linear Program
  linearProgram <- createLinearProgram(g, lambda, weights, forced_in_A, forced_in_B)  
  
  # Create the gurobi model
  model <- list()
  model$A <- linearProgram$constraint_matrix
  model$obj <- linearProgram$objective_function
  model$modelsense <- "max"
  model$rhs <- linearProgram$right_hand_side
  model$sense <- linearProgram$direction
  model$vtype <- 'B'
  params <- list(OutputFlag=0+verbose) 
  
  # Solve model
  
  ptm <- proc.time()
  res <- gurobi(model, params)
  ptm <- proc.time()-ptm
  
  elapsed_time <- round(ptm[3], digits=4)
  
  # print(paste("Solved in", elapsed_time, "seconds."))
  # cat(paste(elapsed_time, ",", sep=""))
  
  return(res$x)
}


#############################################
###                                       ###
###  WRITE AN LP FILE FROM A LINEAR PROG  ###
###                                       ###
#############################################

# If distanceMatrix is not defined, the linear program is max c(A) + lambda*c(B)
# If distanceMatrix is defined, the linear program is max [(1-mu) (c(A) + lamda*c(B)) + mu * sum(z_ij d(i,j)) ]
# z_ij = 1 if i and j are in the same subset, otherwise z_ij = 0
write_lp_file_from_LP <- function(filename, linearProgram, distanceMatrix=NULL, mu=NULL){
  objective_function <- linearProgram$objective_function
  constraint_matrix <- linearProgram$constraint_matrix
  direction <- linearProgram$direction
  right_hand_side <- linearProgram$right_hand_side
  
  n <- length(objective_function)/2
  number_of_constraints <- nrow(constraint_matrix)
  
  variable_names <- paste("x", 1:(2*n), sep="")
  
  if(is.null(distanceMatrix)){
    lines <- array(NA, 1+2+number_of_constraints+1+2*n+1)
  } else{
    lines <- array(NA, 1+2+number_of_constraints+1+2*n+1+ 3*n*(n-1)/2)
  }
  
  
  lines[1] <- "Maximize"
  
  if(is.null(distanceMatrix)){
    objective_line <- paste("\tobj:", paste(paste(objective_function, variable_names), collapse= " + "))
    
  } else{
    objective_line <- paste("\tobj:", paste(paste( (1-mu) * objective_function, variable_names), collapse= " + "))
    objective_distance <- paste(sapply(1:(n-1), function(i) {
      paste(sapply((i+1):n, function(j) {
        paste(-mu*distanceMatrix[j,i], " z", i*n+j, sep="")
      }), collapse = " ")
    }), collapse= " ")
    objective_line <- paste(objective_line, " " ,objective_distance, sep="")
    
  }
  
  lines[2:3] <- c(objective_line, "Subject to")
  
  
  lines[3+(1:number_of_constraints)] <- sapply(1:number_of_constraints, function(i) {
    matline <- ifelse(constraint_matrix[i,]>0, paste("+", constraint_matrix[i,]), constraint_matrix[i,])
    currentLine <- paste("\tc", i, ": ", paste(paste(matline, variable_names)[matline != "0"], collapse=" "), " ", direction[i], " ", right_hand_side[i], sep="")
  })
  if(!is.null(distanceMatrix)){
    
    lines[3+number_of_constraints + 1:(n*(n-1)/2)] <- paste("\tc", number_of_constraints + 1:(n*(n-1)/2), ": ", unlist(sapply(1:(n-1), function(i) {
      sapply((i+1):n, function(j) {
        paste("z", i*n+j, " - x",i, " - x",j , " >= ", -1, sep="")
      })
    })), sep="")
    
    number_of_constraints <- number_of_constraints + n*(n-1)/2
    lines[3+number_of_constraints + 1:(n*(n-1)/2)] <- paste("\tc", number_of_constraints + 1:(n*(n-1)/2), ": ", unlist(sapply(1:(n-1), function(i) {
      sapply((i+1):n, function(j) {
        paste("z", i*n+j, " - x",i+n, " - x",j+n , " >= ", -1, sep="")
      })
    })), sep="")
    
    number_of_constraints <- number_of_constraints + n*(n-1)/2
    
  }
  
  lines[(number_of_constraints+4)] <- "Binary"
  
  lines[(number_of_constraints+4) + 1:(2*n)] <- paste("\t", variable_names, sep="")
  if(!is.null(distanceMatrix)){
    lines[(number_of_constraints+4) + 2*n + 1:(n*(n-1)/2)] <- paste("\t", unlist(sapply(1:(n-1), function(i) {
      sapply((i+1):n, function(j) { paste("z", i*n+j , sep="") })
    })), sep="")
  }
  
  lines[length(lines)] <- "End"
  
  fileConn<-file(filename)
  writeLines(lines, fileConn)
  close(fileConn)
  
}




#############################################
###                                       ###
###   GENERAL FUNCTION TO SOLVE THE VSP   ###
###                                       ###
#############################################

solve_vsp <- function(g, lambda, weights=array(1, length(V(g))), mu=NULL, distanceMatrix=NULL, forced_in_A=c(), forced_in_B=c(), gurobi=T, verbose=F, command_line=F){
  if(is.null(distanceMatrix) && !command_line){
    xy <- solve_vsp_R(g, lambda, weights, forced_in_A, forced_in_B, verbose)
  } else {
    lpFile <- paste(temporaryDirectory, "model.lp", sep="")
    solFile <- paste(temporaryDirectory, "model.sol", sep="")
    system(paste("rm", lpFile))
    system(paste("rm", solFile))
    
    linearProgram <- createLinearProgram(g, lambda, weights, forced_in_A, forced_in_B)
    n <- length(V(g))
    write_lp_file_from_LP(lpFile, linearProgram, distanceMatrix, mu)
    if(gurobi){
      system(paste("gurobi_cl ResultFile=", solFile, " ", lpFile, sep=""))
    } else{
      system(paste(fscipExecutable, scipParamsFile, lpFile, "-fsol", solFile, sep=" "))
    }
    xy <- readSolutionFromFile(solFile, n, gurobi)
  }
  return(xy > 0)
}


#############################################
###                                       ###
###     SOLVE ITERATIVE VERSION OF VSP    ###
###                                       ###
#############################################


solve_vsp_iterative <- function(g, lambda, weights=array(1, length(V(g))), MAX_ITER, MIN_SIZE, verbose=F){
  current_iter <- 0
  membership <- array("", length(V(g)))
  V(g)$label <- as.numeric(V(g))
  return(solve_vsp_recursive(g, lambda, weights, current_iter, membership, MAX_ITER, MIN_SIZE, verbose))
}

solve_vsp_recursive <- function(g, lambda, weights, current_iteration, membership, MAX_ITER, MIN_SIZE, verbose){
  # If we can not separate the graph anymore, just return membership
  if(current_iteration >= MAX_ITER || length(V(g)) < MIN_SIZE){
    return(membership)
  } else {
    n <- length(V(g))
    V(g)$name <- 1:n
    # Solve the VSP
    xy <- solve_vsp(g, lambda, weights, verbose=verbose)
    A <- c(array(TRUE, n), array(FALSE, n))
    pilA <- xy[A]
    pilB <- xy[!A]
    sep <- !pilA & !pilB
    
    # Get all the nodes in both piles and in separator
    # The indexes given in V(g)$label are the indexes of the nodes in the original graph so that we can modify membership
    nodes_in_A <- V(g)$label[pilA]
    nodes_in_B <- V(g)$label[pilB]
    nodes_in_sep <- V(g)$label[sep]
    
    # Update membership
    membership[nodes_in_A] <- paste(membership[nodes_in_A], "A", sep="")
    membership[nodes_in_B] <- paste(membership[nodes_in_B], "B", sep="")
    membership[nodes_in_sep] <- paste(membership[nodes_in_sep], "S", sep="") 
    
    components_in_A <- 0
    components_in_B <- 0
    
    # Remove the separator from the graph to analyze each connected components
    graph_without_separator <- delete.vertices(g, sep)
    clusters <- clusters(graph_without_separator)
    
    # For each cluster
    #   Create the corresponding graph
    #   Add the cluster's number in membership to differentiate the all the connected components
    #   Call recursively solve_vsp_recursive
    for(k in 1:clusters$no){
      vertices_subset <- clusters$membership==k
      connected_component <- induced.subgraph(graph_without_separator, vertices_subset)
      
      nodes_in_connected_component <- V(connected_component)$label
      
      if(min(nodes_in_connected_component) %in% nodes_in_A){
        components_in_A <- components_in_A + 1 
        membership[nodes_in_connected_component] <- paste(membership[nodes_in_connected_component], components_in_A, "-", sep="")
      } else {
        components_in_B <- components_in_B + 1 
        membership[nodes_in_connected_component] <- paste(membership[nodes_in_connected_component], components_in_B, "-", sep="")
      }
      
      membership <- solve_vsp_recursive(connected_component, lambda, weights[V(connected_component)$label], current_iteration+1, membership, MAX_ITER, MIN_SIZE, verbose)
    }
    
    return(membership)
  }
  
}

#############################################
###                                       ###
###     SOLVE VSP WITH BETA CONSTRAINT    ###
###                                       ###
#############################################

solve_vsp_beta <- function(g, beta, weights=array(1, length(V(g))), forced_in_A=NULL, forced_in_B=NULL, verbose=F) {
  lambda <- 1
  linearProgram <- createLinearProgram(g, lambda, weights, forced_in_A, forced_in_B)  
  n <- length(linearProgram$objective_function)/2
  
  beta_constraints <- as.numeric(c(array(0, n), array(1, n)))
  beta_constraints <- rbind(beta_constraints, 1-beta_constraints)
  linearProgram$constraint_matrix <- rbind(linearProgram$constraint_matrix, beta_constraints)
  linearProgram$direction <- c(linearProgram$direction, "<=", "<=")
  linearProgram$right_hand_side <- c(linearProgram$right_hand_side, beta, beta)
  
  model <- list()
  model$A <- linearProgram$constraint_matrix
  model$obj <- linearProgram$objective_function
  model$modelsense <- "max"
  model$rhs <- linearProgram$right_hand_side
  model$sense <- linearProgram$direction
  model$vtype <- 'B'
  params <- list(OutputFlag=0+verbose) 
  
  ptm <- proc.time()
  res <- gurobi(model, params)
  ptm <- proc.time()-ptm
  
  elapsed_time <- round(ptm[3], digits=4)
  
  # print(paste("Solved in", elapsed_time, "seconds."))
  # cat(paste(elapsed_time, ",", sep=""))
  
  return(res$x)
}

vsp_file <- function(filename, g, lambda, weights=array(1, length(V(g))), mu=NULL, distanceMatrix = NULL, forced_in_A=NULL, forced_in_B=NULL){
  linearProgram <- createLinearProgram(g, lambda, weights, forced_in_A, forced_in_B)
  write_lp_file_from_LP(filename, linearProgram, distanceMatrix, mu)
}





