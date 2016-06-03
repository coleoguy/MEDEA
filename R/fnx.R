selection <- function(pop, N, s) {
  # create an empty vector for absolute fitness
  abs.fit <- vector(mode = "numeric", length = 16)
  abs.fit[1:16] <- 1 - s * c(4, 3, 3, 2, 3, 2, 2, 1, 3, 2, 2, 1, 2, 1, 1, 0)
  # convert absolute fitness to relative or marginal fitness
  w <- abs.fit / max(abs.fit[pop != 0]) * (pop / N)
  # use marginal fitness to draw from population
  as.vector(rmultinom(1:16, size = N, prob = w))
}

# this fnx takes a population table and generates the gamete pool
makeGametes <- function(pop, r, N) {
  # sperm pool x[1:4]: AB, Ab, aB, ab
  # egg pool x[5:13] AB,  Ab, aB, ab, Ab', ab', a'B, a'b, a'b'
  # in these cases a prime symbol indicates "poisoned" for that
  # allele.
  #
  #
  #create an empty vector for gametes
  x <- vector(length = 13, mode = "numeric")
  # sperm
  # AB
  x[1] <- .5 * sum(
      pop[1],
      pop[2] * (1 / 2),
      pop[3] * (1 / 2),
      pop[4] * (1 / 2) * (1 - r),
      pop[5] * (1 / 2),
      pop[7] * (1 / 2) * r,
      pop[9] * (1 / 2),
      pop[10] * (1 / 2) * r,
      pop[13] * (1 / 2) * (1 - r)
    )
  # Ab
  x[2] <- .5 * sum(
      pop[2] * (1 / 2),
      pop[4] * (1 / 2) * r,
      pop[5] * (1 / 2),
      pop[6],
      pop[7] * (1 / 2) * (1 - r),
      pop[8] * (1 / 2),
      pop[10] * (1 / 2) * (1 - r),
      pop[13] * (1 / 2) * r,
      pop[14] * (1 / 2)
    )
  # aB
  x[3] <- .5 * sum(
      pop[3] * (1 / 2),
      pop[4] * (1 / 2) * r,
      pop[7] * (1 / 2) * (1 - r),
      pop[9] * (1 / 2),
      pop[10] * (1 / 2) * (1 - r),
      pop[11],
      pop[12] * (1 / 2),
      pop[13] * (1 / 2) * r,
      pop[15] * (1 / 2)
    )
  # ab
  x[4] <- .5 * sum(
      pop[4] * (1 / 2) * (1 - r),
      pop[7] * (1 / 2) * r,
      pop[8] * (1 / 2),
      pop[10] * (1 / 2) * r,
      pop[12] * (1 / 2),
      pop[13] * (1 / 2) * (1 - r),
      pop[14] * (1 / 2),
      pop[15] * (1 / 2),
      pop[16]
  )
  # eggs
  # AB
  x[5] <- .5 * sum(
      pop[1],
      pop[2] * (1 / 2),
      pop[3] * (1 / 2),
      pop[4] * (1 / 2) * (1 - r),
      pop[5] * (1 / 2),
      pop[7] * (1 / 2) * r,
      pop[9] * (1 / 2),
      pop[10] * (1 / 2) * r,
      pop[13] * (1 / 2) * (1 - r)
    )
  # Ab
  x[6] <- .5 * sum(
      pop[6],
      pop[8] * (1 / 2),
      pop[14] * (1 / 2)
    )
  # aB
  x[7] <- .5 * sum(
      pop[11],
      pop[12] * (1 / 2),
      pop[15] * (1 / 2)
    )
  # ab
  x[8] <- .5 * pop[16]
  # Ab'
  x[9] <- 0.5 * sum(
    pop[2] * (1 / 2),
    pop[4] * (1 / 2) * r,
    pop[5] * (1 / 2),
    pop[7] * (1 / 2) * (1 - r),
    pop[10] * (1 / 2) * (1 - r),
    pop[13] * (1 / 2) * r
  )
  # a'B
  x[10] <- 0.5 * sum(
    pop[3] * (1 / 2),
    pop[4] * (1 / 2) * r,
    pop[7] * (1 / 2) * (1 - r),
    pop[9] * (1 / 2),
    pop[10] * (1 / 2) * (1 -r),
    pop[13] * (1 / 2) * r
  )
  # ab'
  x[11] <- 0.5 * sum(
    pop[12] * (1 / 2),
    pop[15] * (1 / 2)
  )
  # a'b
  x[12] <- 0.5 * sum(
    pop[8] * (1 / 2),
    pop[14] * (1 / 2)
  )
  # a'b'
  x[13] <- 0.5 * sum(pop[4] * (1 / 2) * (1 - r),
                     pop[7] * (1 / 2) * r,
                     pop[10] * (1 / 2) * r,
                     pop[13] * (1 / 2) * (1 - r)
  )
  # now we return this vector
  return(x)
}

# this fnx takes a gamete pool and reconstitutes a
makePop <- function(gametes, N) {
  # first we calculate probs based on frequency of gametes in pool
  # then we just draw from a multinomial distribution
  # in all cases crosses are listed sire by dame
  
  g <- gametes/sum(gametes)
  z.freq <- c(
    #-1-ABxAB
    g[1] * g[5],
    #-2-ABAb
    g[1] * g[6] + g[1] * g[9],
    #-3-ABaB
    g[1] * g[7] + g[1] * g[10],
    #-4-ABab
    g[1] * g[8] + g[1] * g[11] + g[1] * g[12] + g[1] * g[13],
    #-5-AbAB
    g[2] * g[5],
    #-6-AbAb
    g[2] * g[6],
    #-7-AbaB
    g[2] * g[7] + g[2] * g[10],
    #-8-Abab
    g[2] * g[8] + g[2] * g[12],
    #-9-aBAB
    g[3] * g[5],
    #-10-aBAb
    g[3] * g[6] + g[3] * g[9],
    #-11-aBaB
    g[3] * g[7],
    #-12-aBab
    g[3] * g[8] + g[3] * g[11],
    #-13-abAB
    g[4] * g[5],
    #-14-abAb
    g[4] * g[6],
    #-15-abaB
    g[4] * g[7],
    #-16-abab
    g[4] * g[8]
  )
  as.vector(rmultinom(1:16, N, prob = z.freq))
}

# this fnx takes gametes table and assigns mutations
Mutation <- function(gametes, gain, loss) {
  x <- gametes
  x <- round(x)
  # male losses
  #AB -> aB
  z <- sum(runif(x[1]) < loss)
  x[3] <- x[3] + z
  x[1] <- x[1] - z
  #AB -> Ab
  z <- sum(runif(x[1]) < loss)
  x[2] <- x[2] + z
  x[1] <- x[1] - z
  #Ab -> ab
  z <- sum(runif(x[2]) < loss)
  x[4] <- x[4] + z
  x[2] <- x[2] - z
  #aB -> ab
  z <- sum(runif(x[3]) < loss)
  x[4] <- x[4] + z
  x[3] <- x[3] - z
  # female losses
  #AB -> aB
  z <- sum(runif(x[5]) < loss)
  x[7] <- x[7] + z
  x[5] <- x[5] - z
  #AB -> Ab
  z <- sum(runif(x[5]) < loss)
  x[6] <- x[6] + z
  x[5] <- x[5] - z
  #Ab -> ab
  z <- sum(runif(x[6]) < loss)
  x[8] <- x[8] + z
  x[6] <- x[6] - z
  #aB -> ab
  z <- sum(runif(x[7]) < loss)
  x[8] <- x[8] + z
  x[7] <- x[7] - z
  #Ab' -> ab'
  z <- sum(runif(x[9]) < loss)
  x[11] <- x[11] + z
  x[9] <- x[9] - z
  #a'B -> a'b
  z <- sum(runif(x[10]) < loss)
  x[12] <- x[12] + z
  x[10] <- x[10] - z
  # male gains
  #Ab -> AB
  z <- sum(runif(x[2]) < gain)
  x[1] <- x[1] + z
  x[2] <- x[2] - z
  #aB -> AB
  z <- sum(runif(x[3]) < gain)
  x[1] <- x[1] + z
  x[3] <- x[3] - z
  #ab -> Ab
  z <- sum(runif(x[4]) < gain)
  x[2] <- x[2] + z
  x[4] <- x[4] - z
  #ab -> aB
  z <- sum(runif(x[4]) < gain)
  x[3] <- x[3] + z
  x[4] <- x[4] - z
  # female gains
  #Ab -> AB
  z <- sum(runif(x[6]) < gain)
  x[5] <- x[5] + z
  x[6] <- x[6] - z
  #aB -> AB
  z <- sum(runif(x[7]) < gain)
  x[5] <- x[5] + z
  x[7] <- x[7] - z
  #ab -> Ab
  z <- sum(runif(x[8]) < gain)
  x[6] <- x[6] + z
  x[8] <- x[8] - z
  #ab -> aB
  z <- sum(runif(x[8]) < gain)
  x[7] <- x[7] + z
  x[8] <- x[8] - z
  #Ab' -> AB
  z <- sum(runif(x[9]) < gain)
  x[5] <- x[5] + z
  x[9] <- x[9] - z
  #a'B -> AB
  z <- sum(runif(x[10]) < gain)
  x[5] <- x[5] + z
  x[10] <- x[10] - z
  #ab' -> Ab'
  z <- sum(runif(x[11]) < gain)
  x[9] <- x[9] + z
  x[11] <- x[11] - z
  #ab' -> aB
  z <- sum(runif(x[11]) < gain)
  x[7] <- x[7] + z
  x[11] <- x[11] - z
  #a'b -> Ab
  z <- sum(runif(x[12]) < gain)
  x[6] <- x[6] + z
  x[12] <- x[12] - z
  #a'b -> a'B
  z <- sum(runif(x[12]) < gain)
  x[10] <- x[10] + z
  x[12] <- x[12] - z
  #a'b' - Ab'
  z <- sum(runif(x[13]) < gain)
  x[9] <- x[9] + z
  x[13] <- x[13] - z
  #a'b' - a'B
  z <- sum(runif(x[13]) < gain)
  x[10] <- x[10] + z
  x[13] <- x[13] - z
  return(x)
}

# this function monitors the allele of interest
trackAllele <- function(results, allele) {
  x <- vector()
  for (j in 1:ncol(results)) {
    if (allele == "A" | allele == "a") {
      x[j] <-
        ((2 * sum(results[c(1, 2, 5, 6), j]) + sum(results[c(3, 4, 7, 8, 9, 10, 13, 14), j]))) / (2 *
                                                                                                    sum(results[, j]))
    }
    if (allele == "B" | allele == "b") {
      x[j] <-
        ((2 * sum(results[c(1, 3, 9, 11), j]) + sum(results[c(2, 4, 5, 7, 10, 12, 13, 15), j]))) / (2 *
                                                                                                      sum(results[, j]))
    }
  }
  if (allele == "b" | allele == "a")
    x <- 1 - x
  return(x)
}

# this fnx check to see if we A,a,B, or b are present at end of sim
checkAllele <- function(results, allele) {
  # get final result
  final <- results[, ncol(results)]
  # get all alleles
  x <-
    unlist(strsplit(as.character(row.names(results))[final != 0], split = ""))
  # test
  return(allele %in% x)
}

modelMedea <- function(max, N, pop.init, s, sims, loss, gain, bottleneck, bottleneck.freq, mon.both, fate, mutation, verbose){
  #################
  ################# Here is the actual function
  #################
  ## SETUP FOR RECORDING RESULTS
  if(mon.both == T){
    final.result <- list(length=sims)
  }
  if(fate == T){
    final.result <- matrix(,2,sims)
    row.names(final.result) <- c("A", "B")
  }
  if((fate+mon.both) == 2){
    stop("Both fate and mon.both can not be TRUE choose to track either both alleles through all generations or choose to track only the presence of Medea factors at the end of the simulation.")
  }
  if((fate+mon.both) == 0){
    stop("Set either fate or mon.both to true to return a result from the simulation.")
  }
  ## SETUP FOR SIMULATION
  # the population table
  pop <- vector(length=16, mode="numeric")
  names(pop) <- c("ABAB", "ABAb", "ABaB", "ABab", "AbAB", "AbAb", "AbaB", "Abab", 
                  "aBAB", "aBAb", "aBaB", "aBab", "abAB", "abAb", "abaB", "abab")
  # now our gamete table
  gametes <- vector(length=13, mode="numeric")
  names(gametes) <- c("AB", "Ab", "aB", "ab", "AB", "Ab", "aB", "ab", "Ab'", "a'B", "ab'", "a'b", "a'b'")
  # a matrix to store results of a single simulation
  results <- data.frame(matrix(, 16, max))
  row.names(results) <- names(pop)
  # loop for iterations
  for(i in 1:sims){
    # setup starting conditions
    results[, ] <- NA
    pop[1:16] <- pop.init
    # increment variable
    counter <- 1
    # record starting condition
    results[, counter] <- pop
    ## Simulate generations
    while(counter < max){
      if(verbose == T) print(counter)
      counter <- counter + 1
      pop[1:16] <- selection(pop, N, s)
      gametes[1:13] <- makeGametes(pop, r = 0.5, N = N)
      #current solution is too slow to slow to use:
      if(mutation == T){
        gametes[1:13] <- Mutation(gametes=gametes, gain=gain, loss=loss)
      }
      ##TODO CHANGE 7 TO A ARGUMENT
      if(counter %% bottleneck.freq == 0){
        Ne <- N - (bottleneck * N)
      }else{
        Ne <- N
      }
      pop[1:16] <- makePop(gametes, N = Ne)
      results[, counter] <- pop
    }
    # now we record the appropriate values from this sim
    if(mon.both == T){
      x <- rbind(trackAllele(results, allele = "A"),
                 trackAllele(results, allele = "B"))
      row.names(x) <- c("A", "B")
      colnames(x) <- 1:max
      final.result[[i]] <- x
    }
    if(fate == T){
      final.result[1:2, i] <- c(checkAllele(results, "A"), 
                                checkAllele(results, "B"))
    }
    paste("Simulation", i, "complete\n")
  }
  return(final.result)
}

