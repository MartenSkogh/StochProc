# Assignment 3 
# Stochastic Processes and Bayesian Inference
install.packages("mcsm")
library("mcsm")

writeLines("==================================")
writeLines("\n\tAssignment 3\n\tMårten Skogh\n")
writeLines("==================================")

#par(mfrow=c(2,1))

# 1. (a) ---------------------
writeLines("\n===== Question 1. (a) =====\n")
pois_param <- 36
area <- 0.4 * 0.4

prob_6_or_more <- 1 - ppois(6, pois_param * area)

print(sprintf("Probability of six or more trees in area %f is: %f", area, prob_6_or_more))


# 1. (b) -------------------
writeLines("\n===== Question 1. (b) =====\n")
# Non-overlapping areas are uncorrelated

NB <- 0:4
NAC <- 4 - NB

area_B <- 0.2**2
area_AC <- 0.4**2 - area_B

probs <- ppois(NAC, pois_param * area_AC)**2 * ppois(NB, pois_param * area_B)
prob <- sum(probs)

print(sprintf("Probability of 4 trees in both areas is: %f", prob))

# 1. (c) -------------------
writeLines("\n===== Question 1. (c) =====\n")
nbr_trees <- rpois(1, pois_param)
x <- runif(nbr_trees)
y <- runif(nbr_trees)

plot(x, y, main="Constant posterior")

# 1. (d) -------------------
writeLines("\n===== Question 1. (d) =====\n")
nbr_trees <- rpois(1, rgamma(1, 36, 1))
x <- runif(nbr_trees)
y <- runif(nbr_trees)

dev.new()
plot(x, y, main="Gamma posterior")

# 1. (e) -------------------
writeLines("\n===== Question 1. (e) =====\n")
min_dist <- function(x, y) {
    X <- rep(x, times = length(x))
    X <- matrix(X, nrow = length(x), byrow=TRUE)
    
    Y <- rep(y, times = length(y))
    Y <- matrix(Y, nrow = length(y), byrow=TRUE)
    
    dist <- sqrt(sweep(-X, 1, x, '+')**2 + sweep(-Y, 1, y, '+')**2)
    dist[dist == 0] <- Inf
    sapply(1:length(x), function(i) {min(dist[,i])})

}

trials <- 10000
sum_avg <- 0
for (i in 1:trials) {
    nbr_trees <- rpois(1, rgamma(1, 36, 1))
    x <- runif(nbr_trees)
    y <- runif(nbr_trees)
    sum_avg <- sum_avg + sum(min_dist(x,y)) / length(x)
}
avg_min_dist <- sum_avg / trials

print(sprintf("Average minimum distance is: %f", avg_min_dist))

# 1. (f) -------------------
writeLines("\n===== Question 1. (f) =====\n")

# 2. (a) -------------------
writeLines("\n===== Question 2. (a) =====\n")
q <- c(1/5, 1, 1/2)

p <- c(  0, 0.5, 0.5,
       0.5,   0, 0.5,
       0.5, 0.5,   0)

p <- matrix(p, ncol = 3, byrow=TRUE)
Q <- p * q - diag(3) * rowSums(p * q)

writeLines("Q:")
print(Q)
v <- solve(rbind(t(Q[,1:2]),c(1,1,1)), c(0,0,1))
writeLines("\nStationary solution:")
print(v)


# 2. (b) -------------------
writeLines("\n===== Question 2. (b) =====\n")
states <- c(1, 3, 2, 3, 1, 2, 1, 3, 1, 2)
durations <- c(6.83, 4.01, 1.63, 0.44, 5.11, 0.29, 2.87, 1.30, 4.76, 1.92)

nbr_times_in_state <- rep(0, 3)
tot_stay <- rep(0, 3)
counts <- matrix(rep(0, 9), ncol = 3, byrow = TRUE)

for (i in 1:length(states)) {
    nbr_times_in_state[states[i]] <- nbr_times_in_state[states[i]] + 1
    tot_stay[states[i]] <- tot_stay[states[i]] + durations[i]
    if (i < length(states)) {
        counts[states[i],states[i+1]] <- counts[states[i], states[i+1]] + 1
    }
}


writeLines("Count matrix:")
print(counts)

q_posterior_means <- tot_stay / nbr_times_in_state
P_posterior_means <- matrix(
                            c(  0, 1/2, 1/2,
                              1/2,   0, 1/2,
                              5/8, 3/8,   0),
                            ncol = 3,
                            byrow = TRUE)

writeLines("\nP posterior:")
print(P_posterior_means)

get_trans_matrix <- function(posterior) {
    P <- matrix(c(  0, 0, 0,
                    0, 0, 0,
                    0, 0, 0),
                ncol = 3,
                byrow = TRUE)

     P[1,] <- rdirichlet(1, posterior[1,])
     P[2,] <- rdirichlet(1, posterior[2,])
     P[3,] <- rdirichlet(1, posterior[3,])

     return(P)
}

get_q_vector <- function(tot_stay, nbr_times_in_state) {
    q <- c(0,0,0)
    q[1] <- rgamma(1, tot_stay[1], nbr_times_in_state[1])
    q[2] <- rgamma(1, tot_stay[2], nbr_times_in_state[2])
    q[3] <- rgamma(1, tot_stay[3], nbr_times_in_state[3])
    return(q)
}

q <- get_q_vector(tot_stay, nbr_times_in_state)
P <- get_trans_matrix(P_posterior_means)
Q <- P * q - diag(3) * rowSums(P * q)

writeLines("\nq-vector:")
print(q)
writeLines("\nTransition matrix:")
print(P)
writeLines("\nInfinitesimal generator matrix:")
print(Q)
