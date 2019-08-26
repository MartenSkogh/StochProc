# Assignment 3 
# Stochastic Processes and Bayesian Inference

print("Assignment 3")
print("Mårten Skogh")
print("")

par(mfrow=c(2,1))

# 1. (a) ---------------------
print("===== Question 1. (a) =====")
pois_param <- 36
area <- 0.4 * 0.4

prob_6_or_more <- 1 - ppois(6, pois_param * area)

print(sprintf("Probability of six or more trees in area %f is: %f", area, prob_6_or_more))


# 1. (b) -------------------
print("===== Question 1. (b) =====")
# Non-overlapping areas are uncorrelated

NB <- 0:4
NAC <- 4 - NB

area_B <- 0.2**2
area_AC <- 0.4**2 - area_B

probs <- ppois(NAC, pois_param * area_AC)**2 * ppois(NB, pois_param * area_B)
prob <- sum(probs)

print(sprintf("Probability of 4 trees in both areas is: %f", prob))

# 1. (c) -------------------
print("===== Question 1. (c) =====")
nbr_trees <- rpois(1, pois_param)
x <- runif(nbr_trees)
y <- runif(nbr_trees)

plot(x, y, main="Constant posterior")

# 1. (d) -------------------
print("===== Question 1. (d) =====")
nbr_trees <- rpois(1, rgamma(1, 36, 1))
x <- runif(nbr_trees)
y <- runif(nbr_trees)

plot(x, y, main="Gamma posterior")

# 1. (e) -------------------
print("===== Question 1. (e) =====")
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
print("===== Question 1. (f) =====")

# 2. (a) -------------------
print("===== Question 2. (a) =====")
q <- c(1/5, 1, 1/2)

p <- c(  0, 0.5, 0.5,
       0.5,   0, 0.5,
       0.5, 0.5,   0)

p <- matrix(p, ncol = 3, byrow=TRUE)
Q <- p * q - diag(3) * rowSums(p * q)

print("Q:")
print(Q)
v <- solve(rbind(t(Q[,1:2]),c(1,1,1)), c(0,0,1))
print("Stationary solution:")
print(v)


# 2. (b) -------------------
print("===== Question 2. (b) =====")
states <- c(1, 3, 2, 3, 1, 2, 1, 3, 1, 2)
durations <- c(6.83, 4.01, 1.63, 0.44, 5.11, 0.29, 2.87, 1.30, 4.76, 1.92)

nbr_times_in_state <- rep(0, 3)
tot_stay <- rep(0, 3)
nbr_trans <- matrix(rep(0, 9), ncol = 3, byrow = TRUE)

for (i in 1:length(states)) {
    nbr_times_in_state[states[i]] <- nbr_times_in_state[states[i]] + 1
    tot_stay[state[i]] <- tot_stay[state[i]] + durations[i]
    if (i < length(state)) {
        nbr_trans[states[i],states[i+1]] <- nbr_trans[states[i], states[i+1]] + 1
    }
}

mean_stay <- tot_stay / nbr_times_in_state

print(nbr_trans)
print(mean_stay)
