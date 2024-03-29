# Assignment 1 
# Stochastic Processes and Bayesian Inference
install.packages('expm')
library('expm')

writeLines("==================================")
writeLines("\n\tAssignment 2\n\tMarten Skogh\n")
writeLines("==================================")

# Excersize 1.35 ----------------------------
accidents <- function() {
   poisson_param <- runif(1,0,3)
   rpois(1, poisson_param)
}

number_days <- 100000
nbr_accidents <- replicate(number_days, accidents())

sprintf("Mean number of accidents per day: %f", mean(nbr_accidents))

sprintf("Variance in number of accidents per day: %f", var(nbr_accidents))

hist(nbr_accidents,
     breaks=(-0.5:(max(nbr_accidents)+0.5)),
     main='Frequency of Daily Accidents',
     xlab='Number of Accidents')

# Exercise 2.27 ------------------------------

# (a)
gamble <- function(k,n,p) {
   # Run until win or loss
   while (k > 0 && k < n) {
      r <- runif(1,0,1)
      if (r >= p) {
         k <- k + 1
      } else {
         k <- k - 1
      }
      #print(k)
   }
   if (k == 0) {
      0
   } else {
      1
   }
}

k <- 2   # Starting capital
n <- 5   # Win condition
p <- 1/2 # Win probability

trials <- 100000

simlist <- replicate(trials, gamble(k,n,p))
win_mean <- mean(simlist)

sprintf("Expected risk of ruin: %f", 1 - win_mean)

# (b)
init_state = c(0, 0, 1, 0, 0, 0)

trans_prob <- c(  1,   0,   0,   0,   0,   0,
                1/2,   0, 1/2,   0,   0,   0,
                  0, 1/2,   0, 1/2,   0,   0,
                  0,   0, 1/2,   0, 1/2,   0,
                  0,   0,   0, 1/2,   0, 1/2,
                  0,    0,   0,   0,   0,  1)

P <- matrix(trans_prob, nrow=6, byrow=TRUE)
rownames(P) <- c("0", "1", "2", "3", "4", "5")
colnames(P) <- c("0", "1", "2", "3", "4", "5")

n <- 100

final_prob <- init_state %*% (P %^% n)

sprintf("Longterm probability for ruin: %f", final_prob[1])

# (c)

exact_prob = (n - k)/n
sprintf("Exact probability: %f", exact_prob)


# Exercise 3.52 -------------------------------

# (a) See example 3.28.

# You can't really be in the states 2, 5 and 8 so we don't need those states 
# represented in our transition matrix. 

trans_prob <- c(0, 1, 1, 1, 0, 1, 0,
                0, 0, 2, 1, 0, 1, 0,
                0, 0, 1, 1, 1, 1, 0,
                0, 0, 1, 1, 1, 1, 0,
                0, 0, 0, 1, 1, 1, 1,
                0, 0, 0, 1, 0, 2, 1,
                0, 0, 0, 0, 0, 0, 4) / 4

P <- matrix(trans_prob, ncol=7, byrow=TRUE)
rownames(P) <- c("0", "1", "3", "4", "6", "7", "9")
colnames(P) <- c("0", "1", "3", "4", "6", "7", "9")

Q <- P[1:6,1:6]
F <- solve(diag(6) - Q)
a <- F %*% rep(1,6)

colnames(a) <- c("9")

sprintf("Expected number of steps to win finish the game (0->9): %f", (a[1]))

# (b) See example 3.29.

# We make 3 an absorbing state since we have two possible outcomes: either we 
# pass through 3 or we don't.

P2 <- P[3:7, 3:7] # No prob to land on 0,1
P2 <- P2[c(2,3,4,1,5), c(2,3,4,1,5)] # Reorder to block matrix
P2[4, 1:5] <- c(0, 0, 0, 1, 0) 

Q <- P2[1:3, 1:3]
R <- P2[1:3, 4:5]
F <- solve(diag(3) - Q)
a <- F %*% R

#colnames(a) <- c("3", "9")

sprintf("Probability of passing through 3 if starting on 6: %f", a["6","3"])
sprintf("Probability of not passing through 3 if starting on 6: %f", a["6","9"])