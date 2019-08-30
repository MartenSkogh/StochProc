# Assignment 1 
# Stochastic Processes and Bayesian Inference

writeLines("==================================")
writeLines("\n\tAssignment 2\n\tMarten Skogh\n")
writeLines("==================================")

# 1. (a) ----------------

population <- c(1, 1, 2, 5, 7)
offspring <- c(1, 2, 3, 2, 0, 0, 1, 4,2)

# 1. (b) ----------------
# Using the extinction probability formulas as used in example 4.7 and 4.8.

offspring_dist <- c(2, 2, 3, 1, 1) / 9

polynomial <- function(s, dist) {
    total <- 0
    pow <- 0

    for (d in dist) {
        total <- total + s**pow * d 
        pow <- pow + 1
    }

    total
}

gamma_pgf <- function(s, lambda) {
    exp(lambda * (s - 1))
}

extinction_func <- function(s, func, ...) {
    func(s, ...) - s
}

newton_raphson <- function(func, args=NULL, tol=1E-12, x0=0, N=1000) {
    h <- 0.001
    n <- 1
    x <- x0
    #p <- numeric(N)

    y <- function(x) {
        do.call(func, c(x, args))
    }

    while (n <= N) {
        dy <- (y(x0 + h) - y(x0)) / h
        x <- (x0 - (y(x0) / dy))
        #p[n] <- x
        n <- n + 1
        
        if (abs(x - x0) < tol) {
            x0 <- x
            break
        } 
        x0 <- x
    }
    
    if (n >= N) {
        print("Warning: Did not converge!")
    }

    x0
    #p[1:(i-1)]
}

branch_extinction_prob <- newton_raphson(extinction_func, args=c(gamma_pgf, 15/9), x0=0) 
total_extinction_prob <- branch_extinction_prob**7

sprintf("The probability that single branch will go extinct is: %f", branch_extinction_prob)
sprintf("The probability that all branches will go extinct is: %f", total_extinction_prob)

# 1. (c) ----------------
integrand = function(lambda) {
    prior <- function(lambda) {
        newton_raphson(extinction_func, args=c(gamma_pgf, lambda), x0=0)
    }

    integrand_value <- function(lambda) {
        prior(lambda)^7 * dgamma(lambda, 15, rate=9)
    }

    sapply(lambda, integrand_value)
}

tot_ext_prob <- integrate(integrand, lower=0, upper=Inf)$value

sprintf(sprintf("The probability that all branches will go extinct is: %f", tot_ext_prob)
)

# 1. (d) ----------------
run_branch <- function(maxsteps, offspring_prob) {
    Z <- vector(,maxsteps)
    Z[1] <- 1
    n <- 1

    new_generation <- function(previous_generation, offspring_prob) {
        if (previous_generation <= 0) {
            return(0)
        }
        new_generation_pop <- 0
        for (i in 1:previous_generation) {
            r <- runif(1)
            a <- 0
            n <- 1
            max_n = 100
            while(r > offspring_prob(n - 0.5) && n < max_n) {
                n <- n + 1
            }
            #print(sprintf("n = %f, r = %f", n - 1, r))
            new_generation_pop <- new_generation_pop + n - 1
        }
        return(new_generation_pop)
    }

    while (n < maxsteps) {
        n <- n + 1
        Z[n] <- new_generation(Z[n-1], offspring_prob)
    }

    return(Z)
}

nbr_branches <- 10000
max_gen <- 5
branches <- array(0, dim=c(nbr_branches,max_gen))
nbr_extinct <- 0
for (i in 1:nbr_branches) {
    branches[i,] <- run_branch(max_gen, function(n) {ppois(n, rgamma(1, 15, rate=9))})
    if (branches[i,max_gen] == 0) {
        nbr_extinct <- nbr_extinct + 1
    }
    if (i %% 1000 == 0) { print(sprintf("%.2f %%", i / nbr_branches * 100)) }
}

#branches[,max_gen]

sum(branches[,max_gen] == 0)

print(sprintf("%f",nbr_extinct/nbr_branches))


# 2. (a) -----------------------
data <- read.delim('C:\\Users\\marte\\Documents\\Stochastic Processes\\StochProc\\Regressiondata.txt', header=FALSE, sep=' ')
x <- data[,1]
y <- data[,2]

plot(x,y)

# 2. (b) -----------------------
epsilon <- rnorm(23,0,1)
theta <- c(6, 2 * pi * 0.8, 1)

model <- function(x, theta) {
    theta[1] * sin(x / theta[2]) + theta[3] * runif(1)
}
plot(x,y)
lines(x, model(x, theta), col="red")

# 2. (c) -----------------------

log_likliehood <- function(x, y, theta) {
    sum(-log(sqrt(2 * pi) * theta[3]) - (y - model(x, theta))**2 / (2 * theta[3]**2))
}

range = seq(0,10, length=1000)
ans <- vector(, length(range))
for (i in 1:length(range)) {
    theta <- c(range[i], 5.02, 1.15)
    ans[i] <- log_likliehood(x, y, theta)
}
plot(range,ans)
which.max(ans)

range = seq(0, 10, length=1000)
ans <- vector(, length(range))
for (i in 1:length(range)) {
    theta <- c(5.08, range[i], 1.15)
    ans[i] <- log_likliehood(x, y, theta)
}
plot(range,ans)
which.max(ans)

range = seq(0, 10, length=1000)
ans <- vector(, length(range))
for (i in 1:length(range)) {
    theta <- c(5.08, 5.02, range[i])
    ans[i] <- log_likliehood(x, y, theta)
}
plot(range,ans)
which.max(ans)

# 2. (e) -----------------------
metropolis_hastings <- function(steps) {
    accepted <- 0
    theta <- array(0, dim = c(steps ,3))
    theta[1,] <- c(6, 2 * pi * 0.8, 1) # Initial value
    
    proposal_function <- function(theta) {abs(theta + rnorm(3, 0, 0.1))}

    for (i in 2:steps) {
        # Generate
        proposal_value <- proposal_function(theta[i - 1,])
        # Calculate
        p <- log_likliehood(x, y, proposal_value) / log_likliehood(x, y, theta[i - 1,])
        if (p > runif(1)) { # Accept
            theta[i,] <- proposal_value
            accepted <- accepted + 1
        } else { # Reject
            theta[i,] <- theta[i - 1,]
        }
    }
    return(list(accepted/steps, theta))
}

steps <- 1000
a <- metropolis_hastings(steps)
acc_rate <- a[[1]]
theta_trace <- a[[2]] 

plot(theta_trace[,1], type='l', main='Theta_1')
plot(theta_trace[,2], type='l', main='Theta_2')
plot(theta_trace[,3], type='l', main='Theta_3')


# 2. (f) ------------------
pnorm(0, model(15, c(theta_trace[steps,1], theta_trace[steps,2], 0)), theta_trace[steps,3])