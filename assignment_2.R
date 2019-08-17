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

newton_raphson <- function(func, args=NULL, tol=1E-12, x0=0, N=100) {
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
        if (abs(x - x0) < tol) break
        x0 <- x
    }
    
    x0
    #p[1:(i-1)]
}

branch_extinction_prob <- newton_raphson(extinction_func, args=c(gamma_pgf, 15/9), x0=0) 
total_extinction_prob <- branch_extinction_prob**7

sprintf("The probability that single branch will go extinct is: %f", branch_extinction_prob)
sprintf("The probability that all branches will go extinct is: %f", total_extinction_prob)

# 1. (c) ----------------

integrand <- function(lambda) {
    shape <- 15
    rate <- 9
    lambda * dgamma(lambda, shape, rate=rate)
}

integral <- integrate(integrand, lower=0, upper=Inf)