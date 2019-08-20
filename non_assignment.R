install.packages('expm')
library('expm')

# Example 6.2 -------------------------
ppois(200,2*100) 

# Example 6.5 -------------------------
# i.
pgamma(60,4,1/15) - pgamma(55,4,1/15)

# ii.
integrand <- function(t) {
    (1/15)**3 * t**3 * exp(-t/15) / 2 
}
1/pgamma(60,3,1/15) * integrate(integrand, lower=0, upper=60)$value

# Example 6.6 -------------------------
# i. 
male_prob <- 108/(100 + 108)
lambda <- 2
male_param <- lambda * male_prob
female_param <- lambda * (1 - male_prob)

eight_hour_expect <- female_param * 8
std_deviation <- sqrt(eight_hour_expect)

sprintf("Expected number of female births during 8 hours: %f", eight_hour_expect)
sprintf("Standard deviation of female births during 8 hours: %f", std_deviation)

# ii.
dpois(0, male_param * 3) * (1 - dpois(0, female_param * 3))

# iii. 
dbinom(2, size=5, prob=male_prob)

# Example 6.7 -------------------------
# No numerics

# Example 6.8 -------------------------
# i.
pgamma(3,60,20) - pgamma(2.9,60,20)

# ii.
1 - (2.9/3)**60

# Example 6.11 -------------------------
dpois(5, 4 * pi * 1/2)

# Example 6.13 -------------------------
intensity <- function(t) {
    select_intensity <- function(t) {
        if (t >= 0 && t <= 1) {
            100 + 100 * t
        } else if (t > 1 && t <= 3) {
            200
        } else if (t > 3 && t < 4) {
            500 - 100 * t
        } else {
            0
        }
    }

    vec_intensity <- Vectorize(select_intensity)
    vec_intensity(t)
}

expected_nbr_students <- integrate(intensity, lower=0.5, upper=2.5)$value
1 - ppois(399,expected_nbr_students)

# Example 7.25 ------------------------

# 1. 
rate <- 6
mu <- 2
c <- 5

1 / (sum(((rate / mu)**(0:(c - 1)) / factorial(0:(c - 1)))) 
    + (rate / mu)**5 / factorial(5) * (1 / (1 - (rate / (mu * c)))))

# 2. 


# Example 7.27
death_rate <- 2

Q <- c(-1,  1,  0,  0,  0,  0,
        2, -4,  2,  0,  0,  0,
        0,  2, -5,  3,  0,  0,
        0,  0,  2, -6,  4,  0,
        0,  0,  0,  2, -7,  5,
        0,  0,  0,  0,  2, -2)

Q <- matrix(Q, nrow=6, byrow=TRUE)
colnames(Q) <- as.character(0:5)
rownames(Q) <- as.character(0:5)

max_q <- max(abs(Q))

R <- Q/max_q + diag(6)

error <- Inf
limit <- 0.5E-4
N <- 0

while (error >= limit) {
    N <- N + 1
    error <- 1 - ppois(N, max(abs(Q)) * 1.5)
}

partial_sums <- array(0, c(6,6,N+1))
prob_1_5 <- array(0, c(6,6))

colnames(prob_1_5) <- as.character(0:5)
rownames(prob_1_5) <- as.character(0:5)

for (k in 0:N) {
    #partial_sums[,,k+1] <- (R ** k) * dpois(k, max_q * 1.5)
    #prob_1_5 <- prob_1_5 + partial_sums[,,k+1]
    prob_1_5 <- prob_1_5 + (R %^% k) * dpois(k, max_q * 1.5)
}

comparison <- expm(1.5 * Q)
diff <- comparison - prob_1_5

# Example 8.2 ------------------------
pnorm(2, 0, sqrt(3))

# Example 8.8 ------------------------
(4 / qnorm(0.95, 0, 1))**2


# Example 8.10 -----------------------
pnorm(0.6,0,1) - pnorm(-1.4,0,1)