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