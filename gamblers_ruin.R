gamble <- function(k,n,p) {
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

k <- 20
n <- 60
p <- 1/2
trials <- 1000
simlist <- replicate(trials, gamble(k,n,p))
mean(simlist)