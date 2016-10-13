library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

# Returns a probability vector proportional to (exp(delta, 1, ..., 1))
Pr <- function(n, delta) {
    pi <- rep(1, n)
    pi[1] <- exp(delta)
    pi <- pi / sum(pi)
}

# Find the selective test power for multinomial distribution
# m = number of votes
# n = number of candidates
# distribution is proportional to (exp(delta), 1, ..., 1)
# alpha = level of the test
# num.sample = number of samples to take for simulation purposes
SelectivePower <- function(m, n, delta, alpha, num.sample) {
    # Probability vector of multinomial distribution
    pi <- Pr(n, delta)
    # Generates the samples
    x <- rmultinom(num.sample, m, pi)
    # sec.x is max(Xj; j>1), i.e. X2
    if (n == 2) {
        sec.x <- x[2, ]
    } else {
        sec.x <- apply(x[2:n, ], 2, max)
    }
    # n.x is X1 + X2
    n.x <- sec.x + x[1, ]
    rejected <- length(which(pbinom(sec.x, n.x, 0.5) * 2 <= alpha))
    return(rejected / num.sample)
}

# Find the test power of Gupta & Nagel for multinomial distribution
# m = number of votes
# n = number of candidates
# delta is a vector listing all the deltas we want to plot
# alpha = level of the test
# num.sample = number of sampels to take for simulation purposes
GnPower <- function(m, n, delta, alpha, num.sample) {
    # Initialization; d.star is the biggest d
    d.star <- -Inf
    for (r in 2:n) {
        # Probability vector (1/r, ..., 1/r)
        pi <- rep(1 / r, r)
        # Generates the sample
        x <- rmultinom(num.sample, m, pi)
        # Takes the maximum of each column
        max.x <- apply(x, 2, max)
        # d for a given r is the upper alpha quantile, rounded up
        d <- quantile(max.x - x[1, ], probs=(1 - alpha), type=1, names=FALSE)
        if (d.star < d) {
            d.star <- d
        }
    }
    
    # Initialization; power is a vector of powers corresponding to the deltas
    power <- rep(0, length(delta))
    for (i in 1:length(delta)) {
        # Probability vector
        pi <- Pr(n, delta[i])
        # Generates the sample
        x <- rmultinom(num.sample, m, pi)
        # sec.x is max(Xj; j>1), i.e. X2
        if (n == 2) {
            sec.x <- x[2, ]
        } else {
            sec.x <- apply(x[2:n, ], 2, max)
        }
        rejected <- length(which(x[1, ] - sec.x > d.star))
        power[i] <- rejected / num.sample
    }
    
    return(power)
}

# Plots the power of the selective test vs Gupta & Nagel for multinomial
# distribution
# m = number of votes, i.e. X1 + ... + Xn
# n = number of candidates
# alpha = level of the test
# num.sample = number of samples to take for simulation purposes
PlotMultinomialPower <- function(m, n, alpha, num.sample) {
    # delta ranges from 0 to 2
    delta <- (0:200) / 100
    sel <- sapply(delta, function(d) SelectivePower(m, n, d, alpha / (n - 1) * n, num.sample))
    gn <- GnPower(m, n, delta, alpha, num.sample)
    power <- melt(data.frame(delta, sel, gn), id.vars = "delta")
    return(ggplot(data = power, aes(x = delta,
                                  y = value,
                                  group = variable,
                                  colour = variable,
                                  linetype = variable)) +
               geom_line() +
               coord_cartesian(ylim = c(0, 1)) +
               theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.title = element_blank(),
                     legend.position = "none",
                     plot.title = element_text(size = 10, face = "bold")) +
               ggtitle(sprintf("m = %d, n = %d", m, n)))
}

# Main code
num.sample <- 10000
alpha <- 0.05
p1 <- PlotMultinomialPower(50, 2, alpha, num.sample)
p2 <- PlotMultinomialPower(50, 5, alpha, num.sample)
p3 <- PlotMultinomialPower(50, 10, alpha, num.sample)
p4 <- PlotMultinomialPower(250, 10, alpha, num.sample)
p5 <- PlotMultinomialPower(250, 25, alpha, num.sample)
p6 <- PlotMultinomialPower(250, 50, alpha, num.sample)

# Saving the plot
png("multinomial-power.png", width = 7, height = 5, units = 'in', res = 300)
grid.arrange(arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3,
                         left = textGrob("Power", rot = 90, vjust = 1),
                         bottom = textGrob(expression(delta))))
dev.off()