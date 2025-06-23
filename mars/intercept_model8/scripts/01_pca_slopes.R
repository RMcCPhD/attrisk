
mu <- c(-1.05, -1.79)
sd_lograte <- 0.73
sd_shape <- 2.26
cor_val <- -0.63

Sigma <- matrix(c(sd_lograte^2, cor_val * sd_lograte * sd_shape,
                  cor_val * sd_lograte * sd_shape, sd_shape^2), nrow = 2)

eigen_decomp <- eigen(Sigma)
vec1 <- eigen_decomp$vectors[, 1]   # First principal component

# Normalized direction vector (unit length)
vec1
