
# v01 copied from .../Research/2021/Confidence interval combination/Code/interval_comb_v17d_ExampleD.m 
#     where it constructs pictures for proposition 2 of Magnus and Vasnev (2023) 
#     "On the uncertainty of a combined forecast: The critical role of correlation" in version 17d
# v02 created for revision in January 2025
# v03 is a cleaned version for GitHub with notation matching the paper
#     Magnus and Vasnev (2025) "The role of data and priors in estimating
#     climate sensitivity"
# v04 is translated to R 
#     and added calculations from IPCC_prior_calculations_v7_revision_github.xlsx

# Clear environment
rm(list = ls())

# Define functions
mu_hat_sigma_hat <- function(r, sizeN, V_0, fcst_point) {
  A <- matrix(1, sizeN, sizeN) / sizeN
  alpha <- (sizeN - 1) * r + 1
  beta <- 1 - r
  P_inv <- (1 / alpha) * A + (1 / beta) * (diag(sizeN) - A)
  V_temp <- diag(1 / sqrt(diag(V_0)))
  vones <- rep(1, sizeN)
  num <- t(vones) %*% V_temp %*% P_inv %*% V_temp %*% fcst_point
  denom <- t(vones) %*% V_temp %*% P_inv %*% V_temp %*% vones
  mu_hat <- as.numeric(num / denom)
  err_hat <- fcst_point - mu_hat
  sigma2_hat <- as.numeric(t(err_hat) %*% V_temp %*% P_inv %*% V_temp %*% err_hat) / sizeN
  V_inv <- V_temp %*% ((A / alpha) + ((diag(sizeN) - A) / beta)) %*% V_temp
  i_Vinv_i <- as.numeric(t(vones) %*% V_inv %*% vones)
  tau2_hat <- sigma2_hat / i_Vinv_i
  list(mu_hat = mu_hat, sigma2_hat = sigma2_hat, tau2_hat = tau2_hat, i_Vinv_i = i_Vinv_i)
}

proposition2_lim <- function(x, v) {
  sizeN <- length(x)
  Eps <- 1e-5
  vones <- rep(1, sizeN)
  A <- matrix(1, sizeN, sizeN) / sizeN
  V0_sqrt_inv <- diag(1 / v)
  w <- V0_sqrt_inv %*% vones
  y <- V0_sqrt_inv %*% x
  C_wy <- t(w) %*% (diag(sizeN) - A) %*% y / sizeN
  C_ww <- t(w) %*% (diag(sizeN) - A) %*% w / sizeN
  C_xx <- t(x) %*% (diag(sizeN) - A) %*% x / sizeN
  C_vv <- t(v) %*% (diag(sizeN) - A) %*% v / sizeN
  C_yy <- t(y) %*% (diag(sizeN) - A) %*% y / sizeN
  if (sum(v == vones) == sizeN) {
    mu_hat_lim <- mean(x)
  } else {
    mu_hat_lim <- C_wy / C_ww
  }
  if ((sum(v == vones) == sizeN) && (sum((x / mean(x)) == vones) == sizeN)) {
    sigma2_hat_lim <- 0
  } else if ((sum(v == vones) < sizeN) &&
             (sum((x - cbind(vones, v) %*% solve(t(cbind(vones, v)) %*% cbind(vones, v)) %*% t(cbind(vones, v)) %*% x)^2) < Eps)) {
    sigma2_hat_lim <- C_xx / (sizeN * C_vv)
  } else {
    sigma2_hat_lim <- Inf
  }
  if (sum((x / mean(x)) == vones) == sizeN) {
    tau2_hat_lim <- 0
  } else if ((sum((x / mean(x)) == vones) < sizeN) && (sum(v == vones) < sizeN)) {
    tau2_hat_lim <- (C_ww * C_yy - C_wy^2) / (sizeN * C_ww^2)
  } else {
    tau2_hat_lim <- Inf
  }
  list(mu = mu_hat_lim, sigma2 = sigma2_hat_lim, tau2 = tau2_hat_lim)
}

# Input data
caseN <- 4  # 14 for IPCC5 data, 4 for IPCC6

if (caseN == 14) {
  b_0i <- c(0.68, 0.72, 0.72, 0.99, 0.97, 0.75, 0.89, 1.17, 0.93, 0.99, 1.09, 1.23, 1.15, 1.02)
  sigma_0i <- c(0.21, 0.33, 0.44, 0.18, 0.22, 0.52, 0.43, 0.08, 0.45, 0.40, 0.26, 0.17, 0.28, 0.47)
} else if (caseN == 4) {
  b_0i <- c(1.22, 1.03, 1.2, 1.05)
  sigma_0i <- c(0.36, 0.39, 0.61, 0.36)
} else {
  stop("no such case yet")
}

x <- b_0i
m <- caseN
vones <- rep(1, m)
tmp <- m / sum(sigma_0i^2)
V_0 <- diag(tmp * sigma_0i^2)
v <- sqrt(diag(V_0))

r_set <- seq(0, 0.99, by = 0.0001)
r_count <- length(r_set)

mu_hat_set <- numeric(r_count)
sigma2_hat_set <- numeric(r_count)
tau2_hat_set <- numeric(r_count)
i_Vinv_i_set <- numeric(r_count)

for (i in seq_len(r_count)) {
  res <- mu_hat_sigma_hat(r_set[i], m, V_0, x)
  mu_hat_set[i] <- res$mu_hat
  sigma2_hat_set[i] <- res$sigma2_hat
  tau2_hat_set[i] <- res$tau2_hat
  i_Vinv_i_set[i] <- res$i_Vinv_i
}

fcst_mean <- mean(x)
err_hat <- x - fcst_mean
s2 <- sum(err_hat^2) / m

lim_r <- 1
lim_res <- proposition2_lim(x, v)
lim_mu <- lim_res$mu
lim_sigma2 <- lim_res$sigma2
lim_tau2 <- lim_res$tau2

if (caseN == 14) {
  tmp <- (sigma2_hat_set - 0.2809)^2
} else if (caseN == 4) {
  tmp <- (sigma2_hat_set - 0.0729)^2
} else {
  stop("not implemented yet")
}
I <- which.min(tmp)
r_tmp <- r_set[I]

r_0 <- c(0, 0.5, 0.6, 0.7, 0.8, r_tmp, 0.9, 0.95, 0.99)
mu_hat_0 <- spline(r_set, mu_hat_set, xout = r_0)$y
sigma2_hat_0 <- spline(r_set, sigma2_hat_set, xout = r_0)$y
tau2_hat_0 <- spline(r_set, tau2_hat_set, xout = r_0)$y

my_table5 <- rbind(
  cbind(r_0, mu_hat_0, sqrt(sigma2_hat_0), sqrt(tau2_hat_0)),
  c(lim_r, lim_mu, sqrt(lim_sigma2), sqrt(lim_tau2))
)

colnames(my_table5) <- c("rho", "b_0", "sigma_0", "tau")
my_table5 <- as.data.frame(my_table5)

my_table5_rounded <- round(my_table5, 2)

# output produced by RStudio 2025.05 on MacBook

# Table 5: IPCC5
#       0    1.07    0.27    0.04
#    0.50    1.20    0.32    0.06
#    0.60    1.21    0.35    0.06
#    0.70    1.21    0.40    0.06
#    0.80    1.22    0.49    0.06
#    0.83    1.22    0.53    0.06
#    0.90    1.22    0.68    0.06
#    0.95    1.23    0.96    0.06
#    0.99    1.23    2.13    0.06
#    1.00    1.23     Inf    0.06

# Table 5: IPCC6
#        0    1.11    0.10    0.04
#     0.50    1.10    0.13    0.09
#     0.60    1.09    0.15    0.10
#     0.70    1.09    0.17    0.12
#     0.80    1.07    0.21    0.14
#     0.88    1.06    0.27    0.16
#     0.90    1.05    0.30    0.17
#     0.95    1.03    0.42    0.19
#     0.99    1.01    0.92    0.21
#     1.00    1.00     Inf    0.22


# This R code replicates Excel calculations from IPCC_prior_calculations_v7_revision_github.xlsx


# Define prior parameters based on case
if (caseN == 14) {
  b_2 <- 1.07
  sigma_2 <- 0.53
} else if (caseN == 4) {
  b_2 <- 1.15
  sigma_2 <- 0.27
} else {
  stop("not implemented yet")
}

# Indices for rho = 0.9, 0.95, 0.99
idx <- which(my_table5$rho %in% c(0.9, 0.95, 0.99))

# Extract sigma_0 and b_0 for selected rhos
sigma2_0_sel <- my_table5$sigma_0[idx]^2
b_0_sel <- my_table5$b_0[idx]

# Compute combined posterior
sigma_1 <- sqrt((sigma2_0_sel * sigma_2^2) / (sigma2_0_sel - sigma_2^2))
b_1 <- (sigma2_0_sel * b_2 - sigma_2^2 * b_0_sel) / (sigma2_0_sel - sigma_2^2)

# Construct Table 6
my_table6 <- rbind(
  cbind(my_table5$rho[idx], b_1, sigma_1),
  c(1.0, b_2, sigma_2)
)

colnames(my_table6) <- c("rho", "b_1", "sigma_1")
my_table6 <- as.data.frame(my_table6)

# Rounded version
my_table6_rounded <- round(my_table6, 2)

# output produced by RStudio 2025.05 on MacBook

# Table 6: IPCC5
#  rho  b_1 sigma_1
# 0.90 0.83    0.85
# 0.95 1.00    0.64
# 0.99 1.06    0.55
# 1.00 1.07    0.53

# Table 6: IPCC6
#  rho  b_1 sigma_1
# 0.90 1.66    0.67
# 0.95 1.24    0.36
# 0.99 1.16    0.28
# 1.00 1.15    0.27