library(purrr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(S4Vectors)
library(assertthat)
library(stats)
library(mclust)
library(MASS)
library(mvtnorm)

Mode <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  ux <- unique(x)
  freq <- tabulate(match(x, ux))
  ux[which.max(freq)]
}

log_dmvnorm <- function(x, mean, sigma) {
  d <- length(x)
  # 尝试 Cholesky 分解
  L <- tryCatch(chol(sigma), error = function(e) NA)
  if (any(is.na(L))) return(-Inf)
  log_det <- 2 * sum(log(diag(L)))
  diff <- x - mean
  # 解方程 L y = diff，利用 backsolve（L 为上三角矩阵）
  y <- backsolve(L, diff, upper.tri = TRUE)
  quadform <- sum(y^2)
  -0.5 * d * log(2 * pi) - 0.5 * log_det - 0.5 * quadform
}

iterate_t <- function(Y, df_j, nrep, thin, n, d, gamma, q, init,
                      mu0 = colMeans(Y), lambda0 = diag(0.01, ncol(Y)),
                      alpha = 1, beta = 0.01) {
  num_samples <- floor(nrep / thin) + 1
  df_sim_z <- matrix(0, nrow = num_samples, ncol = n)
  df_sim_mu <- matrix(0, nrow = num_samples, ncol = q * d)
  df_sim_lambda <- vector("list", num_samples)
  df_sim_w <- matrix(0, nrow = num_samples, ncol = n)

  # 初始化
  initmu <- rep(mu0, q)
  df_sim_mu[1, ] <- initmu
  lambda_i <- lambda0
  df_sim_lambda[[1]] <- lambda0
  z <- init  
  df_sim_z[1, ] <- z
  w <- rep(1, n)
  df_sim_w[1, ] <- w
  
  mu0vec <- mu0
  
  for (i in 1:(nrep - 1)) {
    if(i %% 100 ==0){cat(i,"\n")}
    
    # Update mu  
    mu_i <- matrix(0, nrow = q, ncol = d)
    for (k in 1:q) {
      idx <- which(z == k)
      if (length(idx) == 0) {
        n_i <- 0
        Ysums <- rep(0, d)
      } else {
        n_i <- sum(w[idx])
        Ysums <- colSums(Y[idx, , drop = FALSE] * w[idx])
      }
      post_prec <- lambda0 + n_i * lambda_i
      post_cov <- solve(post_prec)
      term <- lambda0 %*% mu0vec + lambda_i %*% Ysums
      mean_i <- as.vector(post_cov %*% term)
      # 从后验分布中采样
      mu_i[k, ] <- MASS::mvrnorm(1, mu = mean_i, Sigma = post_cov)
    }
    
    # Update Lambda
    mu_i_long <- mu_i[z, , drop = FALSE]
    diff_mat <- Y - mu_i_long
    # 计算 Σ = t(diff) %*% diag(w) %*% diff
    sumofsq <- t(diff_mat) %*% (diff_mat * w)
    Vinv <- diag(rep(beta, d))
    S_mat <- solve(Vinv + sumofsq)
    lambda_i <- rWishart(1, n + alpha, S_mat)[,,1]
    sigma_i <- solve(lambda_i)
    
    # Update w  
    w_alpha <- (d + 4) / 2
    mu_i_long <- mu_i[z, , drop = FALSE]
    diff_mat <- Y - mu_i_long
    quad <- rowSums((diff_mat %*% lambda_i) * diff_mat)
    w_beta <- 2 / (quad + 4)
    w <- rgamma(n, shape = w_alpha, scale = w_beta)
    
    # Update z
    for (j in 1:n) {
      z_j_prev <- z[j]
      candidate_clusters <- setdiff(1:q, z_j_prev)
      if (length(candidate_clusters) == 0) next
      z_j_new <- sample(candidate_clusters, 1)
      nbrs <- df_j[[j]]
      if (length(nbrs) > 0) {
        neighbor_factor <- gamma * 2 / length(nbrs)
        h_z_prev <- neighbor_factor * sum(z[nbrs] == z_j_prev) +
          log_dmvnorm(Y[j, ], mean = mu_i[z_j_prev, ], sigma = sigma_i / w[j])
        h_z_new <- neighbor_factor * sum(z[nbrs] == z_j_new) +
          log_dmvnorm(Y[j, ], mean = mu_i[z_j_new, ], sigma = sigma_i / w[j])
      } else {
        h_z_prev <- log_dmvnorm(Y[j, ], mean = mu_i[z_j_prev, ], sigma = sigma_i / w[j])
        h_z_new <- log_dmvnorm(Y[j, ], mean = mu_i[z_j_new, ], sigma = sigma_i / w[j])
      }
      prob_j <- min(exp(h_z_new - h_z_prev), 1)
      z[j] <- sample(c(z_j_prev, z_j_new), 1, prob = c(1 - prob_j, prob_j))
    }

    if ((i + 1) %% thin == 0) {
      sample_index <- (i + 1) / thin + 1  # 第一行为初始状态
      df_sim_mu[sample_index, ] <- as.vector(t(mu_i))  # 行优先展开
      df_sim_lambda[[sample_index]] <- lambda_i
      df_sim_w[sample_index, ] <- w
      df_sim_z[sample_index, ] <- z
    }
  }
  
  list(z = df_sim_z, mu = df_sim_mu, lambda = df_sim_lambda,
       weights = df_sim_w)
}

cluster <- function(
    Y, q, df_j, init = rep(1, nrow(Y)),
    mu0 = colMeans(Y), lambda0 = diag(0.01, ncol(Y)),
    gamma = 3, alpha = 1, beta = 0.01, nrep = 1000, thin = 100
) {
  Y <- as.matrix(Y)
  d <- ncol(Y)
  n <- nrow(Y)
  
  if (q == 1) {
    return(list(z = matrix(rep(1, n), nrow = 1)))
  }
  
  message("Fitting model...")
  
  iterate_t(Y = Y, df_j = df_j, nrep = nrep, thin = thin, n = n, d = d,
    gamma = gamma, q = q, init = init, mu0 = mu0, lambda0 = lambda0, alpha = alpha, beta = beta)
}

.find_neighbors <- function(sce, platform) {
  if (platform == "Visium") {
    offsets <- data.frame(
      x.offset = c(-2, 2, -1, 1, -1, 1),
      y.offset = c(0, 0, -1, -1, 1, 1)
    )
  } else if (platform %in% c("VisiumHD", "ST")) {
    offsets <- data.frame(
      x.offset = c(0, 1, 0, -1),
      y.offset = c(-1, 0, 1, 0)
    )
  } else {
    stop(".find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  sce$spot.idx <- seq_len(ncol(sce))
  spot.positions <- colData(sce)[, c("spot.idx", "array_col", "array_row")]
  
  neighbor.positions <- merge(spot.positions, offsets)
  neighbor.positions$x.pos <- neighbor.positions$array_col + neighbor.positions$x.offset
  neighbor.positions$y.pos <- neighbor.positions$array_row + neighbor.positions$y.offset
  
  neighbors <- merge(as.data.frame(neighbor.positions),
                     as.data.frame(spot.positions),
                     by.x = c("x.pos", "y.pos"), by.y = c("array_col", "array_row"),
                     suffixes = c(".primary", ".neighbor"),
                     all.x = TRUE)
  
  neighbors <- neighbors[order(neighbors$spot.idx.primary,
                               neighbors$spot.idx.neighbor), ]
  df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
  df_j <- purrr::map(df_j, function(nbrs) discard(nbrs, is.na))
  
  sce$spot.neighbors <- vapply(
    df_j,
    function(x) {
      if (length(x) == 0) {
        return(NA_character_)
      } else {
        return(paste0(x, collapse = ","))
      }
    },
    FUN.VALUE = character(1)
  )
  
  n_with_neighbors <- length(keep(df_j, function(nbrs) length(nbrs) > 0))
  message("Neighbors were identified for ", n_with_neighbors, " out of ", ncol(sce), " spots.")
  list(sce, unname(df_j))
}

.init_cluster <- function(Y, q, init = NULL, init.method = c("mclust", "kmeans")) {
  if (is.null(init)) {
    init.method <- match.arg(init.method)
    if (init.method == "kmeans") {
      init <- kmeans(Y, centers = q)$cluster
    } else if (init.method == "mclust") {
      init <- Mclust(Y, q, "EEE", verbose = FALSE)$classification
    }
  }
  init
}

spatialCluster <- function(sce, q, use.dimred = "PCA", d = 15,
                           platform = c("Visium", "VisiumHD", "ST"),
                           init = NULL, init.method = c("mclust", "kmeans"),
                           nrep = 50000, burn.in = 1000, thin = 100, gamma = NULL,
                           mu0 = NULL, lambda0 = NULL, alpha = 1, beta = 0.01) {
  assert_that(nrep >= 100)
  assert_that(burn.in >= 0)
  if (burn.in >= nrep) {
    stop("Please specify a burn-in period shorter than the total number of iterations.")
  }
  Y <- reducedDim(sce, use.dimred)
  d <- min(ncol(Y), d)
  Y <- Y[, seq_len(d)]
  if (length(platform) > 1) {
    platform <- .bsData(sce, "platform", match.arg(platform))
  } else {
    platform <- match.arg(platform)
  }
  .neighbors <- .find_neighbors(sce, platform)
  sce <- .neighbors[[1]]
  df_j <- .neighbors[[2]]  # 注意：这里 df_j 中邻居索引为 1-indexed
  init <- .init_cluster(Y, q, init, init.method)
  if (is.null(init)) {
    stop("Empty initialization. Please use a different initialization method.")
  }
  if (is.null(mu0)) {
    mu0 <- colMeans(Y)
  }
  if (is.null(lambda0)) {
    lambda0 <- diag(0.01, ncol(Y))
  }
  if (is.null(gamma)) {
    if (platform == "Visium") {
      gamma <- 3
    } else if (platform %in% c("VisiumHD", "ST")) {
      gamma <- 2
    }
  }
  results <- cluster(Y, q, df_j, init = init, mu0 = mu0, lambda0 = lambda0,
                     gamma = gamma, alpha = alpha, beta = beta, nrep = nrep,
                     thin = thin)
  sce$cluster.init <- init
  if (!exists("BayesSpace.data", metadata(sce))) {
    metadata(sce)$BayesSpace.data <- list()
  }
  metadata(sce)$BayesSpace.data$platform <- platform
  message("Calculating labels using iterations ", burn.in + 1, " through ", nrep, ".")
  .burn.in <- burn.in %/% thin
  .nrep <- nrep %/% thin
  zs <- results$z[seq(.burn.in + 2, .nrep + 1), ]
  if (.burn.in + 1 == .nrep) {
    labels <- matrix(zs, nrow = 1)
  } else {
    labels <- apply(zs, 2, Mode)
  }
  colData(sce)$spatial.cluster <- unname(labels)
  sce
}