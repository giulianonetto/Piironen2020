# Define auxiliary functions ---------------------------------------------------
generate_data_01 <- function(N = 50, p = 500, rho = 0.3, p_rel = p/2) {
  f <- rnorm(n = N, mean = 0, sd = 1)
  y <- rnorm(n = N, mean = f, sd = 1)
  X <- matrix(nrow = N, ncol = p)
  for (i in 1:p) {
    if (i > p_rel) {
      X[, i] <- rnorm(n = N, mean = 0, sd = 1)
    } else {
      X[, i] <- rnorm(n = N, mean = sqrt(rho)*f, sd = sqrt(1 - rho))
    }
  }
  
  return(data.frame(y, X))
}

rmse <- function(pred, truth) {
  sqrt(mean((truth - pred)^2))
}

generate_data_02 <- function(N = 50, p = 500, rho = 0.3, p_rel = p/2, ...) {
  alist <- list(...)
  if (!("beta" %in% names(alist)) ) {
    print(str_glue("Simulating beta"))
    beta <- c(1, runif(p_rel - 1, 50, 100), rep(0, p - p_rel))
  } else {
    beta <- alist[['beta']]
  }
  #X <- cbind(1, replicate(p - 1, rnorm(N, 0, 1)))
  Sigma <- diag(rep(1, p - 1))
  Sigma[which(lower.tri(diag(1:(p_rel-1))), arr.ind = T)] <- rho
  Sigma[which(upper.tri(diag(1:(p_rel-1))), arr.ind = T)] <- rho
  X <- cbind(
    1,
    MASS::mvrnorm(n = N, mu = rep(0, p-1), Sigma = Sigma)
  )
  y <- rnorm(n = N, mean = X %*% beta, sd = 1)
  return(data.frame(y, X[,-1]))
}

fit_projpred <- function(df) {
  n <- nrow(df)
  D <- ncol(df[, -1])
  p0 <- floor(D/10) # prior guess for the number of relevant variables
  tau0 <- p0/(D-p0) * 1/sqrt(n) # scale for tau
  fit <- brm(y ~ ., family=gaussian(), data=df,
             prior=prior(horseshoe(scale_global = tau0, scale_slab = 1), class=b),
             seed=1, refresh = 0, cores = 4)
  fit <- rstanarm::stan_lm(y ~ ., data = df, cores = 4)
  refmodel <- get_refmodel(fit)
  vs <- varsel(refmodel, method = "l1", nterms_max = 50)
  .size <- NA; .pct <- 0.05
  while(is.na(.size)) {
    .size <- suggest_size(vs, pct = .pct, stat='rmse')
    .pct <- .pct * 2
  }
  # pred <- proj_linpred(vs, nterms = .size, integrated = TRUE)
  return(list(vs = vs, .size = .size))
  
}

fit_homemade_proj <- function(df, fit_l1, y_ref) {
  n <- nrow(df)
  d <- ncol(df[, -1])
  
  # use glmnet to search potential submodels
  coef_l1 <- as.matrix(coef(fit_l1,s = 'lambda.1se'))[-1, ]
  coef_l1 <- coef_l1[abs(coef_l1)>0] %>% sort(decreasing = T)
  if (length(coef_l1) > 100) {
    coef_l1 <- coef_l1[1:100]
  }
  pred_l1 <- predict(fit_l1, newx = as.matrix(df[,-1]), s = 'lambda.1se')
  nsubmodels <- length(coef_l1)
  rmse_sub_models <- vector('numeric', length = nsubmodels)
  best_ix <- NA
  for (i in nsubmodels:1) {
    ix <- names(coef_l1)[1:i]
    # TODO: Error: Selections can't have missing values. below
    .subdf <- df %>% select(all_of(ix))
    
    subfit <- lm(y_ref ~ ., data = data.frame(y_ref, .subdf))
    rmse_submodel <- sqrt(mean(
      (df$y - predict(subfit))^2
    ))
    rmse_sub_models[i] <- rmse_submodel
    r <- rmse_submodel/min(rmse_sub_models[rmse_sub_models>0]) 
    if (r <= 5) {
      #print(str_glue('updating to i={i}'))
      best_ix <- ix
    }
  }
  
  .df <- df %>% select(all_of(best_ix))
  fit <- lm(y ~ . , data = data.frame(y = y_ref, .df))
  
  return(fit)
  
  
}
