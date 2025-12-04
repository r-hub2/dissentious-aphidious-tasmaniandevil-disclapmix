## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(fig.width = 5, fig.height = 5)

## -----------------------------------------------------------------------------
library(ggplot2)

library(disclapmix)
data(danes)

## -----------------------------------------------------------------------------
db <- as.matrix(danes[rep(seq_len(nrow(danes)), danes$n), seq_len(ncol(danes)-1)])
str(db)

## -----------------------------------------------------------------------------
fits <- disclapmix_adaptive(x = db)

## -----------------------------------------------------------------------------
BICs <- sapply(fits, function(x) x$BIC_marginal)
bestfit <- fits[[which.min(BICs)]]
bestfit

## -----------------------------------------------------------------------------
clusters <- seq_along(fits) # or: sapply(fits, function(x) nrow(x$y))
plot(clusters, BICs, type = "b")
points(nrow(bestfit$y), bestfit$BIC_marginal, pch = 4, cex = 4, col = "red")

## -----------------------------------------------------------------------------
db_haplotypes <- as.matrix(danes[, seq_len(ncol(db))])
p <- predict(bestfit, db_haplotypes)
plot(log10(danes$n / sum(danes$n)), log10(p), 
     xlab = "Observed relative frequency", 
     ylab = "Discrete Laplace estimated probability")
abline(0, 1)

## ----eval=FALSE---------------------------------------------------------------
# fit <- disclapmix(x = db, clusters = 2L)

## -----------------------------------------------------------------------------
clusters <- 1L:5L
fits <- lapply(clusters, function(clusters) {
  fit <- disclapmix(x = db, clusters = clusters)
  return(fit)
})

marginalBICs <- sapply(fits, function(fit) {
  fit$BIC_marginal
})

bestfit <- fits[[which.min(marginalBICs)]]

## -----------------------------------------------------------------------------
bestfit
summary(bestfit)

## ----fig.width=10-------------------------------------------------------------
plot(bestfit)

## -----------------------------------------------------------------------------
bestfit$y
bestfit$disclap_parameters

## -----------------------------------------------------------------------------
disclap_estimates <- predict(bestfit, 
                             newdata = as.matrix(danes[, 1:(ncol(danes) - 1)]))

## -----------------------------------------------------------------------------
ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(x = danes$n/sum(danes$n), y = disclap_estimates)) +
  labs(x = "Observed frequency",
       y = "Predicted frequency (discrete Laplace)") +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()

