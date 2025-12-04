## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----include = FALSE----------------------------------------------------------
ggplot2dplyr_available <- requireNamespace("ggplot2", quietly = TRUE) && 
  requireNamespace("dplyr", quietly = TRUE)

## -----------------------------------------------------------------------------
set.seed(1)

## -----------------------------------------------------------------------------
library(disclapmix)
data(danes)
db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])
str(db)

## -----------------------------------------------------------------------------
default_fits <- disclapmix_adaptive(db, label = "PAM", margin = 5L)

## -----------------------------------------------------------------------------
clara_fits <- disclapmix_adaptive(db, label = "CLARA", margin = 5L, init_y_method = "clara")

## -----------------------------------------------------------------------------
# Random observations:
my_init_y_generator <- function(k) {
  # Or cluster::pam(), cluster::clara() or something else
  db[sample(seq_len(nrow(db)), k, replace = FALSE), , drop = FALSE]
}

my_init_y_generator(1)
my_init_y_generator(2)

## -----------------------------------------------------------------------------
custom_fits <- disclapmix_adaptive(db, label = "Custom", 
                                   margin = 1L, # Just demonstrating my_init_y_generator()
                                   init_y_generator = my_init_y_generator)
rm(custom_fits_best) # To avoid using it by accident later

## -----------------------------------------------------------------------------
set.seed(2) # For reproducibility
custom_fits_extra <- replicate(5, 
                               disclapmix_adaptive(db, 
                                                   label = "Custom", 
                                                   margin = 5L, 
                                                   init_y_generator = my_init_y_generator, 
                                                   # Random starting points may need more iterations
                                                   glm_control_maxit = 100L)
)
str(custom_fits_extra, 2)

custom_fits_max_n <- max(sapply(custom_fits_extra, length))
custom_fits_best <- vector("list", custom_fits_max_n)

for (i in seq_len(custom_fits_max_n)) {
  best_fit_i <- NULL
  
  for (j in seq_along(custom_fits_extra)) {
    if (length(custom_fits_extra[[j]]) < i) {
      next
    }
    
    if (is.null(best_fit_i) || 
        best_fit_i$BIC_marginal > custom_fits_extra[[j]][[i]]$BIC_marginal) {
      
      best_fit_i <- custom_fits_extra[[j]][[i]]
    }
  }
  
  custom_fits_best[[i]] <- best_fit_i
}

## -----------------------------------------------------------------------------
fits <- c(default_fits, clara_fits, custom_fits_best)

## -----------------------------------------------------------------------------
d <- data.frame(
  Label = sapply(fits, function(x) x$label),
  BIC = sapply(fits, function(x) x$BIC_marginal),
  Clusters = sapply(fits, function(x) nrow(x$y))
)

## ----eval = ggplot2dplyr_available, fig.width=8, fig.height=5-----------------
library(ggplot2)
ggplot(d, aes(Clusters, BIC, color = Label)) + 
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = unique(d$Clusters)) + 
  theme_bw()

## ----eval = ggplot2dplyr_available--------------------------------------------
library(dplyr)

d %>% 
  group_by(Label) %>% 
  summarise(best_clusters = Clusters[which.min(BIC)])

## ----eval = FALSE-------------------------------------------------------------
# saveRDS(default_fits, "obj-default_fits.Rdata")
# saveRDS(clara_fits, "obj-clara_fits.Rdata")
# saveRDS(custom_fits_best, "obj-custom_fits_best.Rdata")

## -----------------------------------------------------------------------------
fits <- disclapmix_adaptive(db, criteria = c("BIC", "AIC", "AICc"), margin = 5L)
length(fits)
d <- data.frame(
  BIC = sapply(fits, function(x) x$BIC_marginal),
  AIC = sapply(fits, function(x) x$AIC_marginal),
  AICc = sapply(fits, function(x) x$AICc_marginal),
  Clusters = sapply(fits, function(x) nrow(x$y))
)
best_BIC <- d$Clusters[which.min(d$BIC)]
best_AIC <- d$Clusters[which.min(d$AIC)]
best_AICc <- d$Clusters[which.min(d$AICc)]

## ----eval = ggplot2dplyr_available, fig.width=8, fig.height=5-----------------
library(ggplot2)
ggplot(d) + 
  
  geom_vline(aes(xintercept = best_BIC, color = "BIC"), linetype = "dashed") +
  geom_vline(aes(xintercept = best_AIC, color = "AIC"), linetype = "dashed") +
  geom_vline(aes(xintercept = best_AICc, color = "AICc"), linetype = "dashed") +
  
  geom_line(aes(Clusters, BIC, color = "BIC")) +
  geom_line(aes(Clusters, AIC, color = "AIC")) +
  geom_line(aes(Clusters, AICc, color = "AICc")) +
  
  scale_x_continuous(breaks = unique(d$Clusters)) + 
  
  labs(y = "Information criteria value", color = "Information criteria") + 
  
  theme_bw()

