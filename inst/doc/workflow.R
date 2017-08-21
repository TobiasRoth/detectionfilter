## ------------------------------------------------------------------------
# Proportion of  species observed at least once but in less than 10% of sites
mean(apply(commat_obs, 2, sum) < 20)

