library(RcppEigen)
source('New Code/IWTlm/IWTlmpointwise_FoF.R')

##load test model
load('New Code/IWTlm/test_data_IWTlm.rdata')

attach(cov)

test_model=fastIWTlm(co2_mat_log ~ population_mat_log + gdppc_mat_log + end + ff + lc)
