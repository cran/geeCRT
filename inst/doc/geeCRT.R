## ----start,echo=FALSE,results="hide"------------------------------------------

library(geeCRT)


## ---- echo=FALSE, results='asis'----------------------------------------------
knitr::kable(head(sampleSWCRTSmall))

## ----createz, fig.keep="all", fig.width = 7, fig.height=4---------------------

createzCrossSec = function (m) {

    Z = NULL
    n = dim(m)[1]
    
    for (i in 1:n) {
        
        alpha_0 = 1; alpha_1 = 2; n_i = c(m[i, ]); n_length = length(n_i)
        POS = matrix(alpha_1, sum(n_i), sum(n_i))
        loc1 = 0; loc2 = 0
        
        for (s in 1:n_length) {
            
            n_t = n_i[s]; loc1 = loc2 + 1; loc2 = loc1 + n_t - 1
            
            for (k in loc1:loc2) {

                for (j in loc1:loc2) {

                    if (k != j) { POS[k, j] = alpha_0 } else { POS[k, j] = 0 }}}}

        zrow = diag(2); z_c = NULL
        
        for (j in 1:(sum(n_i) - 1)) { 

            for (k in (j + 1):sum(n_i)) {z_c = rbind(z_c, zrow[POS[j,k],])}}
        
        Z = rbind(Z, z_c) }

    return(Z)}

## ----geemaee_small, fig.keep="all", fig.width = 7, fig.height=4---------------

sampleSWCRT = sampleSWCRTSmall

### Individual-level id, period, outcome, and design matrix
id = sampleSWCRT$id; period = sampleSWCRT$period;
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'treatment')])
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]

### design matrix for correlation parameters
Z = createzCrossSec(m) 

### (1) Matrix-adjusted estimating equations and GEE 
### on continous outcome with nested exchangeable correlation structure
 
### MAEE
est_maee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                           X = X, id  = id, Z = Z, 
                           family = "continuous", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)

### GEE
est_uee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                          X = X, id = id, Z = Z, 
                          family = "continuous", 
                          maxiter = 500, epsilon = 0.001, 
                          printrange = TRUE, alpadj = FALSE, 
                          shrink = "ALPHA", makevone = FALSE)

### (2) Matrix-adjusted estimating equations and GEE 
### on binary outcome with nested exchangeable correlation structure

### MAEE
est_maee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                           X = X, id = id, Z = Z, 
                           family = "binomial", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_maee_ind_bin)

### GEE
est_uee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                          X = X, id = id, Z = Z, 
                          family = "binomial", 
                          maxiter = 500, epsilon = 0.001, 
                          printrange = TRUE, alpadj = FALSE, 
                          shrink = "ALPHA", makevone = FALSE)


## ----set-options1, echo=FALSE, fig.keep="all", fig.width = 7, fig.height=4------------------------
options(width = 100)
 

## ---- fig.keep="all", fig.width = 7, fig.height=4-------------------------------------------------
 # MAEE for continuous outcome
 print(est_maee_ind_con)

 # GEE for continuous outcome
 print(est_uee_ind_con)
 
 # MAEE for binary outcome
 print(est_maee_ind_bin)
 
 # GEE for binary outcome
 print(est_uee_ind_bin)


## ---- echo=FALSE, results='asis'------------------------------------------------------------------
knitr::kable(head(sampleSWCRTLarge))

## ----geemaee_large, fig.keep="all", fig.width = 7, fig.height=4-----------------------------------

########################################################################
### Example 2): simulated SW-CRT with larger cluster-period sizes (20~30)
########################################################################
## This will elapse longer. 
sampleSWCRT = sampleSWCRTLarge

### Individual-level id, period, outcome, and design matrix
id = sampleSWCRT$id; period =  sampleSWCRT$period;
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'period5', 'treatment')])
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
### design matrix for correlation parameters
Z = createzCrossSec(m) 

### (1) Matrix-adjusted estimating equations and GEE 
### on continous outcome with nested exchangeable correlation structure
 
### MAEE
est_maee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                           X = X, id  = id, Z = Z, 
                           family = "continuous", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_maee_ind_con)

### GEE
est_uee_ind_con = geemaee(y = sampleSWCRT$y_con, 
                          X = X, id = id, Z = Z, 
                          family = "continuous", 
                          maxiter = 500, epsilon = 0.001, 
                          printrange = TRUE, alpadj = FALSE, 
                          shrink = "ALPHA", makevone = FALSE)
print(est_uee_ind_con)

### (2) Matrix-adjusted estimating equations and GEE 
### on binary outcome with nested exchangeable correlation structure

### MAEE
est_maee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                           X = X, id = id, Z = Z, 
                           family = "binomial", 
                           maxiter = 500, epsilon = 0.001, 
                           printrange = TRUE, alpadj = TRUE, 
                           shrink = "ALPHA", makevone = FALSE)
print(est_maee_ind_bin)

### GEE
est_uee_ind_bin = geemaee(y = sampleSWCRT$y_bin, 
                          X = X, id = id, Z = Z, 
                          family = "binomial", 
                          maxiter = 500, epsilon = 0.001, 
                          printrange = TRUE, alpadj = FALSE, 
                          shrink = "ALPHA", makevone = FALSE)
print(est_uee_ind_bin)



## ----cpgee_gather_small, fig.keep="all", fig.width = 7, fig.height=4------------------------------

########################################################################
### Example 1): simulated SW-CRT with smaller cluster-period sizes (5~10)
########################################################################

sampleSWCRT = sampleSWCRTSmall

### cluster-period id, period, outcome, and design matrix
### id, period, outcome
id = sampleSWCRT$id; period = sampleSWCRT$period; y =  sampleSWCRT$y_bin
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'treatment')])
 
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
clp_mu = tapply(y,list(id,period), FUN=mean)
y_cp = c(t(clp_mu))
 
### design matrix for correlation parameters
trt = tapply(X[, t + 1], list(id, period), FUN=mean)
trt = c(t(trt))

time = tapply(period,list(id, period), FUN = mean); time = c(t(time))
X_cp = matrix(0, n * t, t)

s = 1
for (i in 1:n) { for (j in 1:t) { X_cp[s, time[s]] = 1; s = s + 1 }}
X_cp = cbind(X_cp, trt); id_cp = rep(1:n, each = t); m_cp =  c(t(m))

## ----cpgee_small, fig.keep="all", fig.width = 7, fig.height=4-------------------------------------

### cluster-period matrix-adjusted estimating equations (MAEE) 
### with exchangeable, nested exchangeable and exponential decay correlation structures 
# exponential
est_maee_exc = cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exchangeable", 
                        alpadj = TRUE)
print(est_maee_exc)

# nested exchangeable
est_maee_nex = cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "nest_exch", 
                        alpadj = TRUE)
print(est_maee_nex)

# exponential decay 
est_maee_ed = cpgeeSWD(y  = y_cp, X = X_cp, id = id_cp, 
                       m = m_cp, corstr = "exp_decay", 
                       alpadj = TRUE)
print(est_maee_ed)
 

### cluster-period GEE 
### with exchangeable, nested exchangeable and exponential decay correlation structures

# exchangeable
est_uee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exchangeable",
                        alpadj = FALSE)
print(est_uee_exc)

# nested exchangeable
est_uee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "nest_exch", 
                        alpadj = FALSE)
print(est_uee_nex)

# exponential decay 
est_uee_ed <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                       m = m_cp, corstr = 'exp_decay', 
                       alpadj = FALSE)
print(est_uee_ed)




## ----cpgee_large, fig.keep="all", fig.width = 7, fig.height=4-------------------------------------

########################################################################
### Example 2): simulated SW-CRT with larger cluster-period sizes (20~30)
########################################################################

sampleSWCRT = sampleSWCRTLarge

### cluster-period id, period, outcome, and design matrix
### id, period, outcome
id = sampleSWCRT$id; period =  sampleSWCRT$period; y =  sampleSWCRT$y_bin
X = as.matrix(sampleSWCRT[, c('period1', 'period2', 'period3', 'period4', 'period5', 'treatment')])
 
m = as.matrix(table(id, period)); n = dim(m)[1]; t = dim(m)[2]
clp_mu<-tapply(y,list(id,period), FUN=mean)
y_cp <- c(t(clp_mu))
 
### design matrix for correlation parameters
trt <- tapply(X[, t + 1], list(id, period), FUN=mean)
trt <- c(t(trt))

time <- tapply(period,list(id, period), FUN = mean); time <- c(t(time))
X_cp <- matrix(0, n * t, t)

s = 1
for(i in 1:n){for(j in 1:t){X_cp[s, time[s]] <- 1; s = s + 1}}
X_cp <- cbind(X_cp, trt); id_cp <- rep(1:n, each= t); m_cp <-  c(t(m))

### cluster-period matrix-adjusted estimating equations (MAEE) 
### with exchangeable, nested exchangeable and exponential decay correlation structures 
# exponential
est_maee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                         m = m_cp, corstr = "exchangeable", 
                         alpadj = TRUE)
print(est_maee_exc)

# nested exchangeable
est_maee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                         m = m_cp, corstr = "nest_exch", 
                         alpadj = TRUE)
print(est_maee_nex)

# exponential decay 
est_maee_ed <- cpgeeSWD(y  = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exp_decay", 
                        alpadj = TRUE)
print(est_maee_ed)
 

### cluster-period GEE 
### with exchangeable, nested exchangeable and exponential decay correlation structures

# exchangeable
est_uee_exc <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "exchangeable",
                        alpadj = FALSE)
print(est_uee_exc)

# nested exchangeable
est_uee_nex <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                        m = m_cp, corstr = "nest_exch", 
                        alpadj = FALSE)
print(est_uee_nex)

# exponential decay 
est_uee_ed <- cpgeeSWD(y = y_cp, X = X_cp, id = id_cp, 
                       m = m_cp, corstr = 'exp_decay', 
                       alpadj = FALSE)
print(est_uee_ed)


## ----sim_ep, fig.keep="all", fig.width = 7, fig.height=4------------------------------------------


#### Simulate 2 clusters, 3 periods and cluster-period size of 5

t = 3; n = 2; m = 5
# means of cluster 1
u_c1 = c(0.4, 0.3, 0.2)
u1 <- rep(u_c1, c(rep(m, t)))
# means of cluster 2
u_c2 = c(0.35, 0.25, 0.2)
u2 <- rep(u_c2, c(rep(m, t)))

# List of mean vectors
mu = list(); mu[[1]] = u1; mu[[2]] = u2;

# List of correlation matrices

## correlation parameters
alpha0 = 0.03; alpha1 = 0.015; rho = 0.8

## (1) exchangeable
Sigma = list()
Sigma[[1]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )

Sigma[[2]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )

y_exc = simbinPROBIT(mu = mu, Sigma = Sigma, n = n)

## (2) nested exchangeable
Sigma = list()
cor_matrix = matrix(alpha1, m * t,  m * t)
loc1 = 0; loc2 = 0
for(t in 1:t){loc1 = loc2 + 1; loc2 = loc1 + m - 1
  for(i in loc1:loc2){for(j in loc1:loc2){
         if(i != j){cor_matrix[i, j] = alpha0}else{cor_matrix[i, j] = 1}}}}
 
Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix

y_nex = simbinPROBIT(mu = mu, Sigma = Sigma, n = n)

## (3) exponential decay
 
Sigma = list()
 
### function to find the period of the ith index
 region_ij<-function(points, i){diff = i - points
     for(h in 1:(length(diff) - 1)){if(diff[h] > 0 & diff[h + 1] <= 0){find <- h}}
  return(find)}

 cor_matrix = matrix(0,  m * t,  m * t)
 useage_m = cumsum(m * t); useage_m = c(0, useage_m)

 for(i in 1:(m * t)){i_reg = region_ij(useage_m, i)
      for(j in 1:(m * t)){j_reg = region_ij(useage_m, j)
          if(i_reg == j_reg & i != j){
              cor_matrix[i, j] = alpha0}else if(i == j){cor_matrix[i, j] = 1
 }else if(i_reg != j_reg){cor_matrix[i,j] = alpha0 * (rho^(abs(i_reg - j_reg)))}}}

Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix

y_ed = simbinPROBIT(mu = mu, Sigma = Sigma, n = n)


## ----sim_clf, fig.keep="all", fig.width = 7, fig.height=4-----------------------------------------
##### Simulate 2 clusters, 3 periods and cluster-period size of 5

t = 3; n = 2; m = 5
 
# means of cluster 1
u_c1 = c(0.4, 0.3, 0.2)
u1 <- rep(u_c1, c(rep(m, t)))
# means of cluster 2
u_c2 = c(0.35, 0.25, 0.2)
u2 <- rep(u_c2, c(rep(m, t)))

# List of mean vectors
mu = list()
mu[[1]] = u1; mu[[2]] = u2;

# List of correlation matrices
 
## correlation parameters
alpha0 = 0.03; alpha1 = 0.015; rho = 0.8

## (1) exchangeable
Sigma = list()
Sigma[[1]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )
Sigma[[2]] = diag(m * t) * ( 1 - alpha1) + matrix(alpha1, m * t,  m * t )
y_exc = simbinCLF(mu = mu, Sigma = Sigma, n = n)

## (2) nested exchangeable
Sigma = list()
cor_matrix = matrix(alpha1, m * t,  m * t)
loc1 = 0; loc2 = 0
for(t in 1:t){loc1 = loc2 + 1; loc2 = loc1 + m - 1
    for(i in loc1:loc2){for(j in loc1:loc2){
         if(i != j){cor_matrix[i, j] = alpha0}else{cor_matrix[i, j] = 1}}}}
 
Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix
y_nex = simbinCLF(mu = mu, Sigma = Sigma, n = n)

## (3) exponential decay
 
Sigma = list()

### function to find the period of the ith index
region_ij<-function(points, i){diff = i - points
     for(h in 1:(length(diff) - 1)){if(diff[h] > 0 & diff[h + 1] <= 0){find <- h}}
  return(find)}

cor_matrix = matrix(0,  m * t,  m * t)
useage_m = cumsum(m * t); useage_m = c(0, useage_m)

for(i in 1:(m * t)){i_reg = region_ij(useage_m, i)
      for(j in 1:(m * t)){j_reg = region_ij(useage_m, j)
          if(i_reg == j_reg & i != j){
              cor_matrix[i, j] = alpha0}else if(i == j){cor_matrix[i, j] = 1
 }else if(i_reg != j_reg){cor_matrix[i,j] = alpha0 * (rho^(abs(i_reg - j_reg)))}}}

Sigma[[1]] = cor_matrix; Sigma[[2]] = cor_matrix

y_ed = simbinCLF(mu = mu, Sigma = Sigma, n = n)


## ----info, results='markup', echo=FALSE-----------------------------------------------------------
sessionInfo()

