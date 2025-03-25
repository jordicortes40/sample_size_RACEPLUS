rm(list=ls())

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
# library(copula)
library(bindata)   # Generate binomial data
library(statpsych) # paired proportions confidence interval


#-------------------------------------------------------------------------------
# Parameters of the joint distribution
# Table of empirical concordance
#               RACEPLUS_acierta
# RACE_acierta  No Yes
#          No  253 109
#         Yes   88 653
# This table is calculated in the following path:
# M:\La meva unitat\Consultoria\Can Ruti\RACEPLUS\RACEPLUS - Sample size\Sample size validation - Paired\scripts\Cohort retrospectiva
# This table represents LVO correct guess - YES/NO with RACEPLUS versus 
# LVO correct guess - YES/NO with RACE 
#-------------------------------------------------------------------------------
tab_emp <- matrix(c(253,109,88,653), nrow=2, byrow=TRUE, 
                  dimnames = list(c('RACE - NO',    'RACE - YES'),
                                  c('RACEPLUS - NO','RACEPLUS - YES')))

# To check if theoretical correlation is correct:
# source('estimate_correlation between_RACE_RACEPLUS.R')
# Conclusion: rho = 0.58

# Theoretical correlation:
# https://web.pdx.edu/~newsomj/uvclass/ho_chisq.pdf
chisq.test(tab_emp)$sta/sum(tab_emp) # --> 0.587

#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------
##-- Parameters for sample size calculation ------------------------------------
DELTA <- seq(-0.01,0.03,0.005)          # Difference proportion
p0  <- 0.67                             # Proportion of correct guess RACE
P1  <- p0 + DELTA                       # Proportion of correct guess RACEPLUS
rho <- 0.58                             # Correlation between guesses (previously calculated)
alpha    <- 0.05                        # Probability of type I error
ss_total <- 3500                        # Planned total sample size
P_INT    <- c(0.1, 0.25, 0.33, 0.5)     # Proportion of the sample for doing the interim
SS_INT   <- round(p_int * ss_total)     # Sample size Interim
nsim     <- 10000                       # Number of simulations  


#-------------------------------------------------------------------------------
# A priori distribution
#-------------------------------------------------------------------------------
dapriori <- function(x,n,p1,p2){
  expected_value <-  p2 - p1
  variance       <-  (p1 * (1 - p1))/n +
    (p2 * (1 - p2))/n  
  dnorm(x,expected_value,sqrt(variance))
}
##-- Example
# n <- sum(tab_emp)
# p_correct_RACE       <- sum(tab_emp['RACE - YES',])/n     # ~0.67  # Accuracy with RACE
# p_correct_RACEPLUS   <- sum(tab_emp[,'RACEPLUS - YES'])/n # ~0.69  # Accuracy with RACE-PLUS
# dapriori(0, n = n , p1 = p_correct_RACE, p2 = p_correct_RACEPLUS)
# dapriori(0.02, n = n , p1 = p_correct_RACE, p2 = p_correct_RACEPLUS)
# dapriori(0.05, n = n , p1 = p_correct_RACE, p2 = p_correct_RACEPLUS)
# p_apriori <- matrix(nrow=length(DELTA),ncol=length(DELTA))
# Asumiendo cada delta como real, calculo la distribución a priori precisamente en los deltas
p_correct_RACE       <- sum(tab_emp['RACE - YES',])/n     # ~0.67  # Accuracy with RACE
p_correct_RACEPLUS   <- sum(tab_emp[,'RACEPLUS - YES'])/n # ~0.69  # Accuracy with RACE-PLUS
p_apriori <- dapriori(DELTA, n = sum(tab_emp), p1 = p_correct_RACE, p2 = p_correct_RACEPLUS)  
names(p_apriori) <- paste0('deltaReal=', DELTA)
p_apriori <- p_apriori/sum(p_apriori)  # Standarize each row to sum 1



#-------------------------------------------------------------------------------
# Likelihood
#-------------------------------------------------------------------------------
##-- Objects to save results ---------------------------------------------------
final_mcnemar <- final_dif <- final_lower_ci <- matrix(nrow = nsim, ncol = 1)
int_mcnemar   <- int_dif   <- int_lower_ci   <- matrix(nrow = nsim, ncol = length(ss_int))  

ls_final_mcnemar <- ls_final_dif <- ls_final_lower_ci <- rep(list(final_mcnemar),length(DELTA))
ls_int_mcnemar   <- ls_int_dif   <- ls_int_lower_ci   <- rep(list(int_mcnemar),length(DELTA))  


##-- Simulation ----------------------------------------------------------------

set.seed(12345)

for(delta in DELTA){
  cat('delta:',delta,'\n')
  p1 <- p0 + delta
  
  for (k in 1:nsim){
    # cat('itera:',k,'\n')
    
    final_data <- rmvbin(ss_total, c(p0,p1), bincorr=(1-rho)*diag(2)+rho)
    final_table  <- table(final_data[,1], final_data[,2])
    
    final_mc_test <- mcnemar.test(x = final_table)
    final_mcnemar[k,1] <- final_mc_test$p.value
    final_dif[k,1]     <- diff(colSums(final_data)/ss_total)
    final_lower_ci[k,1]<- ci.prop.ps(alpha, final_table[2,2], final_table[2,1], final_table[1,2], final_table[1,1])[1,"LL"]
    
    for(ssi in SS_INT){
      int_data   <- final_data[1:ssi,]
      int_table  <- table(int_data[,1], int_data[,2])
      
      int_mc_test <- mcnemar.test(x = int_table)
      int_mcnemar[k,SS_INT==ssi] <- int_mc_test$p.value
      int_dif[k,SS_INT==ssi]     <- diff(colSums(int_data)/ssi)
      int_lower_ci[k,SS_INT==ssi]<- ci.prop.ps(alpha, int_table[2,2], int_table[2,1], int_table[1,2], int_table[1,1])[1,"LL"]
    }
  }
  
  ls_final_mcnemar[[which(DELTA==delta)]]  <- final_mcnemar
  ls_final_dif[[which(DELTA==delta)]]      <- final_dif
  ls_final_lower_ci[[which(DELTA==delta)]] <- final_lower_ci
  
  ls_int_mcnemar[[which(DELTA==delta)]]    <- int_mcnemar
  ls_int_dif[[which(DELTA==delta)]]        <- int_dif
  ls_int_lower_ci[[which(DELTA==delta)]]   <- int_lower_ci
  
}


##-- Plots for checking --------------------------------------------------------
par(mfrow=c(3,3))
##-- Difference proportion at the end
for(i in 1:9){
  hist(ls_final_dif[[i]][,1],xlim=c(-0.05,0.05), main = paste0('Final. delta = ',DELTA[i]), freq = FALSE)
}
##-- Difference proportion at the interimn
for(i in 1:9){
  hist(ls_int_dif[[i]][,1],xlim=c(-0.05,0.05), main = paste0('Interim (n =',SS_INT[1],') delta = ',DELTA[i]), freq = FALSE)  
}
##-- Difference proportion at the end
for(i in 1:9){
  hist(ls_int_dif[[i]][,4],xlim=c(-0.05,0.05), main = paste0('Interim (n=',SS_INT[4],') delta = ',DELTA[i]), freq = FALSE)  
}



#-------------------------------------------------------------------------------
# Scenario: success --> to prove RACEPLUS > RACE
#-------------------------------------------------------------------------------
n_cutpoints0 <- 52  # number of cutpoints to assess
min_dif <- 0
p_confidence <- 0.8

##-- Common cutpoints
# cutpoint0 <- seq(min(unlist(ls_int_dif)), 
#                  0.05,
#                  max(unlist(ls_int_dif)),
#                  length = n_cutpoints0)
cutpoint0 <- seq(-0.05, 0.05, length = n_cutpoints0)
cutpoint <- cutpoint0[-c(1,length(cutpoint0))]
n_cutpoints <- n_cutpoints0 - 2
prop_success <- matrix(nrow=n_cutpoints, ncol = length(SS_INT)) 

for(ssi in SS_INT){
  cat('ssi:',ssi,'\n')
  j <- which(SS_INT==ssi)
  
  for(k in 1:n_cutpoints){
    
    p_succes <- 0
    
    for(delta in DELTA){
      int_dif        <- ls_int_dif[[which(DELTA==delta)]]
      final_lower_ci <- ls_final_lower_ci[[which(DELTA==delta)]]
      
      apriori <- p_apriori[DELTA==delta]
      
      sel <- int_dif[,j] >= cutpoint[k]
      likelihood <- sum(final_lower_ci[sel,1] > min_dif)/ length(final_lower_ci[sel,1])
      # if(is.na(likelihood)) likelihood <- 0
      
      p_succes <- p_succes +  likelihood * apriori
      cat('ssi:',ssi,'delta:',delta,'cutpoint:',cutpoint[k],'likelihood:',likelihood,'apriori:',apriori,'\n')
    }
    
    prop_success[k,j] <- p_succes

  }
}
colnames(prop_success) <- paste0('ss_int=',ss_int)
rownames(prop_success) <- paste0('cutp=',round(cutpoint,3))


##-- Plot
par(mfrow=c(1,1))
plot(cutpoint, prop_success[,1], type='l', 
     xlim=c(0,0.05),
     ylim=c(0.4,1),
     xlab='Cutpoint for proportion difference in the interim analysis',
     ylab='Probability of succes at the end of the study',
     main = 'Prob. of demostrating superiority of RACEPLUS over RACE scale',
     sub = 'Calculations according to interim sample size and the cutpoint',
     col.sub = "orange")
abline(h=p_confidence,lty=2,lwd=2,col=6)
for (i in 2:length(ss_int)){
  lines(cutpoint,prop_success[,i], type='l', col=i)
}
legend('topleft', legend = colnames(prop_success), col=1:5, lty=1)  
abline(v=0.025)
save.image("Report/Rdata/sample_size_simulation.RData")



#