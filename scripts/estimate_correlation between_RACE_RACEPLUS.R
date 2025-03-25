p11_emp <- sum(tab_emp['RACE - YES','RACEPLUS - YES'])/sum(tab_emp)
p_correct_RACE       <- sum(tab_emp['RACE - YES',])/    sum(tab_emp) # ~0.67  # Accuracy with RACE
p_correct_RACEPLUS   <- sum(tab_emp[,'RACEPLUS - YES'])/sum(tab_emp) # ~0.70  # Accuracy with RACE-PLUS


# Create one pair of correlated binomial values
set.seed(12345)
##-- First iteration
RHO <- seq(0.1,0.9,0.1)
size <- 10000
p11_teo <- c()
for(rho in RHO){
  trials <- rmvbin(size, c(p_correct_RACE,p_correct_RACEPLUS), bincorr=(1-rho)*diag(2)+rho)
  colSums(trials)
  tab_teo <- table(trials[,1],trials[,2])
  p11_teo[RHO==rho] <- prop.table(tab_teo)[2,2]
}
plot(RHO,p11_teo,type='b')
abline(h=p11_emp,lty=2,col='grey')
abline(v=RHO,lty=2,col='grey')

##-- 2nd iteration
RHO <- seq(0.5,0.6,0.01)
size <- 100000
p11_teo <- c()
for(rho in RHO){
  trials <- rmvbin(size, c(p_correct_RACE,p_correct_RACEPLUS), bincorr=(1-rho)*diag(2)+rho)
  colSums(trials)
  tab_teo <- table(trials[,1],trials[,2])
  p11_teo[RHO==rho] <- prop.table(tab_teo)[2,2]
}
plot(RHO,p11_teo,type='b')
abline(h=p11_emp,lty=2,col='grey')
abline(v=RHO,lty=2,col='grey')

##-- 3rd iteration
RHO <- seq(0.57,0.59,0.005)
size <- 500000
p11_teo <- c()
for(rho in RHO){
  trials <- rmvbin(size, c(p_correct_RACE,p_correct_RACEPLUS), bincorr=(1-rho)*diag(2)+rho)
  colSums(trials)
  tab_teo <- table(trials[,1],trials[,2])
  p11_teo[RHO==rho] <- prop.table(tab_teo)[2,2]
}
plot(RHO,p11_teo,type='b')
abline(h=p11_emp,lty=2,col='grey')
abline(v=RHO,lty=2,col='grey')