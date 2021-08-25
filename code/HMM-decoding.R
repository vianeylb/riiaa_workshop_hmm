########################################################################
############################  HMM Decoding ##############################
########################################################################

# Vamos a ver como implementar los algoritmos de Decoding. Para ello, vamos a: 
# (i) simular datos
# (ii) leave one time-series out CV: FB y Viterbi  

library(tidyverse)
library(mvtnorm)
library(ggplot2)
library(MCMCpack)
library(reshape2)
source('extrafunctions.R')

########################
#### 1. Simulacion #####
########################
#  HMM con 4 estados (J=4) y con observaciones multivariadas Normales

# probabilidades iniciales
del <- c(0.2, 0.1, 0.6, 0.1)

# matriz de transición
tpm <- matrix(c(0.95, 0.04, 0.005, 0.005, 0.01, 0.98, 0.005, 0.005, 0.01, 0.15,
                0.75, 0.09, 0.02, 0.01, 0.11, 0.86), byrow = T, nrow = 4)

# Simulamos la secuencia de estados para T=50 y n=500
states <- matrix(NA, nrow = 500, ncol = 50)
for (i in 1:50) {
  for (j in 1:500) {
    if (j == 1) {
      states[1, i] <- sample(x = 1:4, size = 1, prob = del)
    } else {
      states[j, i] <- sample(x = 1:4, size = 1, prob = tpm[states[j -
                                                                    1, i], ])
    }
  }
}

# Definimos las medias y las matrices de covarianza de cada estado

mus <- list()
vars <- list()
# state 1
mus[[1]] <- c(0, -7, 10)
vars[[1]] <- riwish(3, matrix(c(3, 0.3, 2, 0.3, 1, 0.5, 2, 0.1, 2), 3, 3))
# state 2
mus[[2]] <- c(0, -7, 11)
vars[[2]] <- riwish(3, matrix(c(3, 0.3, 2, 0.3, 1, 0.5, 0.4, 0.1, 3), 3, 3))
# state 3
mus[[3]] <- c(-2, -2, 3)
vars[[3]] <- riwish(3, matrix(c(3, 0.1, 2, 0.3, 1, 0.5, 2, 0.1, 2), 3, 3))
# state 4
mus[[4]] <- c(5, 5, 1)
vars[[4]] <- riwish(3, matrix(c(1, 0.3, 0.2, 0.3, 2, 0.1, 0.5, 0.1, 2), 3, 3))



# Simulamos las observaciones
obs <- list()
for (k in 1:50) obs[[k]] <- matrix(NA, nrow = 500, ncol = 3)
for (i in 1:50) {
  for (j in 1:500) {
    obs[[i]][j, ] <- rmvnorm(n = 1, mean = mus[[states[j, i]]], sigma = vars[[states[j,
                                                                                     i]]])
  }
}


## Guadamos todo en un data frame 
tracks <- paste("T", 1:50, sep = "")
data <- data.frame(obs[[1]], state = states[, 1], timeseries = rep(tracks[1],
                                                                   500),time=seq(1,500))
for (j in 2:50) 
{data <- rbind(data, data.frame(obs[[j]], state = states[, j], timeseries = rep(tracks[j], 500),time=seq(1,500)))}

data$state=as.factor(data$state)

ggplot(data %>% filter(timeseries=='T1'))+ geom_point(aes(x=time,y=state,col=state))+theme_bw()


###############################
#### 2. Cross-validation  #####
###############################

# a. Ajustar el modelo con todas las series temporales menos la i-ésima y luego predecir los estados de la serie i
# b. Estimar el error de clasificacion
# Eso lo vamos a hacer con local y global decoding

## a)

# Estimacion de las probas iniciales
delta_est <- list()
for (i in 1:50) {
  delta_est[[i]] <- table(states[1, -i])/49
}
# Estimacion de las medias
mu_est <- list()
for (i in 1:50) {
  mu_est[[i]] <- data %>% filter(timeseries != tracks[i]) %>% group_by(state) %>%
    summarize(x1m = mean(X1), x2m = mean(X2), x3m = mean(X3))
}
# estimacion de las matrices de covarianza
cov_est <- list()
for (j in 1:50) {
  temp1 <- filter(data, timeseries != tracks[i] & state == 1)
  c1 <- cov(temp1[, c("X1", "X2", "X3")])
  temp2 <- filter(data, timeseries != tracks[i] & state == 2)
  c2 <- cov(temp2[, c("X1", "X2", "X3")])
  temp3 <- filter(data, timeseries != tracks[i] & state == 3)
  c3 <- cov(temp3[, c("X1", "X2", "X3")])
  temp4 <- filter(data, timeseries != tracks[i] & state == 4)
  c4 <- cov(temp4[, c("X1", "X2", "X3")])
  cov_est[[j]] <- data.frame(state = rep(1:4, each = 3), rbind(c1, c2, c3,
                                                               c4))
}

# estimacion de la matriz de transicion 
tpm_est <- list()
for (j in 1:50) {
  tpm_est[[j]] <- diag(4)
  n <- length(states[, j])
  smat <- matrix(0, nrow = 4, ncol = 4)
  for (k in c(1:50)[-j]) {
    w1 <- which(states[1:(n - 1), k] == 1)
    s1 <- states[w1, k] - states[w1 + 1, k]
   
    w2 <- which(states[1:(n - 1), k] == 2)
    s2 <- states[w2, k] - states[w2 + 1, k]
    
    w3 <- which(states[1:(n - 1), k] == 3)
    s3 <- states[w3, k] - states[w3 + 1, k]
    
    w4 <- which(states[1:(n - 1), k] == 4)
    s4 <- states[w4, k] - states[w4 + 1, k]
   
    for(m in 1:4)
    {
      smat[1, m] = smat[1, m]+ length(which(s1==(1-m)))
      smat[2, m] = smat[2, m]+ length(which(s2==(2-m)))
      smat[3, m] = smat[3, m]+ length(which(s3==(3-m)))
      smat[4, m] = smat[4, m]+ length(which(s4==(4-m)))
      
    }
  }
  wsm <- rowSums(smat)
  tpm_est[[j]][1,]<-smat[1,]/wsm[1]
  tpm_est[[j]][2,]<-smat[2,]/wsm[2]
  tpm_est[[j]][3,]<-smat[3,]/wsm[3]
  tpm_est[[j]][4,]<-smat[4,]/wsm[4]  
}

###############################
#### b.1 local decoding  #####
###############################

data$predstateFB <- NA
ACC.FB=numeric(50)
CE=numeric(50)

for(u in 1:50){
  tempdat <- filter(data, timeseries==tracks[u])
  # calculating allprobs matrix - observations evaluated under each state-dependent MVN
  
  mod=list()
  mod$m=4
  mod$delta=delta_est[[u]]
  mod$tpm=tpm_est[[u]]
  
  
  allprobs <- matrix(NA, nrow=500, 4)
  for(j in 1:500) {
    for(k in 1:4) {
      allprobs[j,k] <- dmvnorm(x=tempdat[j,c("X1", "X2", "X3")],
                               mean = unlist(mu_est[[u]][k,2:4]),
                               sigma = as.matrix(filter(cov_est[[u]],
                                                        state==k)[,c("X1", "X2", "X3")]))
    }
    
    
  }
  
  Pstates=norm.HMM.state_probs(x=tempdat , mod ,allprobs)
  Pred.states=apply(Pstates,2,which.max)
  data$predstateFB[which(data$timeseries==tracks[u])]=Pred.states
  
  # estimacion de accuracy
  
  AccFB[u]=mean(tempdat$state==Pred.states)
  quienes= which(tempdat$state==Pred.states)
  CE[u]=cross.entropy(Pstates,tempdat$state,quienes)
  }



# Matriz de confusion
acast(data, state~predstateFB)

###############################
#### b.2 global decoding  #####
###############################

data$predstateV <- NA
ACC.V=numeric(50)

for(u in 1:50){
  tempdat <- filter(data, timeseries==tracks[u])
  # calculating allprobs matrix - observations evaluated under each state-dependent MVN
  allprobs <- matrix(NA, nrow=500, 4)
  for(j in 1:500) {
    for(k in 1:4) {
      allprobs[j,k] <- dmvnorm(x=tempdat[j,c("X1", "X2", "X3")],
                               mean = unlist(mu_est[[u]][k,2:4]),
                               sigma = as.matrix(filter(cov_est[[u]],
                                                        state==k)[,c("X1", "X2", "X3")]))
    }
  }
  #HMM.viterbi - Viterbi Algorithm
  data$predstateV[which(data$timeseries==tracks[u])] <- HMM.viterbi(x=tempdat, m = 4,
                                                                   gamma = tpm_est[[u]],
                                                                   allprobs = allprobs,
                                                                   delta = delta_est[[u]])
  # estimacion de accuracy
  AccV[u]=mean(tempdat$state==data$predstateV[which(data$timeseries==tracks[u])])
  
  }


# Matriz de confusion
acast(data, state~predstateV)

# Resultados
Algoritmo=c(rep('Viterbi',50),rep('FB',50))
Acc=c(AccV,AccFB)
CE.=c(rep(NA,50),CE)
Resultados=data.frame(Algoritmo,Acc,CE.)


# Accuracy
ggplot(Resultados)+geom_boxplot(aes(x=Algoritmo,y=Acc,col=Algoritmo))+theme_bw()
# CE (FB)
ggplot(Resultados %>% filter(Algoritmo=="FB"))+geom_histogram(aes(CE.))+theme_bw()


