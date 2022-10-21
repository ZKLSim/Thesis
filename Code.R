####################################################
################ Loading Data ######################

# load("Mortality.RData")
library(StMoMo)
library(demography)
library(tidyverse)
library(ggplot2)
library(pracma)
library(splines)
library(fpp3)

####################################################

## Function used to extract data 
data_clean <- function(mort, min_age, max_age, min_year, max_year){
  result <- extract.years(extract.ages(mort, ages = min_age:max_age, combine.upper = FALSE),
                          years = min_year:max_year)
  return(result)
}

## Function used for the SDF model

SDF <- function(hmd_data, h, min_age, max_age, min_year, cut, factor, type, imp_fac){
  
  data <- data_clean(hmd_data, min_age, max_age, min_year, cut)
  mort_data <- log(data$rate[[type]])
  
  data_true <- data_clean(hmd_data, min_age, max_age, cut+h, cut+h)
  mort_data_true <- log(data_true$rate[[type]])
  
  Age <- min_age:max_age
  age_all <- imp_fac
  factor <- factor
  h<-h
  
  In_fit <- matrix(0,nrow=factor,ncol=1)
  colnames(In_fit) <- "In-Sample Fit"
  Out_ARIMA <- matrix(0,nrow=factor,ncol=1)
  colnames(Out_ARIMA) <- "ARIMA"
  Out_VAR <- matrix(0,nrow=factor,ncol=1)
  colnames(Out_VAR) <- "VAR"
  
  for(i in seq(1,factor,1)){
    
    if(i == 1){
      # Estimating a,b,k
      ax <- rowMeans(mort_data)
      ax_sp <- smooth.spline(Age, ax, cv = TRUE)
      ax_sp <- ax_sp$y
      err <- (mort_data - ax_sp) %>% as.matrix()
      kappa <- err[1:i,]
      beta <- mrdivide(err, kappa)
      
      bx1_sp <- smooth.spline(Age[2:length(Age)], beta[2:length(Age)], cv = TRUE)
      beta_hat1 <- rbind(as.matrix(beta[i]), as.matrix(bx1_sp$y))
      
      # In
      In_sp <- ax_sp + as.matrix(beta_hat1) %*% kappa
      In_fit[i,] <- (mort_data - In_sp)^2 %>% unlist() %>% mean()
      
      #if(fore_method == "ARIMA"){}
      
      # ARIMA
      k_arima <- auto.arima(kappa, max.p = 10, max.q = 10, start.p = 2, start.q = 2, stepwise = F) %>% forecast(h = h)
      k_fore <- (k_arima$mean)[h]
      Out_ss <- ax_sp + as.matrix(beta_hat1) %*% k_fore
      Out_ARIMA[i,] <- (mort_data_true - Out_ss)^2 %>% unlist() %>% mean()
      
      # VAR
      y <- t(kappa)
      x <- cbind(ones(length(y)-1,1), y[1:(length(y)-1)])
      b <- inv(t(x) %*% x) %*% (t(x) %*% y[2:length(y)])
      intercept <- b[1]
      B <- b[2]
      pred <- matrix(0,i,h+1)
      pred[,1] <- kappa[length(y)]
      
      for(k in seq(1,h,1)){
        pred[,k+1] <- intercept + B %*% pred[,k]
      }
      
      Out_ss_VAR <- ax_sp + as.matrix(beta_hat1) %*% pred[,h+1]
      Out_VAR[i,] <- (mort_data_true - Out_ss_VAR)^2 %>% unlist() %>% mean()
    }
    
    else{
      beta_hat <- matrix(0,length(Age),i)
      Age_2 <- matrix(0,length(Age),1)
      Age_2[(1:i)] <- Age[(age_all[1:i]+1)]
      Age_2[-(1:i)] <- Age[-(age_all[1:i]+1)]
      rate_age <- mort_data
      rate_age[(1:i),] <- mort_data[(age_all[1:i]+1),]
      rate_age[-(1:i),] <- mort_data[-(age_all[1:i]+1),]
      rate_true_age <- mort_data_true
      rate_true_age[(1:i)] <- mort_data_true[(age_all[1:i]+1)]
      rate_true_age[-(1:i)] <- mort_data_true[-(age_all[1:i]+1)]
      alpha_x_age <- rowMeans(rate_age)
      
      ax_sp_age <- smooth.spline(Age_2, alpha_x_age, cv = TRUE)
      ax_sp_age <- ax_sp_age$y
      
      err_age <- (rate_age - ax_sp_age) %>% as.matrix()
      kappa_age <- err_age[1:i,]
      beta_age <- mrdivide(err_age, kappa_age)
      
      # In-sample fit
      
      for(g in seq(1,i,1)){
        bx_sp <- smooth.spline(Age_2[(i+1):length(Age_2)], beta_age[((i+1):length(Age_2)),g], cv = TRUE)
        beta_hat[,g] <- rbind(as.matrix(beta_age[1:i,g]), as.matrix(bx_sp$y))
      }
      In_sp <- (ax_sp_age) + as.matrix(beta_hat) %*% kappa_age
      In_fit[i,] <- (rate_age - In_sp)^2 %>% unlist() %>% mean()
      
      # ARIMA
      k_fore <- matrix(0, i, 1)
      
      for(p in seq(1,i,1)){
        k_arima <- auto.arima(kappa_age[p,], max.p = 10, max.q = 10, start.p = 2, start.q = 2, stepwise = F) %>% forecast(h = h)
        k_fore[p,] <- (k_arima$mean)[h]
      }
      
      Out_ss <- as.numeric(rate_age[,ncol(rate_age)]-as.matrix(beta_hat) %*% kappa_age[,ncol(kappa_age)]) + as.matrix(beta_hat) %*% k_fore[1:i,]
      Out_ARIMA[i,] <- (rate_true_age - Out_ss)^2 %>% unlist() %>% mean()
      
      # VAR
      y <- t(kappa_age)
      x <- cbind(ones(ncol(kappa_age)-1,1), y[1:(ncol(kappa_age)-1),])
      b <- inv(t(x) %*% x) %*% (t(x) %*% y[2:ncol(kappa_age),])
      pred <- matrix(0,i,h+1)
      
      intercept <- as.numeric(b[1,])
      
      B <- t(b[2:(i+1),])
      pred[,1] <- kappa_age[,ncol(kappa_age)]
      
      for(k in seq(1,h,1)){
        pred[,k+1] <- intercept + B %*% pred[,k]
      }
      
      Out_ss_VAR <- as.numeric(rate_age[,ncol(rate_age)]-as.matrix(beta_hat) %*% kappa_age[,ncol(kappa_age)]) + as.matrix(beta_hat) %*% pred[1:i,h+1]
      Out_VAR[i,] <- (rate_true_age - Out_ss_VAR)^2 %>% unlist() %>% mean()
      
    }
  }
  output <- list(In_fit,Out_ARIMA,Out_VAR)
  return(output)
}

## Function used for the Extended model

Extended_LC <- function(hmd_data,h,min_age,max_age,min_year,cut,factor, type){
  data <- data_clean(hmd_data, min_age, max_age, min_year, cut)
  mort_data <- log(data$rate[[type]])
  
  data_true <- data_clean(hmd_data, min_age, max_age, cut+h, cut+h)
  mort_data_true <- log(data_true$rate[[type]])
  
  Age <- min_age:max_age
  betas <- matrix(0,nrow=length(Age), ncol=factor)
  colnames(betas) <- paste0("beta_", 1:factor)
  kappas <- matrix(0, nrow = factor, ncol = ncol(mort_data))
  
  tab_std_in <- matrix(0,nrow=factor,ncol=1)
  rownames(tab_std_in) <- paste0(1:factor, " factor")
  tab_std_out <- matrix(0,nrow=factor,ncol=1)
  rownames(tab_std_out) <- paste0(1:factor, " factor")
  
  alpha_x <- rowMeans(mort_data)
  male_scratch <- mort_data - alpha_x  
  SVD_decom <- svd(male_scratch, nu = nrow(male_scratch),  nv = ncol(male_scratch))
  D <- diag(SVD_decom$d)
  
  for(n in seq(1,factor)){
    betas[,n] <- -1*SVD_decom$u[,n] * D[n,n]
    kappas[n,] <- -1 * SVD_decom$v[,n]
    kappas[n,] <- kappas[n,]*sum(betas[,n])
    kappas[n,] <- kappas[n,]-mean(kappas[n,])
    betas[,n] <- betas[,n]/sum(betas[,n])
    
    # In
    
    if(n==1){
      Fit_in <- alpha_x + betas[,1] %*% t(kappas[1,])
      tab_std_in[n,] <- (mort_data - Fit_in)^2 %>% unlist() %>% mean()
    }
    
    else{  
      Fit_in <- alpha_x + as.matrix(betas[,1:n]) %*% kappas[1:n,]
      tab_std_in[n,] <- (mort_data - Fit_in)^2 %>% unlist() %>% mean()
    }
  }
  
  # Out-of-Sample
  # ARIMA
  k_fore <- matrix(0,nrow = factor, ncol = 1)  
  
  for(p in seq(1,factor,1)){
    k_arima <- auto.arima(kappas[p,], max.p = 10, max.q = 10, start.p = 2, start.q = 2, stepwise = F) %>% forecast(h = h)
    k_fore[p,] <- (k_arima$mean)[h]
    
    if(p==1){
      Fit_out <- alpha_x + betas[,1] %*% t(k_fore[p,])
      tab_std_out[p,] <- (mort_data_true - Fit_out)^2 %>% unlist() %>% mean()
    }
    
    else{
      Fit_out <- alpha_x + as.matrix(betas[,1:p]) %*% k_fore[1:p,]
      tab_std_out[p,] <- (mort_data_true - Fit_out)^2 %>% unlist() %>% mean()
    }
  }
  
  output = list(tab_std_in,tab_std_out)
  return(output)
}

## Function used for the Standard LC model

Standard_LC <- function(hmd_data,h,min_age,max_age,min_year,cut, type){
  
  data <- data_clean(hmd_data, min_age, max_age, min_year, cut)
  mort_data <- log(data$rate[[type]])
  
  data_true <- data_clean(hmd_data, min_age, max_age, cut+h, cut+h)
  mort_data_true <- log(data_true$rate[[type]])
  
  alpha_x <- rowMeans(mort_data)
  male_scratch <- mort_data - alpha_x  
  SVD_decom <- svd(male_scratch, nu = nrow(male_scratch),  nv = ncol(male_scratch))
  D <- diag(SVD_decom$d)
  
  
  betax_1 <- -1 * SVD_decom$u[,1] * D[1,1]
  betax_1_sum <- sum(betax_1)
  betax_1_adj <- betax_1/betax_1_sum
  
  kappa_1 <- -1 * SVD_decom$v[,1]
  kappa_1_hat <- kappa_1 * betax_1_sum
  kappa_1_mean <- mean(kappa_1_hat)
  kappa_1_adj <- kappa_1_hat - kappa_1_mean
  
  # In
  
  fit_in <- alpha_x + betax_1_adj %*% t(kappa_1_adj)
  In_Standard_LC <- ((fit_in - mort_data)^2) %>% unlist() %>% mean()
  
  # Out
  
  mu <- mean(kappa_1_adj[2:length(kappa_1_adj)] - kappa_1_adj[1:(length(kappa_1_adj)-1)])
  
  k_for <- matrix(kappa_1_adj[length(kappa_1_adj)] + h*mu)
  fit_out <- alpha_x + betax_1_adj %*% t(k_for)
  Out_Standard_LC <- ((fit_out - mort_data_true)^2) %>% unlist() %>% mean()
  
  In_Standard_LC <- as.matrix(In_Standard_LC)
  Out_Standard_LC <- as.matrix(Out_Standard_LC)
  colnames(In_Standard_LC) <- "In Sample Fit"
  colnames(Out_Standard_LC) <- "Out-of Sample Fit"
  
  output <- list(In_Standard_LC,Out_Standard_LC)
  return(output)
}

####################################################

### Empirical Findings and Analysis

## Figure 1
min_age <- 0
max_age <- 100
min_year <- 1816
max_year <- 2019
type <- "male"

data <- data_clean(FRdata, min_age, max_age, min_year, max_year)
M_xt_male <- log(data$rate[[type]])
as_tibble(M_xt_male) %>%
  add_column(Age = min_age:max_age, .before = "1950") %>%
  pivot_longer(cols = -Age, names_to = "Year", values_to = "Mortality") %>%
  ggplot(aes(x = Year, y = Mortality, group = Age, col = Age)) +
  geom_line(size = 1.1)+
  scale_color_gradientn(colours = rainbow(10))+
  scale_x_discrete(breaks = seq(from = min_year, to = max_year, by = 50))+
  theme(text = element_text(size = 15))

## Figure 2
data <- data_clean(FRdata, min_age, max_age, 1816, max_year)
mort_data <- log(data$rate[[type]])
data_plot <- cbind(Age,mort_data[,c("1856","1918","1935","2005","2019")])

p1 <- as_tibble(mort_data) %>% 
  add_column(Age = 0:100, .before = "1816") %>% 
  pivot_longer(cols = -Age, names_to = "Year", values_to = "logmx")


as_tibble(data_plot) %>% 
  pivot_longer(cols = -Age, names_to = "Year", values_to = "logmx") %>% 
  ggplot() +
  geom_line(data = p1, aes(y = logmx, x=Age, group = Year), col = "grey",size = 0.8)+
  geom_line(aes(y = logmx, x=Age, group = Year, col = Year),size = 0.8) +
  ylab("log mortality")+
  theme(text = element_text(size = 15))+
  ggtitle("France Log Mortality (1816-2019)")

## Model Performance code

# Input value
min_age <- 0
max_age <- 100
Age <- 0:100
min_year <- 1950
max_year <- 2019
type <- "male"

data <- data_clean(FRdata, min_age, max_age, min_year, max_year)

Mxt_all <- log(data$rate[[type]])
num_year <- 2000-min_year+1
M_xt_male <- Mxt_all[,1:num_year]
M_xt_male_true <- Mxt_all[,-(1:num_year)]
#
h <- 15
factor <- 6
cut <- 2000

# The Standard LC

MSE_SLC <- matrix(0,nrow = 1, ncol = h)

for(w in seq(1,h,1)){
  
  SLC <- matrix(0,nrow = 1, ncol = (max_year-cut-w+1)) 
  SLC_In <- matrix(0,nrow = 1, ncol = (max_year-cut-w+1))
  
  for(o in seq(0,max_year-cut-w,1)){
    SLC[,o+1] <- Standard_LC(FRdata,w,min_age,max_age,min_year,cut+o,type)[[2]]
    SLC_In[,o+1] <- Standard_LC(FRdata,w,min_age,max_age,min_year,cut+o,type)[[1]]
  }
  
  MSE_SLC[,w] <- SLC %>% rowMeans()
  MSE_SLC_In[,w] <- SLC_In %>% rowMeans()
}

# The Extended LC

MSE_ELC <- matrix(0,nrow = factor, ncol = h)
MSE_ELC_In <- matrix(0,nrow = factor, ncol = h)

for(w in seq(1,h,1)){
  
  print(w)
  ELC <- matrix(0,nrow = factor, ncol = (max_year-cut-w+1)) 
  ELC_In <- matrix(0,nrow = factor, ncol = (max_year-cut-w+1)) 
  
  for(o in seq(0,max_year-cut-w,1)){
    Hey <- Extended_LC(FRdata,w,min_age,max_age,min_year,cut+o,factor,type)
    ELC[,o+1] <- Hey[[2]]
    ELC_In[,o+1] <- Hey[[1]]
  }
  
  MSE_ELC[,w] <- ELC %>% rowMeans()
  MSE_ELC_In[,w] <- ELC_In %>% rowMeans()
}

# SDF 

imp_fac <- c(0,99,33,66,21,18)

MSE_SDF <- matrix(0,nrow = factor, ncol = h)
MSE_SDF_VAR <- matrix(0,nrow = factor, ncol = h)
MSE_SDF_In <- matrix(0,nrow = factor, ncol = h)

for(w in seq(1,h,1)){
  print(w)
  
  hh <- matrix(0,nrow = factor, ncol = (max_year-cut-w+1))
  hh_VAR <- matrix(0,nrow = factor, ncol = (max_year-cut-w+1)) 
  hh_In <- matrix(0,nrow = factor, ncol = (max_year-cut-w+1))
  
  for(o in seq(0,max_year-cut-w,1)){
    Hey <- Expand(FRdata, w, min_age, max_age, min_year, cut+o, factor, type, imp_fac)
    hh[,o+1] <- Hey[[2]]
    hh_VAR[,o+1] <- Hey[[3]]
    hh_In[,o+1] <- Hey[[1]]
  }
  
  MSE_SDF[,w] <- hh %>% rowMeans()
  MSE_SDF_VAR[,w] <- hh_VAR %>% rowMeans()
  MSE_SDF_In[,w] <- hh_In %>% rowMeans()
}

# The Hyndman_Ullah Model

MSE_HU <- matrix(0,nrow = factor, ncol = h)
MSE_HU_In <- matrix(0,nrow = factor, ncol = h)

for(u in seq(1,h,1)){
  print(u)
  rob <- matrix(0, nrow = factor, ncol = (max_year-cut-u+1))
  tab_rob_in <- matrix(0, nrow = factor, ncol = (max_year-cut-u+1))
  for(hu in seq(0,max_year-cut-u,1)){
    
    # Hyndman_Ullah In
    data_rob_in <- data_clean(FRdata, min_age, max_age, min_year, cut+hu)
    rate_real_rob <- log(data_rob_in$rate[[type]])
    
    data_rob_out <- data_clean(FRdata, min_age, max_age, cut+hu+u, cut+hu+u)
    rate_rob_out <- log(data_rob_out$rate[[type]])
    
    for(q in seq(1,factor,1)){
      # in sample
      fit <- demography::fdm(data_rob_in, order = q, series = type, method = "M")
      tab_rob_in[q,] <- (fit$fitted[[2]] - rate_real_rob)^2 %>% unlist() %>% mean()
      
      # Hyndman_Ullah Out
      fore <- forecast::forecast(fit, h = u)
      rob[q,hu+1] <- ((log(fore$rate$male))[,u] - rate_rob_out)^2 %>% unlist() %>% mean()
    }
  }
  MSE_HU[,u] <- rob %>% rowMeans()
  MSE_HU_In[,u] <- tab_rob_in %>% rowMeans()
}

# Figure 3

nf <- 6 #where nf takes value from 1 to 6 which representing 1 factor to 6 factors
as_tibble(cbind(MSE_SDF[nf,], MSE_HU[nf,], t(MSE_SLC), MSE_ELC[nf,])) %>% 
  add_column(`h step` = 1:(h), .before="V1") %>% 
  mutate(SDF = V1, Hyndman_Ullah = V2, Standard_LC = V3, Extended_LC = V4) %>% 
  dplyr::select(-V1,-V2,-V3,-V4) %>% 
  pivot_longer(cols = -`h step`,names_to = "Methods", values_to = "MSE") %>% 
  ggplot(aes(x = `h step`, y = MSE, group = Methods, col = Methods)) +
  geom_line(size = 0.8)+
  scale_x_continuous(breaks = seq(from = 1, to = h, by = 1))+
  ggtitle(paste(nf, " factors"))+
  theme(text = element_text(size = 17))

# 30 Years forecasting code
i=5
h=30
imp_fac <- c(0,99,33,66,21)
data <- data_clean(FRdata, min_age, max_age, min_year, max_year)
mort_data <- log(data$rate[[type]])

Age <- min_age:max_age
age_all <- imp_fac
factor <- 5
h<-h

beta_hat <- matrix(0,length(Age),i)
Age_2 <- matrix(0,length(Age),1)
Age_2[(1:i)] <- Age[(age_all[1:i]+1)]
Age_2[-(1:i)] <- Age[-(age_all[1:i]+1)]
rate_age <- mort_data
rate_age[(1:i),] <- mort_data[(age_all[1:i]+1),]
rate_age[-(1:i),] <- mort_data[-(age_all[1:i]+1),]

alpha_x_age <- rowMeans(rate_age)

ax_sp_age <- smooth.spline(Age_2, alpha_x_age, cv = TRUE)
ax_sp_age <- ax_sp_age$y

err_age <- (rate_age - ax_sp_age) %>% as.matrix()
kappa_age <- err_age[1:i,]
beta_age <- mrdivide(err_age, kappa_age)

for(g in seq(1,i,1)){
  bx_sp <- smooth.spline(Age_2[(i+1):length(Age_2)], beta_age[((i+1):length(Age_2)),g], cv = TRUE)
  beta_hat[,g] <- rbind(as.matrix(beta_age[1:i,g]), as.matrix(bx_sp$y))
}

k_fore <- matrix(0, i, h)

for(p in seq(1,i,1)){
  k_arima <- auto.arima(kappa_age[p,], max.p = 10, max.q = 10, start.p = 2, start.q = 2, stepwise = F) %>% forecast(h = h)
  k_fore[p,] <- k_arima$mean
}

For30Y <- as.numeric(rate_age[,ncol(rate_age)]-as.matrix(beta_hat) %*% kappa_age[,ncol(kappa_age)]) + as.matrix(beta_hat) %*% k_fore[1:i,]

# Figure 4
colnames(For30Y) <- paste0(2020:(2020+h-1))

plot50_19 <- as_tibble(mort_data) %>%
  add_column(Age = min_age:max_age, .before = "1950") %>%
  pivot_longer(cols = -Age, names_to = "Year", values_to = "Mortality")

as_tibble(For30Y) %>%
  add_column(Age = Age_2, .before = "2020") %>%
  arrange(Age) %>%
  pivot_longer(cols = -Age, names_to = "Year", values_to = "Mortality") %>% 
  mutate(Year = as.numeric(Year))%>% 
  ggplot(aes(x = Age, y = Mortality, group = Year, col = Year)) +
  geom_line(size = 0.8)+
  scale_color_gradientn(colours = rainbow(10))+
  geom_line(data = plot50_19, color = "black" ,alpha = 0.09 ,aes(x = Age, y = Mortality, group = as.numeric(Year)))+
  ylab("log(mxt)") +
  theme(text = element_text(size = 15))


# Figure 5

kappa_empirical <- kappa_age
nf <- 5 #where nf takes value from 1 to 5 which representing 1st kappa to 5th kappa

first <- auto.arima(kappa_empirical[nf,], max.p = 10, max.q = 10, start.p = 2, start.q = 2, stepwise = F)
fore1 <- first %>% forecast(h = 30)
`95low` <- fore1$lower[,2]
`95up` <- fore1$upper[,2]
point <- fore1$mean
Year30 <- max_year:(max_year+29)
interval1 <- cbind(Year30, `95low`, `95up`, point) %>% as_tibble()
fitted_first <- kappa_empirical[nf,]
arima_first <- first$fitted


as_tibble(fitted_first) %>%
  add_column(arima_first = arima_first, Year = min_year:max_year, .before = "value") %>%
  mutate(SDF_fitted = value, Arima_fitted = arima_first) %>% 
  dplyr::select(-value,-arima_first) %>% 
  pivot_longer(cols = -Year, names_to = "Fit", values_to = "Fitted Latent Factor") %>% 
  ggplot() +
  geom_line(aes(x = Year, y = `Fitted Latent Factor`, group = Fit, col = Fit),size = 0.8)+
  geom_line(data = interval1 ,size = 0.8 ,aes(y = `95low`, x = Year30), linetype="dotted")+
  geom_line(data = interval1 ,size = 0.8 ,aes(y = `95up`, x = Year30), linetype="dotted")+
  geom_line(data = interval1 ,size = 0.8 ,aes(y = point, x = Year30), linetype="twodash")+
  ylab("Fifth Latent Factor") +
  theme(text = element_text(size = 15))


# Figure 6

data <- data_clean(FRdata, min_age, max_age, min_year, max_year)
mort_data <- log(data$rate[[type]])
factor <- 5 #Comparing 5 factors classical mortality model with 5 factors SDF mdoel

Age <- min_age:max_age
betas <- matrix(0,nrow=length(Age), ncol=factor)
colnames(betas) <- paste0("beta_", 1:factor)
kappas <- matrix(0, nrow = factor, ncol = ncol(mort_data))

alpha_x <- rowMeans(mort_data)
male_scratch <- mort_data - alpha_x  
SVD_decom <- svd(male_scratch, nu = nrow(male_scratch),  nv = ncol(male_scratch))
D <- diag(SVD_decom$d)

for(n in seq(1,factor)){
  betas[,n] <- -1*SVD_decom$u[,n] * D[n,n]
  kappas[n,] <- -1 * SVD_decom$v[,n]
  kappas[n,] <- kappas[n,]*sum(betas[,n])
  kappas[n,] <- kappas[n,]-mean(kappas[n,])
  betas[,n] <- betas[,n]/sum(betas[,n])
}

nf <- 1 #where nf takes value from 1 to 5 which representing 1st kappa to 5th kappa

as_tibble(kappas[nf,]) %>% 
  add_column(Year = min_year:max_year) %>% 
  ggplot()+
  geom_line(aes(x = Year, y = value),size = 0.8)+
  ylab("Factor") +
  theme(text = element_text(size = 15))

# Correlation code for table 4

cor_matrix <- matrix(0,nrow = 5, ncol = 5)

for(i in seq(1,5)){
  
  for(j in seq(1,5)){
    cor_matrix[i,j] <- cor(kappa_empirical[i,],kappas[j,])
  }
  
}

# Robust Test of the SDF model

# Robustness test
SDF1 <- c(0,99,33,66,21,18)
SDF2 <- c(5,96,38,61,22,19)
SDF3 <- c(8,91,41,70,28,15)
factor <- 5

# where `imp_fac_diff` take values of SDF2 when second age set is used, take values of SDF3 when third age set is used.

MSE_SDF_diff <- matrix(0,nrow = factor, ncol = h)
MSE_SDF_In_diff <- matrix(0,nrow = factor, ncol = h)

for(w in seq(1,h,1)){
  print(w)
  
  hh_diff <- matrix(0,nrow = factor, ncol = (max_year-cut-w+1))
  hh_In_diff <- matrix(0,nrow = factor, ncol = (max_year-cut-w+1))
  
  for(o in seq(0,max_year-cut-w,1)){
    Hey_diff <- Expand(FRdata, w, min_age, max_age, min_year, cut+o, factor, type, imp_fac_diff)
    hh_diff[,o+1] <- Hey_diff[[2]]
    hh_In_diff[,o+1] <- Hey_diff[[1]]
  }
  
  MSE_SDF_diff[,w] <- hh_diff %>% rowMeans()
  MSE_SDF_In_diff[,w] <- hh_In_diff %>% rowMeans()
}









