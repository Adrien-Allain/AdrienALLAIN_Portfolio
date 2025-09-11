#Library

library(readr)
library(lubridate)
library(ggplot2)
library(tseries)
library(forecast)
library(RJDemetra)
library(dplyr)
library(zoo)
library(tidyr)
library(smooth)
library(TSA)
library(seastests)
library(leaps)

#Import bases
google_complet_train <- read_delim("data/google_complet_train.csv")
google_complet_train <- google_complet_train[,-1]

obs_google <- read_csv("data/obs_google.csv")
obs_google <- obs_google[,-1]

data_entrainement <- read_delim("data/data_entrainement.csv")
data_entrainement <- data_entrainement[,-1]
realisee <- read_delim("data/realisee.csv")
realisee <- realisee[,-1]

data_entrainement_diff <- read_delim("data/data_entrainement_diff.csv")
data_entrainement_diff <- data_entrainement_diff[,-1]

inc_ts <- ts(data_entrainement_diff$inc, start = c(2004,2), frequency = 12)
zz <- ts(realisee$inc, start = c(2024, 04), frequency = 12)
y_obs <- realisee[[2]]

#########============= Visualisation =============#########
## Série inc 
ggplot(data_entrainement, aes(x = month, y = inc)) +
  geom_line(color = "#3366CC", size = 1) +
  geom_point(color = "#3366CC", size = 3, alpha = 0.7)  +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    date_minor_breaks = "3 months"
  ) +
  labs(
    x = "Années",
    y = "Nombre de cas"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 16, vjust = 1, hjust = 0.5),
    plot.subtitle = element_text(size = 12, color = "darkgrey"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    panel.background = element_rect(fill = "white")
  )
## Heatmap
###Fonction obtenue par IA pour transformé la ts en df
ts_to_df <- function(ts_obj) {
  freq <- frequency(ts_obj)
  start_year <- start(ts_obj)[1]
  start_period <- start(ts_obj)[2]
  
  if (freq == 12) {  
    dates <- seq.Date(from = as.Date(sprintf("%d-%02d-01", start_year, start_period)), 
                      by = "month", length.out = length(ts_obj))
  } else if (freq == 52) { 
    dates <- seq.Date(from = as.Date(sprintf("%d-01-01", start_year)) + weeks(start_period - 1),
                      by = "week", length.out = length(ts_obj))
  } else {
    stop("Fréquence non prise en charge")
  }
  
  tibble(date = dates,
         year = year(dates),
         month = month(dates, label = TRUE, abbr = TRUE),
         value = as.numeric(ts_obj))
}

df_inc <- ts_to_df(inc_ts)
## Figure
ggplot(df_inc, aes(x = month, y = factor(year), fill = value)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  labs(
    x = "Mois", y = "Année") +
  theme_classic()

## Season plot

ts_complet <- ts(c(inc_ts, zz), 
                  start = start(inc_ts), 
                  frequency = frequency(inc_ts))
ggseasonplot(ts_complet)+
  labs(title ="")+
  theme_bw()

## Series Google

par(mfrow = c(ceiling((ncol(data_entrainement_diff)-2)/2), 2))  

for (name in names(data_entrainement_diff)[3:ncol(data_entrainement_diff)]) {
  plot.ts(data_entrainement_diff[[name]], 
          main = paste("Google:", name),
          ylab = "",
          col = "#3366CC", lwd = 2)
}

#########============= Statistiques descriptives et saisonnalité =============#########

base_quanti <- data_entrainement[,-1]

#Test de Dickey-Fuller

adf.test(inc_ts)

adf_test <- sapply(base_quanti,adf.test)
adf_test

### Différenciation

Cephalees_diff <- diff(data_entrainement$Cephalees)
Frissons_diff <- diff(data_entrainement$Frissons)

adf.test(Frissons_diff)
adf.test(Cephalees_diff)

###Base différenciée
data_entrainement_diff <- data_entrainement[-1,]

data_entrainement_diff$Cephalees <- Cephalees_diff
data_entrainement_diff$Frissons <- Frissons_diff

write.csv(data_entrainement_diff, file = "data/data_entrainement_diff.csv")

#Statistiques descriptives

stats <- function(x){
  list(
    summary = summary(x),
    skewness = TSA::skewness(x),
    kurtosis = TSA::kurtosis(x),
    sd = sd(x),
    var = var(x),
    shapiro.test(x)
  )
}

stats_inc <- stats(inc_ts)
shapiro.test(inc_ts) #Distribution non normale

# Corrélation

corr <- cor(base_quanti_diff, method = "s")
corrplot::corrplot(corr, 
                   method = "shade",     
                   type = "upper",        
                   tl.col = "black",  
                   tl.srt = 45,         
                   tl.cex = 1.1,    
                   cl.cex = 1,          
                   diag = FALSE,          
                   addCoef.col = "black",  
                   number.cex = 1,       
                   col = colorRampPalette(c("#6D9EC1", "white", "#E46726"))(200), 
                   addgrid.col = "lightgray",
                   mar = c(0, 0, 1, 0)
)


# ACF PACF

par(mfrow = c(1,2))
Acf(inc_ts,  lwd = 1, main="")
Pacf(inc_ts,  lwd = 1, main="")


# Périodogramme

dyy <- diff(inc_ts, differences =1) 
par(mfrow = c(1,2))
periodogram(inc_ts, main = "Periodogram de la série brute")
periodogram(dyy)

# Tests saisonnalité

## Seasonal dummies
sd <- seasdum(inc_ts)
show(sd)

## Webel-Ollech test 
wot <- combined_test(inc_ts)
show(wot)


# Décomposition
## X13
myspec <- x13_spec("RSA5c")
mysax13 <- x13(inc_ts, myspec)
summary(mysax13$regarima)
mysax13

plot(mysax13)


#########============= Prévisions Univariées =============#########


## Naïf saisonnier (benchmark)
prev_naive <- snaive(inc_ts, h = 12)

## STL
fitstl = stlm(inc_ts)
prevstl <- forecast(fitstl,12)

## STS
fitsts <- StructTS(inc_ts)
prevsts <- forecast(fitsts, h = 12)

## X13-ARIMA-SEATS
myregx13 <- regarima_x13(inc_ts, spec ="RG5c")
prev_x13_arima <- matrix(myregx13$forecast[1:12])
prev_x13_arima_ts <- ts(prev_x13_arima, start = c(2024, 4), frequency = 12)

## BAGGED
fitbag <- baggedModel(inc_ts)
prevbag <- forecast(fitbag,12)

## SSARIMA
fitssarima <- auto.ssarima(inc_ts, 
                           lags=c(1,12),
                           orders=list(ar=c(3,3),
                                       i=(2),
                                       ma=c(3,3),
                                       select=TRUE))

prevssarima <- forecast(fitssarima, h=12, level = 0.90)


##SARIMA

fit_sarima <- auto.arima(inc_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
prevsarima <- forecast(fit_sarima, h = 12)

### Plot des prévisions
autoplot(prev_naive$mean, series = "Naïve", PI = FALSE) + 
  autolayer(prevstl$mean, series = "STL", PI = FALSE) +
  autolayer(prevsts$mean, series = "STS", PI = FALSE) +
  autolayer(prev_x13_arima_ts, series = "x13_ARIMA", PI = FALSE) +
  autolayer(prevbag$mean, series = "BAGGED", PI = FALSE) +
  autolayer(prevsarima$mean, series = "SARIMA", PI = FALSE) +
  autolayer(prevssarima$mean, series = "SSARIMA", PI = FALSE) +
  labs(
    x = "Temps", y = "Valeur des prévisions") +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),labels = function(x) format(as.yearmon(x), "%b %Y")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Lissage expo

## HW
m <- HoltWinters(inc_ts,seasonal="mul")
prev_HW <- forecast(m, h = 12)

## ETS
fitets <- ets(inc_ts, allow.multiplicative.trend = TRUE)
prevets <- forecast(fitets,12)

# BATS
decomp_tbats <- tbats(inc_ts)
prevtbats <- forecast(decomp_tbats,12)

## ADAM ETS
fitadam1 <- auto.adam(inc_ts, model="ZZZ", lags=c(1,12), select=TRUE)
prevadam1 <- forecast(fitadam1,12)


## ADAM ETS SARIMA

fitadam3 <- auto.adam(inc_ts, model="ZZN", lags=c(1,12), orders=list(ar=c(3,3), i=(2), ma=c(3,3), select=TRUE))
prevadam_sarima <- forecast(fitadam3, h=12, level = 0.90)

### Plot des prévisions lissage expo

autoplot(prev_HW$mean, series = "Holt-Winters") +
  autolayer(prevets$mean, series = "ETS") +
  autolayer(prevtbats$mean, series = "TBATS") +
  autolayer(prevadam1$mean, series = "ADAM ETS") +
  autolayer(prevadam_sarima$mean, series = "ADAM ETS SARIMA")+
  labs(
    x = "Temps", y = "Valeur des prévisions")+
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
                     labels = function(x) format(as.yearmon(x), "%b %Y")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


# Plot prévision univariéss avec observations

autoplot(zz, series = "Réalisée",linetype = "dashed" ) +
  autolayer(prev_HW$mean, series = "Holt-Winters", linetype = "dashed") +
  autolayer(prevets$mean, series = "ETS") +
  autolayer(prevtbats$mean, series = "TBATS") +
  autolayer(prevadam1$mean, series = "ADAM ETS") +
  autolayer(prevadam_sarima$mean, series = "ADAM ETS SARIMA") +
  autolayer(prev_naive$mean, series = "Naïve") +
  autolayer(prevstl$mean, series = "STL") +
  autolayer(prevsts$mean, series = "STS") +
  autolayer(prev_x13_arima_ts, series = "X13_ARIMA") +
  autolayer(prevbag$mean, series = "BAGGED") +
  autolayer(prevsarima$mean, series = "SARIMA") +
  autolayer(prevssarima$mean, series = "SSARIMA")+
  labs(x = "Temps", y = "Valeur des prévisions", colour = "Méthodes") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
    labels = function(x) format(as.yearmon(x), "%b %Y")
  ) +
  scale_y_continuous(labels = scales::label_number())+
  theme_minimal() +
  scale_color_manual(values = c(
    "Réalisée" = "black",
    "Holt-Winters" = "#1b9e77",
    "ETS" = "#d95f02",
    "TBATS" = "red",
    "ADAM ETS" = "#e7298a",
    "ADAM ETS SARIMA" = "#66a61e",
    "Naïve" = "#e6ab02",
    "STL" = "#a6761d",
    "STS" = "#666666",
    "X13_ARIMA" = "#1f78b4",
    "BAGGED" = "#b2df8a",
    "SARIMA" = "#fb9a99",
    "SSARIMA" = "cyan"
  ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#########============= Prévisions Multivariées =============#########

#########============= 1 variable explicative =============#########


grippe_google <- data_entrainement_diff$Grippe

grippe_google <- ts(grippe_google, start = c(2004, 2), frequency = 12)

obs_grippe_google <- obs_google$Grippe

## SARIMAX1
model_sarimax1 <- auto.arima(inc_ts, xreg = grippe_google , seasonal = TRUE, stationary = TRUE)
summary(model_sarimax1)

prev_sarimax1 <- forecast(model_sarimax1, xreg = obs_grippe_google, h = 12)

## ARX1
model_arx1=auto.arima(inc_ts, max.q = 0, xreg = grippe_google, seasonal= FALSE, stationary=TRUE)
prev_arx1 <- forecast::forecast(model_arx1, xreg = obs_grippe_google, h = 12)

## LM1

model_lm1 <- lm(inc~ Grippe, data=data_entrainement)

obs_grippe_google <-  as.data.frame(obs_grippe_google) |> 
  rename("Grippe"="obs_grippe_google")

prev_lm1 <- predict(model_lm1, newdata =obs_grippe_google)
prev_lm1 <- ts(prev_lm1, start = c(2024,04), frequency = 12)

### Plot des prévisions 1 variable explicative
autoplot(prev_sarimax1$mean, series = "SARIMAX1")+
  autolayer(prev_arx1$mean, series = "ARX1") +
  autolayer(prev_lm1, series = "LM1")+
  labs(
    x = "Temps", y = "Valeur des prévisions") +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),labels = function(x) format(as.yearmon(x), "%b %Y")) +
  scale_y_continuous(labels = scales::label_number_auto()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#########============= Plusieures variable explicative =============#########

#########============= Séléction des variables =============#########

#LEAPS(BESSS)

leaps <- leaps::regsubsets(inc ~ ., data=base_quanti_diff, nbest=1, method=c("exhaustive"))
summary(leaps)

# Choosing the optimal model
res.sum <- summary(leaps)
data.frame(
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic)
)

par(mfrow=c(2,2))
plot(leaps, scale="adjr2", main = "Adjusted R2")
plot(leaps, scale="Cp", main = "Cp")
plot(leaps, scale="bic", main = "BIC")

#STEPWISE

model_null <- lm(inc~1, data = base_quanti_diff)
summary(model_null)

model_lm_complet <- lm(inc~ .,data = base_quanti_diff)
summary(model_lm_complet)

model_forward <- step(model_null, scope = list(lower = model_null, upper = model_lm_complet), direction = "forward")
summary(model_forward)

model_both <- step(model_null,scope = list(upper = model_lm_complet),direction = "both")
summary(model_both)

model_back <- step(model_lm_complet, direction = "backward")
summary(model_back)


# Base avec les variables sélectionnées

xreg <- data_entrainement_diff[, c("Grippe", "Rhume", "Mal_de_gorge", "Toux", "Fièvre")]
xreg <- as.matrix(xreg)

xreg_obs <- obs_google[,  c("Grippe", "Rhume", "Mal_de_gorge", "Toux", "Fièvre")]
xreg_obs <- as.matrix(xreg_obs)


## SARIMAX2

model_sarimax2 <- auto.arima(inc_ts, xreg = xreg, seasonal = TRUE, stationary = TRUE)
prev_sarimax2 <- forecast(model_sarimax2, xreg = xreg_obs, h = 12)



## ARX2

model_ARX2=auto.arima(inc_ts, max.q = 0, xreg = xreg, seasonal= FALSE,stationary=TRUE)
prev_arx2 <- forecast::forecast(model_ARX2, xreg = xreg_obs, h = 12)



## RLM2
model_lm2 <- lm(inc~ Grippe+ Rhume+ Mal_de_gorge+Toux+Fièvre, data = data_entrainement_diff)
new_data <- as.data.frame(xreg_obs)
prev_lm2 <- predict(model_lm2, newdata = new_data)
prev_lm2 <- ts(prev_lm2, start = c(2024,04), frequency = 12)


## NNETAR

set.seed(123)
fit_nnetar <- nnetar(inc_ts, size = 6,repeats = 100, lambda = 0, xreg=xreg)

### Graphe des prévisions plusieures variables

autoplot(prev_sarimax2$mean, series = "SARIMAX2")+
  autolayer(prev_arx2$mean, series = "ARX2") +
  autolayer(prev_lm2, series = "LM2")+
  labs(
    x = "Temps", y = "Valeur des prévisions") +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),labels = function(x) format(as.yearmon(x), "%b %Y")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### NNETAR plot
autoplot(prev_nnetar$mean, series = "NNETAR")+
  labs(
    x = "Temps", y = "Valeur des prévisions") +
  scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),labels = function(x) format(as.yearmon(x), "%b %Y")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#########============= Evaluation des prévisions =============#########

#plot univariés
autoplot(zz, series = "Réalisée",linetype = "dashed" ) +
  autolayer(prev_HW$mean, series = "Holt-Winters", linetype = "dashed") +
  autolayer(prevets$mean, series = "ETS") +
  autolayer(prevtbats$mean, series = "TBATS") +
  autolayer(prevadam1$mean, series = "ADAM ETS") +
  autolayer(prevadam_sarima$mean, series = "ADAM ETS SARIMA") +
  autolayer(prev_naive$mean, series = "Naïve") +
  autolayer(prevstl$mean, series = "STL") +
  autolayer(prevsts$mean, series = "STS") +
  autolayer(prev_x13_arima_ts, series = "X13_ARIMA") +
  autolayer(prevbag$mean, series = "BAGGED") +
  autolayer(prevsarima$mean, series = "SARIMA") +
  autolayer(prevssarima$mean, series = "SSARIMA")+
  labs(x = "Temps", y = "Valeur des prévisions", colour = "Méthodes") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
    labels = function(x) format(as.yearmon(x), "%b %Y")
  ) +
  scale_y_continuous(labels = scales::label_number())+
  theme_minimal() +
  scale_color_manual(values = c(
  "Réalisée" = "black",
  "Holt-Winters" = "#1b9e77",
  "ETS" = "#d95f02",
  "TBATS" = "red",
  "ADAM ETS" = "#e7298a",
  "ADAM ETS SARIMA" = "#66a61e",
  "Naïve" = "#e6ab02",
  "STL" = "#a6761d",
  "STS" = "#666666",
  "X13_ARIMA" = "#1f78b4",
  "BAGGED" = "#b2df8a",
  "SARIMA" = "#fb9a99",
  "SSARIMA" = "cyan"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#plot multivariés
autoplot(zz, series = "Réalisée",linetype = "dashed")+
  autolayer(prev_nnetar$mean, series = "NNETAR")+
  autolayer(prev_sarimax2$mean, series = "SARIMAX2")+
  autolayer(prev_arx2$mean, series = "ARX2")+
  autolayer(prev_lm2, series = "LM2")+
  autolayer(prev_sarimax1$mean, series = "SARIMAX1")+
  autolayer(prev_arx1$mean, series = "ARX1")+
  autolayer(prev_lm1, series = "LM1")+
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
    labels = function(x) format(as.yearmon(x), "%b %Y")
  ) +
  labs(x = "Temps", y = "Valeur des prévisions") +
  scale_y_continuous(labels = scales::label_number())+
  theme_minimal() +
  scale_color_manual(values = c(
    "Réalisée" = "black",
    "NNETAR" = "#1b9e77",
    "SARIMAX2" = "darkblue",
    "SARIMAX1" = "red",
    "ARX1" = "cyan",
    "ARX2" = "#66a61e",
    "LM1" = "#e6ab02",
    "LM2" = "#666666"
  ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#Métriques de performances

#Fonction accuracy() pour tous les modèles --> MAE MAPE RMSE

accuracy(prev_arx1, y_obs) #Exemple


#########============= Comparaison entre les meilleurs modèles =============#########


#R2_OOS

mse_meilleurs_mod <- list(
  Naïf = mean((y_obs - as.numeric(prev_naive$mean))^2),
  BAGGED = mean((y_obs - as.numeric(prevbag$mean))^2),
  ADAM_ETS = mean((y_obs - as.numeric(prevadam1$mean))^2),
  X13_ARIMA = mean((y_obs - as.numeric(prev_x13_arima_ts))^2),
  SARIMAX1= mean((y_obs - as.numeric(prev_sarimax1$mean))^2),
  SARIMAX2 = mean((y_obs - as.numeric(prev_sarimax2$mean))^2),
  NNETAR = mean((y_obs - as.numeric(prev_nnetar$mean))^2)
)


#Calcul du R2 pour les modèles multivariés

R2_OOS_meilleurs_mod <- c(
  Naïf = 1 - mse_meilleurs_mod[["Naïf"]] / mse_meilleurs_mod[["Naïf"]],
  SARIMAX1 = 1 - mse_meilleurs_mod[["SARIMAX1"]] / mse_meilleurs_mod[["Naïf"]],
  SARIMAX2 = 1 - mse_meilleurs_mod[["SARIMAX2"]] / mse_meilleurs_mod[["Naïf"]],
  NNETAR = 1 - mse_meilleurs_mod[["NNETAR"]] / mse_meilleurs_mod[["Naïf"]],
  BAGGED = 1 - mse_meilleurs_mod[["BAGGED"]] / mse_meilleurs_mod[["Naïf"]],
  ADAM_ETS = 1 - mse_meilleurs_mod[["ADAM_ETS"]] / mse_meilleurs_mod[["Naïf"]],
  X13_ARIMA = 1 - mse_meilleurs_mod[["X13_ARIMA"]] / mse_meilleurs_mod[["Naïf"]]
)


#DM test

erreurs_naïf <- y_obs - prev_naive$mean
erreurs_x13 <- y_obs - prev_x13_arima_ts
erreurs_bagged <- y_obs - prevbag$mean
erreur_adam_ets <- y_obs - prevadam1$mean
erreurs_sarimax1 <- y_obs - prev_sarimax1$mean
erreurs_sarimax2 <- y_obs - prev_sarimax2$mean
erreurs_nnetar <- y_obs - prev_nnetar$mean

dm.test(erreurs_x13, erreurs_naïf, alternative = "l", power=1) 
dm.test(erreurs_bagged, erreurs_naïf, alternative = "l", power=1)
dm.test(erreur_adam_ets, erreurs_naïf, alternative = "l", power=1)
dm.test(erreurs_sarimax1, erreurs_naïf, alternative = "l", power=1)
dm.test(erreurs_sarimax2, erreurs_naïf, alternative = "l", power=1)
dm.test(erreurs_nnetar, erreurs_naïf, alternative = "l", power=1)

# Plot 
autoplot(zz, series = "Réalisée", linetype = "dashed")+
  autolayer(prev_nnetar$mean, series = "NNETAR")+
  autolayer(prev_sarimax2$mean, series = "SARIMAX2")+
  autolayer(prev_sarimax1$mean, series = "SARIMAX1")+
  autolayer(prev_x13_arima_ts, series = "X13-ARIMA")+
  autolayer(prevadam1$mean, series = "ADAM-ETS")+
  autolayer(prevbag$mean, series = "BAGGED")+
  autolayer(prev_naive$mean, series = "NAÏVE")+
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
    labels = function(x) format(as.yearmon(x), "%b %Y")
  ) +
  labs(x = "Temps", y = "Valeur des prévisions") +
  scale_y_continuous(labels = scales::label_number())+
  theme_minimal() +
  scale_color_manual(values = c(
    "Réalisée" = "black",
    "NNETAR"     = "magenta",
    "SARIMAX2"   = "#08306B",   
    "SARIMAX1"   = "#e41a1c",   
    "X13-ARIMA"  = "#00CED1",   
    "NAÏVE"      = "#006400",    
    "BAGGED"     = "#e6ab02", 
    "ADAM-ETS"   = "#999999"    
  ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# CSPE
cspe_prev_x13arima <- cumsum(erreurs_x13^2)
cspe_prev_naïf <- cumsum(erreurs_naïf^2)
cspe_prev_bagged <- cumsum(erreurs_bagged^2)
cspe_prev_sarimax1 <- cumsum(erreurs_sarimax1^2)
cspe_prev_sarimax2 <- cumsum(erreurs_sarimax2^2)
cspe_prev_nnetar <- cumsum(erreurs_nnetar^2)
cspe_prev_adam_ets <- cumsum(erreurs_adam_ets^2)

##DF
cspe_df <- data.frame(
  Période = time(prev_naive$mean),
  Naïf = cspe_prev_naïf,
  X13_ARIMA = cspe_prev_x13arima,
  BAGGED = cspe_prev_bagged,
  SARIMAX1 = cspe_prev_sarimax1,
  SARIMAX2 = cspe_prev_sarimax2,
  NNETAR = cspe_prev_nnetar,
  ADAM_ETS = cspe_prev_adam_ets
)

cspe_long <- cspe_df |> 
  pivot_longer(-Période, names_to = "Modèle", values_to = "CSPE") |> 
  mutate(Période = lubridate::ym(format(zoo::as.yearmon(Période), "%Y-%m")))

### plot
ggplot(cspe_long, aes(x = Période, y = CSPE, color = Modèle)) +
  geom_line() +
  labs(title = "Erreurs de Prévision Cumulées au Carré (CSPE)",
       x = "Période",
       y = "CSPE",
       color = "Modèle") +
  scale_color_manual(values = c(
    "Naïf" = "blue",
    "NNETAR"     = "#1b9e77",
    "SARIMAX2"   = "#08306B",   
    "SARIMAX1"   = "#e41a1c",   
    "X13_ARIMA"  = "#00CED1",   
    "NAÏVE"      = "#006400",    
    "BAGGED"     = "#e6ab02", 
    "ADAM_ETS"   = "#999999" 
  )) +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


#########============= SARIMAX3 =============#########

#Suppresion de mal_de_gorge
xreg2 <- data_entrainement_diff[, c("Grippe", "Rhume", "Toux", "Fièvre")]
xreg2 <- as.matrix(xreg2)
xreg_obs2 <- obs_google[,  c("Grippe", "Rhume", "Toux", "Fièvre")]
xreg_obs2 <- as.matrix(xreg_obs2)

#Estimation/prévision
model_sarimax3 <- auto.arima(inc_ts, xreg = xreg2, seasonal = TRUE, stationary = TRUE)
prev_sarimax3 <- forecast(model_sarimax3, xreg = xreg_obs2, h = 12)

accuracy(prev_sarimax3, y_obs)

#Plot SARIMAX2, SARIMAX3, Obs
autoplot(zz, series = "Réalisée",linetype = "dashed") +
  autolayer(prev_sarimax2$mean, series = "SARIMAX2") +
  autolayer(prev_sarimax3$mean, series = "SARIMAX3") +
  scale_color_manual(values = c(
    "Réalisée" = "black",
    "SARIMAX2" = "darkblue",
    "SARIMAX3" = "red"
  ))+
  labs(
    x = "Temps", y = "Valeur des prévisions") +
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
    labels = function(x) format(as.yearmon(x), "%b %Y")
  ) +
  labs(x = "Temps", y = "Valeur des prévisions") +
  scale_y_continuous(labels = scales::label_number())+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#R2 oos SARIMAX3
mse_sarimax3 <- mean((y_obs - as.numeric(prev_sarimax3$mean))^2)
R2_sarimax3 = 1 - mse_sarimax3 / mse_meilleurs_mod[["Naïf"]]

#DM test SARIMAX3
erreurs_sarimax3 <- y_obs - prev_sarimax3$mean
dm.test(erreurs_sarimax3, erreurs_naïf, alternative = "l", power=1)

#Residus SARIMAX3
par(mfrow =c(2,2))
hist(model_sarimax3$residuals, main = "Histogramme des résidus du modèle SARIMAX3", xlab = "Résidus", col = "lightblue", border = "black")
Pacf(model_sarimax3$residuals, lwd = 2,col = "black", lag.max = 24, main="ACF des résidus du modèle SARIMAX3")
Acf(model_sarimax3$residuals, lwd = 2,col = "black", lag.max = 24, main="ACF des résidus du modèle SARIMAX3")
plot(model_sarimax3$residuals, type = "l", main = "Résidus du modèle SARIMAX3", ylab = "Résidus", xlab = "Temps")



#########============= Graphique note de synthèse =============#########
autoplot(zz, series = "Réalisée", linetype = "dashed")+
  autolayer(prev_sarimax3$mean, series = "SARIMAX3")+
  autolayer(prev_nnetar$mean, series = "NNETAR")+
  autolayer(prev_sarimax1$mean, series = "SARIMAX1")+
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
    labels = function(x) format(as.yearmon(x), "%b %Y")
  ) +
  labs(x = "Temps", y = "Valeur des prévisions") +
  scale_y_continuous(labels = scales::label_number())+
  theme_classic() +
  scale_color_manual(values = c(
    "Réalisée" = "black",
    "BAGGED"  = "red",
    "SARIMAX3"   = "blue",
    "NNETAR"= "darkgreen",
    "SARIMAX1" = "orange",
    "X13-ARIMA" = "purple"
    
  ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

autoplot(zz, series = "Réalisée", linetype = "dashed")+
  autolayer(prevbag$mean, series = "BAGGED")+
  autolayer(prev_x13_arima_ts, series = "X13-ARIMA")+
  scale_x_continuous(
    breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 1/12),
    labels = function(x) format(as.yearmon(x), "%b %Y")
  ) +
  labs(x = "Temps", y = "Valeur des prévisions") +
  scale_y_continuous(labels = scales::label_number())+
  theme_classic() +
  scale_color_manual(values = c(
    "Réalisée" = "black",
    "BAGGED"  = "red",
    "SARIMAX3"   = "blue",
    "NNETAR"= "darkgreen",
    "SARIMAX1" = "orange",
    "X13-ARIMA" = "purple"
    
  ))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
