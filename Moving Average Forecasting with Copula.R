library("copula")
library("VineCopula")
library("lubridate")
library("xts")
library("stringr")
library("dplyr")
library("zoo")
library("lmtest")
library("vars")
library("car")
library("ggpubr")
library("tseries")
library("fitdistrplus")
library("distr")
library("MonteCarlo")
library("car")
library("MASS")
library("caret")
library('HMMcopula')
library("tibble")
library("dplyr")
library("EnvStats")
library("plotly")
library("TSstudio")
library("urca")
library("tsDyn")
library("MTS")
library("lmtest")
library("magrittr")
#library("VARtests")
#install.packages("VARtests")
source("Copula_Class.R")

# Setting:
setwd("C:/Users/Dr. Hoang/Desktop/Copula Forecasting/Codes and Data/")
options(scipen = 30)

Copula <- Copula_Class$new(Copula_Family = c(1,2,3,4,5,6))

########################################## FUNCTIONS ##########################################
data_edit <- function(index_input, headers_input, headers_kept){
        
        colnames(index_input) <- headers_input
        index_input['Close'] <- na.locf(as.double(unlist(index_input['Close'])))
        index_input <- subset(index_input, select = headers_kept)
        return (index_input)
}

# Function that generates random values from selected Parametric Bivariate Copula:
BiCop_Simulation <- function(n, family, par, par2){
        
        result = BiCopSim(n, family, par, par2)
        return(list('result1' = result[1], 'result2' = result[2]))
        
}

# Monte Carlos Simulation:
Monte_Carlo_Simulation <- function(nreps, forecast_intervals, series){
        
        Marginal_Dist_Calibration_1 = fitdist(series[,1], distr = "logis")
        Marginal_Dist_Calibration_2 = fitdist(series[,2], distr = "logis")
        
        df_sim_value1 = data.frame(dummy_col = c(1:forecast_intervals))
        df_sim_value2 = data.frame(dummy_col = c(1:forecast_intervals))
        
        Copula_Calibration = Copula$Copula_Calibration(series,1,2)
        
        Copula_Object = BiCop(Copula_Calibration$family, Copula_Calibration$par, Copula_Calibration$par2, 
                              Copula_Calibration$tau)
        
        for (i in 1:nreps){
                Copula_Generator = BiCopSim(forecast_intervals, Copula_Calibration$family, Copula_Calibration$par, Copula_Calibration$par2)
                
                Simulated_Values_1 = qlogis(Copula_Generator[,1], location = Marginal_Dist_Calibration_1$estimate['location'], 
                                            scale = Marginal_Dist_Calibration_1$estimate['scale'])
                
                Simulated_Values_2 = qlogis(Copula_Generator[,2], location = Marginal_Dist_Calibration_2$estimate['location'], 
                                            scale = Marginal_Dist_Calibration_2$estimate['scale'])
                
                df_sim_value1 = cbind(df_sim_value1, Simulated_Values_1)
                df_sim_value2 = cbind(df_sim_value2, Simulated_Values_2)
                
        }
                df_sim_value1 <- df_sim_value1[-c(1)]
                df_sim_value2 <- df_sim_value2[-c(1)]
                
                return (list('Simulated_Value_1' = df_sim_value1, 'Simulated_Value_2' = df_sim_value2))
}

# Markov Copula Simulation
Markov_Copula_Sim <- function(paired_series, nreps, forecast_num, num_regimes, copula_family, iter_num, error){
        
        df_sim_value1 = data.frame(dummy_col = c(1:forecast_num))
        df_sim_value2 = data.frame(dummy_col = c(1:forecast_num))
        
        Marginal_Dist_Calibration_1 = fitdist(paired_series[,1], distr = "logis")
        Marginal_Dist_Calibration_2 = fitdist(paired_series[,2], distr = "logis")
        
        Estimated_Copula = EstHMMCop(paired_series, num_regimes, copula_family, iter_num, error)
        
        for (i in 1:nreps){
                
                Copula_Simulation = SimHMMCop(Estimated_Copula$Q, copula_family, Estimated_Copula$tau, forecast_num, Estimated_Copula$dof)
                
                Simulated_Values_1 = qlogis(Copula_Simulation$SimData[,1], location = Marginal_Dist_Calibration_1$estimate['location'], 
                                            scale = Marginal_Dist_Calibration_1$estimate['scale'])
                
                Simulated_Values_2 = qlogis(Copula_Simulation$SimData[,2], location = Marginal_Dist_Calibration_2$estimate['location'], 
                                            scale = Marginal_Dist_Calibration_2$estimate['scale'])
                
                df_sim_value1 = cbind(df_sim_value1, Simulated_Values_1)
                df_sim_value2 = cbind(df_sim_value2, Simulated_Values_2)
        }
        
        df_sim_value1 <- df_sim_value1[-c(1)]
        df_sim_value2 <- df_sim_value2[-c(1)]
        
        return (list('Simulated_Value_1' = df_sim_value1, 'Simulated_Value_2' = df_sim_value2))
}

# Granger Causality test for 2-way interaction of stock/index return:

# H0: X does not  Granger cause Y

# 
# A small p (??? 0.05) => reject the null hypothesis

Independence_Test <- function(df, df_training, granger_order, n = 1){
        for (i in names(df_training)) {
                n = n+1
                if (n > length(names(df_training))) {break}
                else
                        remained_stocks <- names(df_training)[seq(from=n, to=length(names(df_training)), by=1)]
                for (j in remained_stocks){
                        Ind_Test_BiCop = BiCopIndTest(pobs(df_training[i]), pobs(df_training[j]))
                        df[nrow(df)+1,] = c(i, j, grangertest(df_training[j], df_training[i], order = granger_order)[2,"Pr(>F)"], Ind_Test_BiCop$p, TRUE, TRUE)
                        df[nrow(df)+1,] = c(j, i, grangertest(df_training[i], df_training[j], order = granger_order)[2,"Pr(>F)"], Ind_Test_BiCop$p, TRUE, TRUE)
                        
                }
                
        }
        
        df %<>% mutate(
                GRANGER_5pct = case_when(
                        GRANGER_p_value > 0.05 ~ FALSE,
                        GRANGER_p_value <= 0.05 ~ TRUE
                ),
                Ind_Test_5pct = case_when(
                        Ind_Test_p_value > 0.05 ~ FALSE,
                        Ind_Test_p_value <= 0.05 ~ TRUE
                )
        )
        
        return(df)
}

########################################## FUNCTIONS ##########################################

headers <- c("Date", "Open", "High", "Low", "Close", "Volume")
headers_select <- c("Date", "Close")

# Load indexes' historical data:
AAPL <- data_edit(read.csv(file = 'AAPL.csv'), headers, headers_select)
GOOG <- data_edit(read.csv(file = 'GOOG.csv'), headers, headers_select)
MSFT <- data_edit(read.csv(file = 'MSFT.csv'), headers, headers_select)
NVDA <- data_edit(read.csv(file = 'NVDA.csv'), headers, headers_select)
NASDAQ <- data_edit(read.csv(file = 'NASDAQ_100.csv'), headers, headers_select)
SP500 <- data_edit(read.csv(file = 'SP500.csv'), headers, headers_select)

# Extract Adj.Price column from each stock:

AAPL_Price = as.data.frame(AAPL[2])
GOOG_Price = as.data.frame(GOOG[2])
MSFT_Price = as.data.frame(MSFT[2])
NVDA_Price = as.data.frame(NVDA[2])
NASDAQ_Price = as.data.frame(NASDAQ[2])
SP500_Price = as.data.frame(SP500[2])

colnames(AAPL_Price)[1] <- "AAPL"
colnames(GOOG_Price)[1] <- "GOOG"
colnames(MSFT_Price)[1] <- "MSFT"
colnames(NVDA_Price)[1] <- "NVDA"
colnames(NASDAQ_Price)[1] <- "NASDAQ"
colnames(SP500_Price)[1] <- "SP500"

df_Price = data.frame(AAPL_Price, GOOG_Price, MSFT_Price, NVDA_Price, NASDAQ_Price, SP500_Price)
rownames(df_Price) <- as.POSIXct(AAPL$Date, format = '%m/%d/%Y')

df_Price_training = df_Price[3000:4000,]
df_Price_forecast = df_Price[4001:4431,]
# fcst_period = [5, 10, 20]
# training_bgn_value = 3000
# training_end_value = 4000
# for i in in fcst_period:
        # df_Price_training = df_Price[training_bgn_value + i:training_end_value + i,]
        # df_Price_forecast = df_Price[training_end_value + i + 1:, training_end_value + 2*i]
        # training_bgn_value = training_bgn_value + i
        # training_end_value = training_end_value + i
        # Marginal Distribution calibration
        # Copula calibration
        # Construct VAR[i]

# Calculate log_return of each stock/index:
AAPL_Return <- as.data.frame(sapply(AAPL[2], function(x) diff(log(x))))
GOOG_Return <- as.data.frame(sapply(GOOG[2], function(x) diff(log(x))))
MSFT_Return <- as.data.frame(sapply(MSFT[2], function(x) diff(log(x))))
NVDA_Return <- as.data.frame(sapply(NVDA[2], function(x) diff(log(x))))
NASDAQ_Return <- as.data.frame(sapply(NASDAQ[2], function(x) diff(log(x))))
SP500_Return <- as.data.frame(sapply(SP500[2], function(x) diff(log(x))))

# Rename a column to log_return:
colnames(AAPL_Return)[1] <- "AAPL"
colnames(GOOG_Return)[1] <- "GOOG"
colnames(MSFT_Return)[1] <- "MSFT"
colnames(NVDA_Return)[1] <- "NVDA"
colnames(NASDAQ_Return)[1] <- "NASDAQ"
colnames(SP500_Return)[1] <- "SP500"

# Create a data frame with all stocks' log_returns:
df <- data.frame(AAPL_Return, GOOG_Return, MSFT_Return, NVDA_Return, NASDAQ_Return, SP500_Return)
rownames(df) <- as.POSIXct(AAPL[-1,]$Date, format = '%m/%d/%Y')

# Creating training and forecasting population:
df_training = df[3000:4000,]
 
# Granger Causality test for 2-way interaction of stock/index return:
# Ho: X does not  Granger cause Y

# Testing asymptotic independence in bivariate extreme
# Ho: 2 variables are asymptotically independent

# A small p (??? 0.05) => reject the null hypothesis

df_test <- data.frame(Y = character(),
                      X=character(), 
                      GRANGER_p_value=numeric(), 
                      Ind_Test_p_value=numeric(),
                      GRANGER_5pct=logical(),
                      Ind_Test_5pct=logical())

indepentdence_test = Independence_Test(df_test, df_Price_training,3)

# AAPL_MSFT, AMZN_NASDAQ, AMZN_SP500, MSFT_SP500, NVDA_SP500

# Pre-scanning for potential marginal distribution:
descdist(df_training[,1], discrete = FALSE) # check for possible marginal distribution for AAPL Return series
descdist(df_training[,2], discrete = FALSE) # check for possible marginal distribution for GOOG Return series
descdist(df_training[,3], discrete = FALSE) # check for possible marginal distribution for MSFT Return series
descdist(df_training[,4], discrete = FALSE) # check for possible marginal distribution for NVDA Return series
descdist(df_training[,5], discrete = FALSE) # check for possible marginal distribution for NASDAQ Return series
descdist(df_training[,6], discrete = FALSE) # check for possible marginal distribution for SP500 Return series

# Marginal Distribution of each return:
AAPL_Distribution_Logis = fitdist(df_training[,1], distr = "logis") 
AAPL_Distribution_Normal = fitdist(df_training[,1], distr = "norm")

GOOG_Distribution_Logis = fitdist(df_training[,2], distr = "logis")
GOOG_Distribution_Normal = fitdist(df_training[,2], distr = "norm")

MSFT_Distribution_Logis = fitdist(df_training[,3], distr = "logis")
MSFT_Distribution_Normal = fitdist(df_training[,3], distr = "norm")

NVDA_Distribution_Logis = fitdist(df_training[,4], distr = "logis")
NVDA_Distribution_Normal = fitdist(df_training[,4], distr = "norm")

NASDAQ_Distribution_Logis = fitdist(df_training[,5], distr = "logis")
NASDAQ_Distribution_Normal = fitdist(df_training[,5], distr = "norm")

SP500_Distribution_Logis = fitdist(df_training[,6], distr = "logis")
SP500_Distribution_Normal = fitdist(df_training[,6], distr = "norm")

# visual check for fitness via parametric CDF vs. empirical CDF:
cdfcomp(list(AAPL_Distribution_Logis, AAPL_Distribution_Normal), legendtext = c("logistic", "norm"), main = "AAPL's return fits")
cdfcomp(list(GOOG_Distribution_Logis, GOOG_Distribution_Normal), legendtext = c("logistic", "norm"), main = "GOOG's return fits")
cdfcomp(list(MSFT_Distribution_Logis, MSFT_Distribution_Normal), legendtext = c("logistic", "norm"), main = "MSFT's return fits")
cdfcomp(list(NVDA_Distribution_Logis, NVDA_Distribution_Normal), legendtext = c("logistic", "norm"), main = "NVDA's return fits")
cdfcomp(list(NASDAQ_Distribution_Logis, NASDAQ_Distribution_Normal), legendtext = c("logistic", "norm"), main = "NASDAQ's return fits")
cdfcomp(list(SP500_Distribution_Logis, SP500_Distribution_Normal), legendtext = c("logistic", "norm"), main = "SP500's return fits")

# Goodness of fit statistics of marginal distributions:
AAPL_GoF = gofstat(list(AAPL_Distribution_Logis, AAPL_Distribution_Normal), 
                   fitnames = c("logistic", "normal"))
typeof(AAPL_GoF)
GOOG_GoF = gofstat(list(GOOG_Distribution_Logis, GOOG_Distribution_Normal), 
                   fitnames = c("logistic", "normal"))

MSFT_GoF = gofstat(list(MSFT_Distribution_Logis, MSFT_Distribution_Normal), 
                   fitnames = c("logistic", "normal"))

NVDA_GoF = gofstat(list(NVDA_Distribution_Logis, NVDA_Distribution_Normal), 
                   fitnames = c("logistic", "normal"))

NASDAQ_GoF = gofstat(list(NASDAQ_Distribution_Logis, NASDAQ_Distribution_Normal), 
                     fitnames = c("logistic", "normal"))

SP500_GoF = gofstat(list(SP500_Distribution_Logis, SP500_Distribution_Normal), 
                    fitnames = c("logistic", "normal"))

# Generating pair pseudo observation:
AAPL_MSFT_pobs = pobs(df_training[,c("AAPL","MSFT")])
AAPL_NVDA_pobs = pobs(df_training[,c("AAPL","NVDA")])
AAPL_NASDAQ_pobs = pobs(df_training[,c("AAPL","NASDAQ")])
AAPL_SP500_pobs = pobs(df_training[,c("AAPL","SP500")])

GOOG_MSFT_pobs = pobs(df_training[,c("GOOG","MSFT")])
GOOG_NVDA_pobs = pobs(df_training[,c("GOOG","NVDA")])
GOOG_NASDAQ_pobs = pobs(df_training[,c("GOOG","NASDAQ")])

MSFT_SP500_pobs = pobs(df_training[,c("MSFT","SP500")])
MSFT_NASDAQ_pobs = pobs(df_training[,c("MSFT","NASDAQ")])

NVDA_NASDAQ_pobs = pobs(df_training[,c("NVDA","NASDAQ")])
NVDA_SP500_pobs = pobs(df_training[,c("NVDA","SP500")])

NASDAQ_SP500_pobs = pobs(df_training[,c("NASDAQ","SP500")])

# Best fitted Copulas:
AAPL_MSFT_Best_Copula = BiCopSelect(AAPL_MSFT_pobs[,1], AAPL_MSFT_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                    rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")
AAPL_NVDA_Best_Copula = BiCopSelect(AAPL_NVDA_pobs[,1], AAPL_NVDA_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                    rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")
AAPL_NASDAQ_Best_Copula = BiCopSelect(AAPL_NASDAQ_pobs[,1], AAPL_NASDAQ_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                    rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")
AAPL_SP500_Best_Copula = BiCopSelect(AAPL_SP500_pobs[,1], AAPL_SP500_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                      rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")

GOOG_MSFT_Best_Copula = BiCopSelect(GOOG_MSFT_pobs[,1], GOOG_MSFT_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                     rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")
GOOG_NVDA_Best_Copula = BiCopSelect(GOOG_NVDA_pobs[,1], GOOG_NVDA_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                    rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")
GOOG_NASDAQ_Best_Copula = BiCopSelect(GOOG_NASDAQ_pobs[,1], GOOG_NASDAQ_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                    rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")

MSFT_SP500_Best_Copula = BiCopSelect(MSFT_SP500_pobs[,1], MSFT_SP500_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                     rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")
MSFT_NASDAQ_Best_Copula = BiCopSelect(MSFT_NASDAQ_pobs[,1], MSFT_NASDAQ_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                     rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")

NVDA_NASDAQ_Best_Copula = BiCopSelect(NVDA_NASDAQ_pobs[,1], NVDA_NASDAQ_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                     rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")
NVDA_SP500_Best_Copula = BiCopSelect(NVDA_SP500_pobs[,1], NVDA_SP500_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                     rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")

NASDAQ_SP500_Best_Copula = BiCopSelect(NASDAQ_SP500_pobs[,1], NASDAQ_SP500_pobs[,2], familyset =  c(1, 2, 3, 4, 5, 6), selectioncrit = "AIC", 
                                     rotations = FALSE, indeptest = TRUE, level = 0.05, method = "mle")

# Goodness of fit of best copulas:
AAPL_MSFT_Best_Copula_GOF = BiCopGofTest(AAPL_MSFT_pobs[,1], AAPL_MSFT_pobs[,2], obj = AAPL_MSFT_Best_Copula, B = 100)
MSFT_SP500_Best_Copula_GOF = BiCopGofTest(MSFT_SP500_pobs[,1], MSFT_SP500_pobs[,2], obj = MSFT_SP500_Best_Copula, B = 100)
NVDA_SP500_Best_Copula_GOF = BiCopGofTest(NVDA_SP500_pobs[,1], NVDA_SP500_pobs[,2], obj = NVDA_SP500_Best_Copula, B = 100)

# Forecast paired return by calibrating for the best copula:
AAPL_MSFT_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("AAPL","MSFT")])
AAPL_NVDA_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("AAPL","NVDA")])
AAPL_NASDAQ_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("AAPL","NASDAQ")])
AAPL_SP500_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("AAPL","SP500")])

GOOG_MSFT_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("GOOG","MSFT")])
GOOG_NVDA_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("GOOG","NVDA")])
GOOG_NASDAQ_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("GOOG","NASDAQ")])

MSFT_NASDAQ_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("MSFT","NASDAQ")])
MSFT_SP500_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("MSFT","SP500")])

NVDA_NASDAQ_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("NVDA","NASDAQ")])
NVDA_SP500_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("NVDA","SP500")])

NASDAQ_SP500_Simulation <- Monte_Carlo_Simulation(10000, 431, df_training[,c("NASDAQ","SP500")])

# Regime Switching Approach for copula forecasting:
AAPL_MSFT_Markov_2 = Markov_Copula_Sim(df_training[,c("AAPL","MSFT")], 1000, 431, 2, "t", 100, 0.001)
AAPL_MSFT_Markov_3 = Markov_Copula_Sim(df_training[,c("AAPL","MSFT")], 1000, 431, 3, "t", 100, 0.001)

AAPL_NVDA_Markov_2 = Markov_Copula_Sim(df_training[,c("AAPL","NVDA")], 1000, 431, 2, "t", 100, 0.001)
AAPL_NVDA_Markov_3 = Markov_Copula_Sim(df_training[,c("AAPL","NVDA")], 1000, 431, 3, "t", 100, 0.001)

AAPL_NASDAQ_Markov_2 = Markov_Copula_Sim(df_training[,c("AAPL","NASDAQ")], 1000, 431, 2, "t", 100, 0.001)
AAPL_NASDAQ_Markov_3 = Markov_Copula_Sim(df_training[,c("AAPL","NASDAQ")], 1000, 431, 3, "t", 100, 0.001)

AAPL_SP500_Markov_2 = Markov_Copula_Sim(df_training[,c("AAPL","SP500")], 1000, 431, 2, "t", 100, 0.001)
AAPL_SP500_Markov_3 = Markov_Copula_Sim(df_training[,c("AAPL","SP500")], 1000, 431, 3, "t", 100, 0.001)

MSFT_SP500_Markov_2 = Markov_Copula_Sim(df_training[,c("MSFT","SP500")], 1000, 431, 2, "t", 100, 0.001)
MSFT_SP500_Markov_3 = Markov_Copula_Sim(df_training[,c("MSFT","SP500")], 1000, 431, 3, "t", 100, 0.001)

MSFT_NASDAQ_Markov_2 = Markov_Copula_Sim(df_training[,c("MSFT","NASDAQ")], 1000, 431, 2, "t", 100, 0.001)
MSFT_NASDAQ_Markov_3 = Markov_Copula_Sim(df_training[,c("MSFT","NASDAQ")], 1000, 431, 3, "t", 100, 0.001)

NVDA_SP500_Markov_2 = Markov_Copula_Sim(df_training[,c("NVDA","SP500")], 1000, 431, 2, "t", 100, 0.001)
NVDA_SP500_Markov_3 = Markov_Copula_Sim(df_training[,c("NVDA","SP500")], 1000, 431, 3, "t", 100, 0.001)

NVDA_NASDAQ_Markov_2 = Markov_Copula_Sim(df_training[,c("NVDA","NASDAQ")], 1000, 431, 2, "t", 100, 0.001)
NVDA_NASDAQ_Markov_3 = Markov_Copula_Sim(df_training[,c("NVDA","NASDAQ")], 1000, 431, 3, "t", 100, 0.001)

NASDAQ_SP500_Markov_2 = Markov_Copula_Sim(df_training[,c("NVDA","SP500")], 1000, 431, 2, "t", 100, 0.001)
NASDAQ_SP500_Markov_3 = Markov_Copula_Sim(df_training[,c("NVDA","SP500")], 1000, 431, 3, "t", 100, 0.001)

#AAPL_MSFT:
VAR_AAPL_MSFT <- vars::VAR(as.matrix(df_Price_training[,c("AAPL","MSFT")]), p = 10)

forecast_AAPL_MSFT = predict(VAR_AAPL_MSFT, n.ahead = 431, ci = 0.95)

johan.test.AAPL_MSFT <- ca.jo(df_Price_training[,c("AAPL","MSFT")],
                          ecdet = "const",
                          type = "trace",
                          K = 10,
                          ) 
#summary(johan.test.AAPL_MSFT) 

AAPL_MSFT_forecast <- predict(vec2var(johan.test.AAPL_MSFT, r = 1),
                             n.ahead = 431,
                             ci = 0.95)

#AAPL_NVDA:
VAR_AAPL_NVDA <- vars::VAR(as.matrix(df_Price_training[,c("AAPL","NVDA")]), p = 10)

forecast_AAPL_NVDA = predict(VAR_AAPL_NVDA, n.ahead = 431, ci = 0.95)

johan.test.AAPL_NVDA <- ca.jo(df_Price_training[,c("AAPL","NVDA")],
                              ecdet = "const",
                              type = "trace",
                              K = 10,
) 
#summary(johan.test.AAPL_NVDA) 

AAPL_NVDA_forecast <- predict(vec2var(johan.test.AAPL_NVDA, r = 1),
                              n.ahead = 431,
                              ci = 0.95)

#AAPL_NASDAQ:
VAR_AAPL_NASDAQ <- vars::VAR(as.matrix(df_Price_training[,c("AAPL","NASDAQ")]), p = 10)
forecast_AAPL_NASDAQ = predict(VAR_AAPL_NASDAQ, n.ahead = 431, ci = 0.95)

johan.test.AAPL_NASDAQ <- ca.jo(df_Price_training[,c("AAPL","NASDAQ")],
                              ecdet = "const",
                              type = "trace",
                              K = 10,
) 
#summary(johan.test.AAPL_NASDAQ) 

AAPL_NASDAQ_forecast <- predict(vec2var(johan.test.AAPL_NASDAQ, r = 1),
                              n.ahead = 431,
                              ci = 0.95)

#AAPL_SP500:
VAR_AAPL_SP500 <- vars::VAR(as.matrix(df_Price_training[,c("AAPL","SP500")]), p = 10)
forecast_AAPL_SP500 = predict(VAR_AAPL_SP500, n.ahead = 431, ci = 0.95)

johan.test.AAPL_SP500 <- ca.jo(df_Price_training[,c("AAPL","SP500")],
                                ecdet = "const",
                                type = "trace",
                                K = 10,
) 
#summary(johan.test.AAPL_SP500) 

AAPL_SP500_forecast <- predict(vec2var(johan.test.AAPL_SP500, r = 1),
                              n.ahead = 431,
                              ci = 0.95)

# GOOG_MSFT:
VAR_GOOG_MSFT <- vars::VAR(as.matrix(df_Price_training[,c("GOOG","MSFT")]), p = 10)
forecast_GOOG_MSFT = predict(VAR_GOOG_MSFT, n.ahead = 431, ci = 0.95)

johan.test.GOOG_MSFT <- ca.jo(df_Price_training[,c("GOOG","MSFT")],
                               ecdet = "const",
                               type = "trace",
                               K = 10,
) 
#summary(johan.test.GOOG_MSFT) 

GOOG_MSFT_forecast <- predict(vec2var(johan.test.GOOG_MSFT, r = 1),
                               n.ahead = 431,
                               ci = 0.95)

# GOOG_NVDA:
VAR_GOOG_NVDA <- vars::VAR(as.matrix(df_Price_training[,c("GOOG","NVDA")]), p = 10)
forecast_GOOG_NVDA = predict(VAR_GOOG_NVDA, n.ahead = 431, ci = 0.95)

johan.test.GOOG_NVDA <- ca.jo(df_Price_training[,c("GOOG","NVDA")],
                              ecdet = "const",
                              type = "trace",
                              K = 10,
) 
#summary(johan.test.GOOG_NVDA) 

GOOG_NVDA_forecast <- predict(vec2var(johan.test.GOOG_NVDA, r = 1),
                              n.ahead = 431,
                              ci = 0.95)

# GOOG_NASDAQ:
VAR_GOOG_NASDAQ <- vars::VAR(as.matrix(df_Price_training[,c("GOOG","NASDAQ")]), p = 10)
forecast_GOOG_NASDAQ = predict(VAR_GOOG_NASDAQ, n.ahead = 431, ci = 0.95)

johan.test.GOOG_NASDAQ <- ca.jo(df_Price_training[,c("GOOG","NASDAQ")],
                              ecdet = "const",
                              type = "trace",
                              K = 10,
) 
#summary(johan.test.GOOG_NASDAQ) 

GOOG_NASDAQ_forecast <- predict(vec2var(johan.test.GOOG_NASDAQ, r = 1),
                              n.ahead = 431,
                              ci = 0.95)

# MSFT_SP500:
VAR_MSFT_SP500 <- vars::VAR(as.matrix(df_Price_training[,c("MSFT","SP500")]), p = 10)
forecast_MSFT_SP500 = predict(VAR_MSFT_SP500, n.ahead = 431, ci = 0.95)

johan.test.MSFT_SP500 <- ca.jo(df_Price_training[,c("MSFT","SP500")],
                                ecdet = "const",
                                type = "trace",
                                K = 10,
) 
#summary(johan.test.MSFT_SP500) 

MSFT_SP500_forecast <- predict(vec2var(johan.test.MSFT_SP500, r = 1),
                                n.ahead = 431,
                                ci = 0.95)

# MSFT_NASDAQ:
VAR_MSFT_NASDAQ <- vars::VAR(as.matrix(df_Price_training[,c("MSFT","NASDAQ")]), p = 10)
forecast_MSFT_NASDAQ = predict(VAR_MSFT_NASDAQ, n.ahead = 431, ci = 0.95)

johan.test.MSFT_NASDAQ <- ca.jo(df_Price_training[,c("MSFT","NASDAQ")],
                               ecdet = "const",
                               type = "trace",
                               K = 10,
) 
#summary(johan.test.MSFT_NASDAQ) 

MSFT_NASDAQ_forecast <- predict(vec2var(johan.test.MSFT_NASDAQ, r = 1),
                               n.ahead = 431,
                               ci = 0.95)

# NVDA_SP500:
VAR_NVDA_SP500 <- vars::VAR(as.matrix(df_Price_training[,c("NVDA","SP500")]), p = 10)
forecast_NVDA_SP500 = predict(VAR_NVDA_SP500, n.ahead = 431, ci = 0.95)

johan.test.NVDA_SP500 <- ca.jo(df_Price_training[,c("NVDA","SP500")],
                               ecdet = "const",
                               type = "trace",
                               K = 10,
) 
#summary(johan.test.NVDA_SP500)

NVDA_SP500_forecast <- predict(vec2var(johan.test.NVDA_SP500, r = 1),
                                n.ahead = 431,
                                ci = 0.95)

# NVDA_NASDAQ:
VAR_NVDA_NASDAQ <- vars::VAR(as.matrix(df_Price_training[,c("NVDA","NASDAQ")]), p = 10)
forecast_NVDA_NASDAQ = predict(VAR_NVDA_NASDAQ, n.ahead = 431, ci = 0.95)

plot(resid(VAR_NVDA_NASDAQ)[,1])
plot(resid(VAR_NVDA_NASDAQ)[,2])

johan.test.NVDA_NASDAQ <- ca.jo(df_Price_training[,c("NVDA","NASDAQ")],
                               ecdet = "const",
                               type = "trace",
                               K = 10,
) 
#summary(johan.test.NVDA_NASDAQ)

NVDA_NASDAQ_forecast <- predict(vec2var(johan.test.NVDA_NASDAQ, r = 1),
                               n.ahead = 431,
                               ci = 0.95)

# NASDAQ_SP500:
VAR_NASDAQ_SP500 <- vars::VAR(as.matrix(df_Price_training[,c("NASDAQ","SP500")]), p = 10, type = 'const')

plot(resid(VAR_NASDAQ_SP500)[,1])
plot(resid(VAR_NASDAQ_SP500)[,2])

descdist(resid(VAR_NASDAQ_SP500)[,1], discrete = FALSE)
NASDAQ_Resi_Distribution_Logis = fitdist(resid(VAR_NASDAQ_SP500)[,2], distr = "logis") 
NASDAQ_Resi_Distribution_Normal = fitdist(resid(VAR_NASDAQ_SP500)[,2], distr = "norm")
NASDAQ_Resi_Distribution_Cauchy = fitdist(resid(VAR_NASDAQ_SP500)[,2], distr = "cauchy")
cdfcomp(list(NASDAQ_Resi_Distribution_Logis, NASDAQ_Resi_Distribution_Normal, NASDAQ_Resi_Distribution_Cauchy), legendtext = c("logistic", "norm", "cauchy"), main = "NASDAQ's return fits")
NASDAQ_GoF = gofstat(list(NASDAQ_Resi_Distribution_Logis, NASDAQ_Resi_Distribution_Normal, NASDAQ_Resi_Distribution_Cauchy), 
                   fitnames = c("logistic", "normal", "cauchy"))

VAR_NASDAQ_SP500_pobs = pobs(resid(VAR_NASDAQ_SP500))

epdfPlot(resid(VAR_NASDAQ_SP500)[,2])

VAR_NASDAQ_SP500_Residual_Copula = BiCopSelect(VAR_NASDAQ_SP500_pobs[,1], VAR_NASDAQ_SP500_pobs[,2], selectioncrit = "AIC", 
                                    rotations = TRUE, indeptest = TRUE, level = 0.05, method = "mle")

NASDAQ_SP500_Resi_Cop_Sim = BiCopSim(431, VAR_NASDAQ_SP500_Residual_Copula$family, VAR_NASDAQ_SP500_Residual_Copula$par, VAR_NASDAQ_SP500_Residual_Copula$par2)

VAR_NASDAQ_SP500_Residual_Copula_fore = predict(VAR_NASDAQ_SP500_Residual_Copula, n.ahead = 431, ci = 0.95)

Simulated_Values_2 = qlogis(NASDAQ_SP500_Resi_Cop_Sim[,2], location = NASDAQ_Resi_Distribution_Logis$estimate['location'], 
                            scale = NASDAQ_Resi_Distribution_Logis$estimate['scale'])

descdist(Simulated_Values_2, discrete = FALSE)
plot(Simulated_Values_2)
epdfPlot(Simulated_Values_2)

johan.test.NASDAQ_SP500 <- ca.jo(df_Price_training[,c("NASDAQ","SP500")],
                                ecdet = "const",
                                type = "trace",
                                K = 10,
) 
#summary(johan.test.NASDAQ_SP500)

NASDAQ_SP500_forecast <- predict(vec2var(johan.test.NASDAQ_SP500, r = 1),
                                n.ahead = 431,
                                ci = 0.95)

# VAR_ALL
VAR_ALL <- vars::VAR(as.matrix(df_Price_training[,c("AAPL", "GOOG", "MSFT", "NVDA", "NASDAQ", "SP500")]), p = 10)
forecast_ALL = predict(VAR_ALL, n.ahead = 431, ci = 0.95)

johan.test.VAR_ALL <- ca.jo(df_Price_training[,c("AAPL", "GOOG", "MSFT", "NVDA", "NASDAQ", "SP500")],
                                 ecdet = "const",
                                 type = "trace",
                                 K = 10,
)


# Visual check of forecast result
Price1 = df_Price[4000,1]
Price2 = df_Price[4000,1]
Price3 = df_Price[4000,1]
#MSFT_Price4 = df_Price[4000,4]
#MSFT_Price5 = df_Price[4000,4]
#SP500_Price_2 = df_Price[4000,1]
#SP500_Price_3 = df_Price[4000,1]

for (i in 1:431){
        Price1 = append(Price1, tail(Price1,1)*(1+rowMeans(AAPL_MSFT_Simulation$Simulated_Value_1[i,])))
}

for (i in 1:431){
        Price2 = append(Price2, tail(Price2,1)*(1+rowMeans(AAPL_MSFT_Markov_2$Simulated_Value_1[i,])))
}

for (i in 1:431){
        Price3 = append(Price3, tail(Price3,1)*(1+rowMeans(AAPL_MSFT_Markov_3$Simulated_Value_1[i,])))
}

#for (i in 1:431){
#        MSFT_Price4 = append(MSFT_Price4, tail(MSFT_Price4,1)*(1+rowMeans(MSFT_SP500_Simulation$Simulated_Value_1[i,])))
#}

#for (i in 1:431){
#        MSFT_Price5 = append(MSFT_Price5, tail(MSFT_Price5,1)*(1+rowMeans(MSFT_NASDAQ_Simulation$Simulated_Value_1[i,])))
#}

#for (i in 1:431){
#        AAPL_Price_3 = append(AAPL_Price_3, tail(AAPL_Price_3,1)*(1+rowMeans(AAPL_SP500_Markov_3$Simulated_Value_1[i,])))
#}

plot(x= 1:4432, y=df_Price[['AAPL']], col="black", lwd=2, type = "l")
lines(x = 4000:4431, y = Price1, col="red", lwd=2, type = "l") # Copula single regime
lines(x = 4000:4431, y = Price2, col="yellow", lwd=2, type = "l") # Copula 2 regimes
lines(x = 4000:4431, y = Price3, col="grey", lwd=2, type = "l") # Copula 3 regimes
lines(x = 4000:4430, y = forecast_AAPL_MSFT$fcst$AAPL[,1], col="green", lwd=2, type = "l") # VAR
lines(x = 4000:4430, y = AAPL_MSFT_forecast$fcst$AAPL[,1], col="blue", lwd=2, type = "l") # VECM

AAPL_MSFT_forecast

#, ylim=c(0,5000)

#[1] AAPL
#[2] GOOG
#[3] MSFT
#[4] NVDA
#[5] NASDAQ
#[6] SP500
