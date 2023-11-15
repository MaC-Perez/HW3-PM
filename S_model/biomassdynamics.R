## surplus production model

## Load TMB TRUE
library(TMB)
library(tidyverse)
library(readxl)
library(TMBhelper)

## Data and parameters

catch <- read_xlsx("data/flounder_index_data.xlsx","catch")
index <- read_xlsx("data/flounder_index_data.xlsx","survey")

index_years <- index$Year - min(catch$Year)

# we have time series of catches and survey index 
data <- list(catches = catch$Catch,
               index = index$Index,
               index_years = index_years)
data

# carrying capacity and growth rate positive
# m shape param
# observation error for catch and survey positive both
#exploitation rate u for every catch year 
parameters <- list(logK=10, 
                   logr=log(0.4), 
                   m=2, 
                   logSigmaC=log(0.05), 
                   logSigmaI=log(0.2),
                   #logq = -2, 
                   upar = rep(-1.6,length(data$catches)))
parameters

# to create de template that we already have
## Make C++ file
#TMB::template("biomassdynamics.cpp")
setwd("C:/Users/macristina.perez/Documents/UMassD/Classes/2023/Fall_2023/Population Modeling/HW/HW3/Perez/S_model")

## Compile and load the model
dll_name <- "biomassdynamics"
if(is.loaded(dll_name)) {
  dyn.unload(dynlib(dll_name))
}
compile("biomassdynamics.cpp")
dyn.load(dynlib(dll_name))


## Make a function object //map arg fixing params in their initial values, we dont want to estimate m and logSigmaC
obj <- MakeADFun(data, parameters, DLL= dll_name,
                 map = list(m=factor(NA),logSigmaC=factor(NA)),
                 control=list(eval.max=10000,iter.max=10000,rel.tol=1e-15))

## Call function minimizer
opt <- nlminb(obj$par, obj$fn, obj$gr, control= list(eval.max = 10000,
                                                        iter.max = 10000))
#goes down to 0 ... ok

## Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
sdr
summary(sdr)

pl <- obj$env$parList(opt$par) 

q <- summary(sdr, type = "report") %>% 
  as.data.frame() %>% 
  janitor::clean_names() %>% 
  rownames_to_column(var="type") %>% 
  filter(str_detect(type, "logq"))

summary(sdr, type = "report") %>% 
  as.data.frame() %>% 
  janitor::clean_names() %>% 
  rownames_to_column(var="type") %>% 
  filter(str_detect(type, "biomass")) %>% 
  separate(type, c("type","yr")) %>% 
  mutate(yr = as.numeric(ifelse(is.na(yr),0,yr))) %>% 
  ggplot() +
  aes(x = yr, y = estimate) +
  geom_line() +
  geom_point(data = index, aes(x = Year-1935, y = Index/exp(q[1,2])),
             col = "darkgreen")



###########################################################
## Question 5

dll_name_pt <- "biomassdynamics"
if(is.loaded(dll_name_pt)) {
  dyn.unload(dynlib(dll_name_pt))
}
compile("biomassdynamics.cpp")
dyn.load(dynlib(dll_name_pt))

#om values
# estimated values from the original model, values used to generate the new data 
om<-sdreport(obj)%>%
  summary() %>%
  as.data.frame() %>%
  rownames_to_column(var = "parameter")%>%
  janitor::clean_names()%>%
  rename(om_val = estimate)%>%
  select(-std_error)

# generate siimulations
sim_data <- obj$simulate(complete=TRUE)# one sim first 
sim_data
sim_data <- map(1:500, function(x)obj$simulate(complete=TRUE))
str(sim_data, max.level=1)

#make a function to fit the model
refit<- function(data,parameters, dll_name){
  model_sim<- MakeADFun(data = data,
                        parameters = parameters,
                        DLL= dll_name, map = list(logSigmaC = factor(NA) ),
                        control=list(eval.max=1000,iter.max=1000,rel.tol=1e-5))
  fit<- nlminb(model_sim$par, model_sim$fn, model_sim$gr)
  sdr<- sdreport(model_sim)%>%
    summary() %>%
    as.data.frame()%>%
    rownames_to_column(var = "parameter")
  return(sdr)
}

sim_estimates <- map_dfr(sim_data, refit,
                         parameters = parameters,
                         dll_name = dll_name, .id = "isim") %>% janitor::clean_names()


par(mfcol=c(2,2))
logr<- sim_estimates %>% filter(parameter == "logr")
hist(logr$estimate, main="logr")

logK<- sim_estimates %>% filter(parameter == "logK")
hist(logK$estimate, main="logK")

B82<- sim_estimates %>% filter(parameter == "biomass.82")
hist(B82$estimate, main="Biomass end")

logq<- sim_estimates %>% filter(parameter == "logq")
hist(logq$estimate, main= "logq")


params <- sim_estimates  %>% filter(parameter %in% c("logr", "logK", "biomass.82","logq"))

plots <- params %>%
  left_join(om) %>%
  mutate(ree = (estimate-om_val)/om_val) %>%
  ggplot() + 
  aes(x = parameter, y = ree) +
  geom_boxplot() 


TMBAIC(opt)
#-147.987



  




