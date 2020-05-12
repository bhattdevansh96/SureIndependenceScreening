install.packages("SIS") #install the SIS package
install.packages("lars") #install the lars package
install.packages("flare") #install the flare package
install.packages("ncvreg") #install the ncvreg package
start = Sys.time() #record start for entire simulation
library(SIS) #import the SIS library
library(lars) #import the lars library
library(flare) #import the flare library
library(ncvreg) #import the ncvreg library

#creating simulation for SIS-SCAD with sample size 200 and 1000 variables. 
errors_SISSCAD_small = c() #create empty array to store errors from each run
model_sizes_SISSCAD_small = c() #create empty array to store selected model size from each simulation
start_time = Sys.time() #record start time for SIS-SCAD small
for(i in 1:200){ #testing on 200 different data sets
  
  set.seed(123*i) #re-create results at a later time but different seed for each data set
  n = 200 #sample size
  p = 1000 #variables
  x = matrix(rnorm(n*p, mean=0, sd=1), n, p) #creating IID standard Gaussian predictors

  set.seed(456*i) #re-create results at a later time but different seed for each data set
  s = 8 #true model size
  #compute non-zero beta parameters
  u = rbinom(s, 1, 0.4)
  z = rnorm(s, mean=0, sd=1)
  a = 4*log(n)/sqrt(n)
  b = ((-1)**(u))*(a + abs(z))
  #use non-zero beta components to create response vector
  y=x[, 1:s]%*%b
 
  #creating SIS-SCAD model and gaussian response. iter=FALSE means not doing ISIS
  modelSIS_small = SIS(x, y, family='gaussian', iter = FALSE, nsis=(n/log(n)))
  
  #create SCAD model using parameter space chosen by SIS
  modelSCAD = cv.ncvreg(x[,modelSIS_small$sis.ix0], y, alpha=1, family="gaussian", penalty="SCAD")
  modelNew = ncvreg(x[,modelSIS_small$sis.ix0], y, alpha=1, family="gaussian", penalty="SCAD", lambda=modelSCAD$lambda.min)
  
  predTest = predict(modelNew, X=x[,modelSIS_small$sis.ix0], type="response")
  a = (modelNew$beta[] != 0) 
  mse = mean((y - predTest)^2) #Using MSE to calculate the estimation errors
  rmse = sqrt(mse) #Square root of MSE (Square-loss function)
  errors_SISSCAD_small[i] = rmse #store the RMSE
 
  model_sizes_SISSCAD_small[i] = table(a)["TRUE"] #Store model size including intercept
  print(c(i, errors_SISSCAD_small[i], model_sizes_SISSCAD_small[i])) #print results

}
end_time = Sys.time() #record end-time for SIS-SCAD small
SISSCAD_total_time_small = end_time - start_time #record total time the SIS-SCAD small simulation took to run
SISSCAD_error_small_median = median(errors_SISSCAD_small) #store median of estimation error 
SISSCAD_model_Size_small_median = median(model_sizes_SISSCAD_small) #store median model size selected 

#creating simulation for SIS-SCAD with sample size 800 and 20000 variables
errors_SISSCAD_Large = c() #create empty array to store errors from each run
model_sizes_SISSCAD_Large = c() #create empty array to store selected model size from each run
start_time = Sys.time() 
for(i in 1:200){ #200 different data sets

  set.seed(123*i) #re-create results at a later time but different seed for each data set
  n = 800 #sample size
  p = 20000 #variables
  x = matrix(rnorm(n*p, mean=0, sd=1), n, p) #creating IID standard Normal predictors
  
  # gaussian response
  set.seed(456*i) #re-create results at a later time but different seed for each data set
  s = 18 #true model size
  #creating non-zero beta parameters
  u = rbinom(s, 1, 0.4)
  z = rnorm(s, mean=0, sd=1)
  a = 5*log(n)/sqrt(n)
  b= ((-1)**(u))*(a + abs(z))
  y=x[, 1:s]%*%b #using non-zero beta parameters to create response vector
  
  #creating SIS model and gaussian response. iter=FALSE means not doing ISIS
  modelSIS_Large = SIS(x, y, family='gaussian', iter = FALSE, nsis=(n/log(n)))
  
  modelSCAD = cv.ncvreg(x[,modelSIS_Large$sis.ix0], y, alpha=1, family="gaussian", penalty="SCAD")
  modelNew = ncvreg(x[,modelSIS_Large$sis.ix0], y, alpha=1, family="gaussian", penalty="SCAD", lambda=modelSCAD$lambda.min)
  
  predTest = predict(modelNew, X=x[,modelSIS_Large$sis.ix0], type="response")
  a = (modelNew$beta[] != 0) 
  mse = mean((y - predTest)^2) #using MSE for estimation errors
  rmse = sqrt(mse) #Square root of MSE (Square-loss function)
  errors_SISSCAD_Large[i] = rmse #store the RMSE
  
  model_sizes_SISSCAD_Large[i] = table(a)["TRUE"] #Store model size including intercept
  print(c(i, errors_SISSCAD_Large[i], model_sizes_SISSCAD_Large[i])) #print results
  
}
end_time = Sys.time() #record end time for SIS-SCAD large simulation
SISSCAD_total_time_Large = end_time - start_time
SISSCAD_error_Large_median = median(errors_SISSCAD_Large)
SISSCAD_model_Size_Large_median = median(model_sizes_SISSCAD_Large)

#creating simulation for LASSO wit sample size 200  and 1000 variables. 
errors_LASSO = c() #create empty array to store errors from each simulation
model_sizes_LASSO = c() #create empty array to store selected model size from each simulation
start_time = Sys.time()
for(i in 1:200){ #200 different data sets 

  set.seed(123*i) #re-create results at a later time but different seed for each data set
  n = 200 #sample size
  p = 1000 #variables
  x = matrix(rnorm(n*p, mean=0, sd=1), n, p) #creating IID standard Gaussian predictors
  
  # gaussian response
  set.seed(456*i) #re-create results at a later time but different seed for each data set
  s = 8
  u = rbinom(s, 1, 0.4)
  z = rnorm(s, mean=0, sd=1)
  a = 4*log(n)/sqrt(n)
  b= ((-1)**(u))*(a + abs(z))
  y=x[, 1:s]%*%b 
  
  #creating LASSO model
  modelLASSO = lars(x, y, type="lasso", normalize=FALSE, use.Gram = FALSE)
  
  predTest = predict.lars(modelLASSO, newx=x, type="fit") #create predictions using test data test
  mse = mean((y - predTest$fit[,9])^2) #compare predictions to real values of test y
  rmse = sqrt(mse) #Square root of MSE (Square-loss function)
  errors_LASSO[i] = rmse #store the RMSE
  
  model_sizes_LASSO[i] = length(modelLASSO$RSS) #Store model size including intercept
  print(c(i, errors_LASSO[i], model_sizes_LASSO[i])) #print results

}
end_time = Sys.time()
LASSO_total_time = end_time - start_time
LASSO_error_median = median(errors_LASSO)
LASSO_model_Size_median = median(model_sizes_LASSO)

#creating simulation for SIS-DS with sample size 800 and 20000 variables. 
errors_SISDS_Large = c() #create empty array to store errors from each simulation
model_sizes_SISDS_Large = c() #create empty array to store selected model size from each simulation
start_time = Sys.time()
for(i in 1:200){

  set.seed(123*i) #re-create results at a later time but different seed for each data set
  n = 800 #sample size
  p = 20000 #variables
  x = matrix(rnorm(n*p, mean=0, sd=1), n, p) #creating IID standard Gaussian random predictors
  
  # gaussian response
  set.seed(456*i) #re-create results at a later time but different seed for each data set
  s = 18
  u = rbinom(s, 1, 0.4)
  z = rnorm(s, mean=0, sd=1)
  a = 5*log(n)/sqrt(n)
  b= ((-1)**(u))*(a + abs(z))
  y=x[, 1:s]%*%b 
  
  #creating SIS-DS model and gaussian response. iter=FALSE means not doing ISIS
  modelSIS_Large = SIS(x, y, family='gaussian', iter = FALSE, nsis=(n/log(n)))
  
  modelDS = slim(x[,modelSIS_Large$sis.ix0], y, method="dantzig")
  
  predTest = predict(modelDS, newdata=x[,modelSIS_Large$sis.ix0]) #create predictions using test data test
  a = modelDS$beta[, 5] != 0
  mse = mean((y - predTest[[1]][,3])^2) #compare predictions to real values of test y
  rmse = sqrt(mse) #Square root of MSE (Square-loss function)
  errors_SISDS_Large[i] = rmse #store the RMSE
  
  model_sizes_SISDS_Large[i] = table(a)["TRUE"] #Store model size including intercept
  print(c(i, errors_SISDS_small[i], model_sizes_SISDS_Large[i])) #print results

}
end_time = Sys.time()
SISDS_total_time_Large = end_time - start_time
SISDS_error_Large_median = median(errors_SISDS_Large)
SISDS_model_Size_Large_median = median(model_sizes_SISDS_Large)

#creating simulation for SIS-DS with sample size 200 and 1000 variables. 
errors_SISDS_small = c() #create empty array to store errors from each simulation
model_sizes_SISDS_small = c() #create empty array to store selected model size from each simulation
start_time = Sys.time()
for(i in 1:200){
  
  set.seed(123*i) #re-create results at a later time but different seed for each data set
  n = 200 #sample size
  p = 1000 #variables
  x = matrix(rnorm(n*p, mean=0, sd=1), n, p) #creating IID standard Gaussian random predictors
  
  # gaussian response
  set.seed(456*i) #re-create results at a later time but different seed for each data set
  s = 8
  u = rbinom(s, 1, 0.4)
  z = rnorm(s, mean=0, sd=1)
  a = 4*log(n)/sqrt(n)
  b= ((-1)**(u))*(a + abs(z))
  y=x[, 1:s]%*%b 
  
  #creating SIS-DS model and gaussian response. iter=FALSE means not doing ISIS
  modelSIS_small = SIS(x, y, family='gaussian', iter = FALSE, nsis=(n/log(n)))
  
  modelDS = slim(x[,modelSIS_small$sis.ix0], y, method="dantzig")
  
  predTest = predict(modelDS, newdata=x[,modelSIS_small$sis.ix0]) #create predictions using test data test
  a = modelDS$beta[, 5] != 0
  mse = mean((y - predTest[[1]][,3])^2) #compare predictions to real values of test y
  rmse = sqrt(mse) #Square root of MSE (Square-loss function)
  errors_SISDS_small[i] = rmse #store the RMSE
  
  model_sizes_SISDS_small[i] = table(a)["TRUE"] #Store model size including intercept
  print(c(i, errors_SISDS_small[i], model_sizes_SISDS_small[i])) #print results

}
end_time = Sys.time()
SISDS_total_time_small = end_time - start_time
SISDS_error_small_median = median(errors_SISDS_small)
SISDS_model_Size_small_median = median(model_sizes_SISDS_small)

#creating simulation for SIS-DS with 200 sample size and 1000 variables. 
errors_DS_small = c() #create empty array to store errors from each simulation
model_sizes_DS_small = c() #creat empty array to store selected model size from each simulation
start_time = Sys.time()
for(i in 1:200){
  
  set.seed(123*i) #re-create results at a later time but different seed for each data set
  n = 200 #sample size
  p = 1000 #variables
  x = matrix(rnorm(n*p, mean=0, sd=1), n, p) #creating IID standard Gaussian random predictors
  
  # gaussian response
  set.seed(456*i) #re-create results at a later time but different seed for each data set
  s = 8
  u = rbinom(s, 1, 0.4)
  z = rnorm(s, mean=0, sd=1)
  a = 4*log(n)/sqrt(n)
  b= ((-1)**(u))*(a + abs(z))
  y=x[, 1:s]%*%b 
  
  modelDS = slim(x, y, method="dantzig")
  
  predTest = predict(modelDS, newdata=x) #create predictions using test data test
  a = modelDS$beta[, 5] != 0
  mse = mean((y - predTest[[1]][,3])^2) #compare predictions to real values of test y
  rmse = sqrt(mse) #Square root of MSE (Square-loss function)
  errors_DS_small[i] = rmse #store the RMSE
  
  model_sizes_DS_small[i] = table(a)["TRUE"] #Store model size including intercept
  print(c(i, errors_DS_small[i], model_sizes_DS_small[i])) #print results

}
end_time = Sys.time()
DS_total_time_small = end_time - start_time
DS_error_small_median = median(errors_DS_small)
DS_model_Size_small_median = median(model_sizes_DS_small)

print(c("SIS-SCAD Small", SISSCAD_model_Size_small_median, SISSCAD_error_small_median, SISSCAD_total_time_small))
print(c("SIS-SCAD LARGE", SISSCAD_model_Size_Large_median, SISSCAD_error_Large_median, SISSCAD_total_time_Large))
print(c("SIS-DS Small", SISDS_model_Size_small_median, SISDS_error_small_median, SISDS_total_time_small))
print(c("SIS-DS LARGE", SISDS_model_Size_Large_median, SISDS_error_Large_median, SISDS_total_time_Large))
print(c("DS", DS_model_Size_small_median, DS_error_small_median, DS_total_time_small))
print(c("LASSO", LASSO_model_Size_median, LASSO_error_median, LASSO_total_time))

end = Sys.time()
totaltime = end - start #total time to run whole simulation
print(totaltime)