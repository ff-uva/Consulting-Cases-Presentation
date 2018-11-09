#***************************************************************
#
#      Project 3: Simulation & Bootstrapping
#   		 A Case Study: Transplant Center
#   		   
#***************************************************************

#***************************************************************
#
#  Read in the data
#
#***************************************************************

#Set working directory

#Read data
r11xplant <- read.table("C:/Users/Think/Documents/Transplant/R11xplant.csv", sep = ",", header = T)

r11donor<-read.table("C:/Users/Think/Documents/Transplant/R11donor.csv", sep = ",", header = T)

uva <- read.table("C:/Users/Think/Documents/Transplant/UVAxplant.csv", sep = ",", header = T)

duke <- read.table("C:/Users/Think/Documents/Transplant/Dukexplant.csv", sep = ",", header = T)

mcv <- read.table("C:/Users/Think/Documents/Transplant/MCVxplant.csv", sep = ",", header = T)

unc <- read.table("C:/Users/Think/Documents/Transplant/UNCxplant.csv", sep = ",", header = T)



#    Source the bootstrapping functions
library(boot) #If you don't have this library, install it by: install.packages('boot')

source("C:/Users/Think/Downloads/TSbootfunctions.R")

source("C:/Users/Think/Downloads/SPM_Panel.R")

source("C:/Users/Think/Downloads/Transplant.plots.R")


#***************************************************************
#
# Part 1: Basic Statistics
#
#***************************************************************

#Step 1.1 Compare the performance of UVa with MCV kidney transplants



#Get the distribution of uva$Kidney, mcv$Kidney, r11donor$Kidney. What do you observe? 
kidney<-data.frame(uva$Kidney[-30],mcv$Kidney[-30],r11donor$Kidney[-30])
uva.pairs(as.matrix(kidney))
#We observed that the distribution for the uva center is much more symmetric than the other two
#The distribution for the mcv center is very asymmetric, with higher density on the extremes
#The mcv center has the highest correlation with region 11 kidney donors


#On average, how many kidney transplants are performed at UVa per year? MCV?
mean(uva$Kidney[-30])
mean(mcv$Kidney[-30])
#Excluding year 2017 data, 67.31 kidney transplants are performed at UVA  per year
#and 84.69 kidney transplants are performed at MCV per year, on average


#Perform a paired t-test between uva and mcv kidney transplants
#What is the result?
t.test(uva$Kidney[-30],mcv$Kidney[-30],paired=T)
#the paired t-test resulted in a p-value of about .006
#we reject the null hypothesis, indicating that there is a significant difference in the means of UVA and MCV kidney transplants


#Step 1.2 Compare the performance of UVa with Duke kidney transplants 

#Get the distribution of uva$Kidney, Duke$Kidney, r11xplant$Kidney. What do you observe?
kidney1<-data.frame(uva$Kidney[-30],duke$Kidney[-30],r11xplant$Kidney[-30])
uva.pairs(as.matrix(kidney1))
#we observed that the dist. for the UVA center is much more symmetric than the other two
#The dist. for the Duke center is roughly uniform, with a peak on the lower end of the range
#The Duke center has the higher correlation with region 11 kidney explants


#On average, how many kidney transplants are performed at UVa per year? Duke?
mean(uva$Kidney[-30])
mean(duke$Kidney[-30])
#Excluding year 2017 data, 67.31 kidney transplants are peformed at UVA per year
#and 87.79 kidney transplants are performed at Duke per year


#Perform a paired t-test between uva and duke kidney transplants
#What is the result?
t.test(uva$Kidney[-30],duke$Kidney[-30],paired=T)
#the paired t-test resulted in a p-value of about .000065
#we reject the null hypothesis, indicating that there is a significant difference in the means of UVA and Duke kidney transplants


#Step 1.3 Use bootstrapping to test the hypothesis: there is not a significant difference between UVa and MCV kidney transplants.

#What are the standard errors of the mean? What're the 95% confidence intervals? Do you accept or reject the null hypothesis?
bs.mean<-function(x,i)
{
  return(mean(x[i]))
}
uvamcv.diff<-ts(uva$Kidney-mcv$Kidney,1988,2016)
bs.uvamcv.diff<-boot(uvamcv.diff,bs.mean,R=10000)
bs.uvamcv.diff
boot.ci(bs.uvamcv.diff,0.95,type=c('bca','perc'))
#The standard errors of the mean is about 5.76
#The 95% confidence intervals are (-28.93, -6.41) [Percentile]; (-29.21, -6.72) [BCa]
#We reject the null hypothesis because the confidence intervals don't include zero


#Step 1.4 Use bootstrapping to test the hypothesis: There is not a significant difference between UVa and Duke kidney transplants.

#What are the standard errors of the mean? What're the 95% confidence intervals? Do you accept or reject the null hypothesis?
bs.mean<-function(x,i)
{
  return(mean(x[i]))
}
uvaduke.diff<-ts(uva$Kidney-duke$Kidney,1988,2016)
bs.uvaduke.diff<-boot(uvaduke.diff,bs.mean,R=10000)
bs.uvaduke.diff
boot.ci(bs.uvaduke.diff,0.95,type=c('bca','perc'))
#The standard error of the mean is about 4.29
#The 95% confidence intervals are (-29.07,-12.21) [Percentile]; (-29.38,-12.59) [BCa]
#We reject the null hypothesis because the confidence intervals don't include zero


#Step 1.5 Get the scatter plot matrix with the above 4 variables (UVA kidney, Duke kidney, R11 trnasplants, R11 donors). Describe what you observe. You can use either uva.pairs() {source("SPM_Panel.R")} or pairs().
uva.pairs(as.matrix(data.frame(uva$Kidney[-30],duke$Kidney[-30],r11xplant$Kidney[-30],r11donor$Kidney[-30])))
#we observe that the Duke center is more correlated with both Region 11 kidney transplants, as well as Region 11 kidney donors, than UVA
#The dist. for kidney transplants at UVA is somewhat symmetric with extreme modes in the center
#The dist. for kidney transplants at Duke is much less symmetric, closer to being uniform
#The UVA and Duke centers are almost equally correlated with Region 11 Donors
#The R11 Transplants and R11 donors for kidney transplants are extremely correlated (almost 1)

#***************************************************************
#
# Part 2: Linear Regression Models
#
#***************************************************************

# Test the hypothesis: There is no difference between the forecasted numbers of kidney 
# transplants that will performed at UVA and at MCV in 2017.

# Step 2.1 Build a linear model: 
# uva$Kidney-mcv$Kidney = b0+b1*r11donor$Kidney+e. Call it uva.kidney.lm.
# Analyze the result: R^2, model utility test, t-tests, etc.
uva.mcv.diff<-(uva$Kidney[-30]-mcv$Kidney[-30])
uva.kidney.lm<-lm(uva.mcv.diff~r11donor$Kidney[-30])
summary(uva.kidney.lm)
#We obtained an Adjusted R^2 value of .4122 and a Multiple R^2 value of .4332, indicating somewhat good fit
#The model  utility test resulted in a p-value of .0001042, indicating that the model is highly significant, even at the .01 level
#The t-test for the r11donor variable results in a p-value of .000104, indicating that the paramater is different from zero
#These results indicate that the r11 donor variable is very useful in predicting the difference between the mean of uva and mcv kidney transplants
#For every increase in r11donor$kidney, uva kidney transplants decreases relative to mcv by -.072


#Step 2.2 Generate the diagnostic plots. Interpret the results of the diagnostic plots. Do you see any problem?
par(mfrow=c(2,2))
plot(uva.kidney.lm)
par(mfrow=c(1,1))
#The Residuals vs Fitted plot shows a lack of fit and heteroscendasticity
#The QQ plot shows that the tails of the distribution are not Gaussian
#The Residuals vs Leverage plot shows there are multiple outliers, none of which are influential


# Step 2.3 Estimate the model with bootstrapping (by residuals). Is b1 significant?
# Get the fitted values from the regression model
uvamcv.fit<-fitted(uva.kidney.lm)
#    Get the residuals from the regression model
uvamcv.res<-residuals(uva.kidney.lm)
#    Get the regression model
uvamcv.mod<-model.matrix(uva.kidney.lm)
#   Bootstrapping LM
uvamcv.boot<-RTSB(uva.mcv.diff,r11donor$Kidney[-30],uvamcv.fit,uvamcv.res,uvamcv.mod,5000)
uvamcv.boot
sqrt(abs(var(uvamcv.boot$t)))

#    95% CI of r11donor
boot.ci(uvamcv.boot, .95, index=2)
#none of the confidence intervals include zero as a value

#    Distribution of b1
par(mfrow = c(1,2))
hist(uvamcv.boot$t[,2], main = "Region 11 Donors",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uvamcv.boot $t[,2])
qqline(uvamcv.boot $t[,2])
par(mfrow = c(1,1))

#    Is b1 significant?
#b1 is significant because all four of the 95% confidence intervals for b1 estimate do not include zero, implying that b1 is not equal to zero
#the histogram for region 11 donor coefficient values forms almost a normal distribution, centered around -.07, implying that b1 is significant


#Step 2.4* (bonus) What about Duke? Repeat the above steps and compare the results. 

# Test the hypothesis: There is no difference between the forecasted numbers of kidney 
# transplants that will be performed at UVA and at Duke in 2017.

# Build a linear model and analyze the results 
uva.duke.diff<-uva$Kidney[-30]-duke$Kidney[-30]
uva.duke.lm<-lm(uva.duke.diff~r11donor$Kidney[-30])
summary(uva.duke.lm)
#We obtained an Adjusted R^2 value of .04548 and a Multiple R^2 value of .079, indicating a worse fit than the previous model
#The model  utility test resulted in a p-value of .1382, indicating that the model is not significant, contrasted to the previous model comparing UVA and MCV
#The t-test for the r11donor variable results in a p-value of .138, indicating that the paramater is not different from zero
#These results indicate that the r11 donor variable is not useful in predicting the difference between the mean of uva and duke kidney transplants
#the r11donor variable is much worse at predicting the difference between UVA and Duke than it is at predicting UVA and MCV kidney transplants

#Generate the diagnostic plots. Interpret the results of the diagnostic plots. Do you see any problem?
par(mfrow=c(2,2))
plot(uva.duke.lm)
par(mfrow=c(1,1))
#The Residuals vs Fitted plot also shows a lack of fit and heteroscendasticity
#The QQ plot shows that the distrbution is Gausian, even at the tails, different from the previous model
#The Residuals vs Leverage plot again shows there are multiple outliers, none of which are influential


#Estimate the model with bootstrapping (by residuals). Is b1 significant?
# Get the fitted values from the regression model
uvaduke.fit<-fitted(uva.duke.lm)
#Get the residuals from the regression model
uvaduke.resid<-residuals(uva.duke.lm)
#Get the regression model
uvaduke.mod<-model.matrix(uva.duke.lm)
#Bootstrapping LM [UVA and Duke]
uvaduke.boot<-RTSB(uva.duke.diff,r11donor$Kidney[-30],uvaduke.fit,uvaduke.resid,uvaduke.mod,5000)
uvaduke.boot
sqrt(abs(var(uvaduke.boot$t)))
#95% CI of r11donor
boot.ci(uvaduke.boot, .95, index=2)
#All 4 confidence intervals include zero as a value, in contrast to the confidence interavals for the previous model

#    Distribution of b1
par(mfrow = c(1,2))
hist(uvaduke.boot$t[,2], main = "Region 11 Donors",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uvaduke.boot $t[,2])
qqline(uvaduke.boot $t[,2])
par(mfrow = c(1,1))

#    Is b1 significant?

#b1 is not significant for the model predicting the differences between UVA and Duke, in contrast to the previous model predicting the differences between UVA and MCV (which was significant)
#this is because all four of the confidence intervals for b1 included zero as a value
#in addition, the histogram for b1 has a much greater frequency for a coefficeint value of zero than for the previous model
#and the normal QQ plot has a sample quantile of zero that is much closer to a zero theoretical quantile than the previous model

#***************************************************************
#
# Part 3: Time Series Models
#
#***************************************************************
#Step 3.1 Generate the ACF and PACF plots of the residuals from your part 2 linear model for UVA and MCV kidney transplants. Interpret the results of the ACF and PACF plots. Do you recommend modeling the residuals? If so, what kind of model should you try based on these plots?  
acf(uvamcv.res)
pacf(uvamcv.res)
#The acf and pacf plots indicate that the the residuals are stationary - exponential decay in the ACF plot
#However, the plots show that there are correlated residuals, indicating that there is still structure to be modeled
#I would recommend modeling the residuals with an AR(1) model because the PACF plot cuts off after lag 1, but the ACF shows sinusoidal decay

#Step 3.2 Based on the above ACF and PACF plots, fit an ar model to the residuals

#    Add the AR model of the residuals to regression linear model based on the acf/pacf plots. 
#    Call this model uvamcv.kidney.lm2. Analyze the regression results.
uva.mcv.diff2<-(uva$Kidney[2:29]-mcv$Kidney[2:29])
uvamcv.kidney.lm2<- lm(uva.mcv.diff2~r11donor$Kidney[2:29]+uva.kidney.lm$residuals[1:28])
summary(uvamcv.kidney.lm2)
#The new model has an Adjusted R^2 of .787, indicating a great fit
#A model utility test resulted in a extremely significant p-value (very close to zero), indicating that the model is very significant
#Lastly, the p-values for both the r11donor variable and the residuals variable are significant, implying that the added residuals model is useful in predicting the difference in UVA and MCV kidney transplants

#Generate diagnostic plots for uvamcv.kidney.lm2. What are your observations?
par(mfrow=c(2,2))
plot(uvamcv.kidney.lm2)
par(mfrow=c(1,1))
#the residuals vs fitted plot indicates that while there is some heteroscendasticity, there is much less than in the previous model
#the normal QQ plot shows that the distribution is relatively normal
#there are also a few outliers, with one having a Leverage of close to .5
#overall, this model is a better fit than the previous model without the added residuals


#Step 3.3 Bootstrap the above time series model. Are the coefficients significant?

#    Get the fitted values from the regression model
uvamcv.lm2.fit<-fitted(uvamcv.kidney.lm2)
#    Get the residuals from the regression model
uvamcv.lm2.resid<-residuals(uvamcv.kidney.lm2)
#    Get the regression model
uvamcv.lm2.mod<-model.matrix(uvamcv.kidney.lm2)
#     Use the RTSB function to obtain the bootstrap
uvamcv.lm2.boot<-RTSB(uva.mcv.diff2,r11donor$Kidney[2:29],uvamcv.lm2.fit,uvamcv.lm2.resid,uvamcv.lm2.mod,5000)
#     The estimates
uvamcv.lm2.boot
summary(uvamcv.kidney.lm2)
uvamcv.lm2.boot$t
sqrt(abs(var(uvamcv.lm2.boot$t)))
boot.ci(uvamcv.lm2.boot, .95, index=2) #for r11donor variable
boot.ci(uvamcv.lm2.boot, .95, index=3) #for residual variable
#    Plot the results for the coeffiecient for region 11 donors
par(mfrow = c(1,2))
hist(uvamcv.lm2.boot$t[,2], main = "Region 11 Donors",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uvamcv.lm2.boot$t[,2])
par(mfrow = c(1,1))
#    Plot the results for the coeffiecient for time series components
par(mfrow = c(1,2))
hist(uvamcv.lm2.boot$t[,3], main = "Time Series Component",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uvamcv.lm2.boot$t[,3])
par(mfrow = c(1,1))

#    Are the coefficients significant?
#The coefficient for region 11 donors is very significant, due to the histogram of its coefficient values; there is a frequency of about zero for a coefficient value of 0
#The confidence intervals for the region 11 donors coefficient also do not contain zero in them
#The coefficient for the time series components is also very significant due to the histogram being centered around a coefficient value of .8
#The confidence intervals for the time sereies component coefficeint also do not contain zero in them.


#Step 3.5* (bonus) What about Duke? Repeat the above steps and compare the results. 

# Generate the ACF and PACF plots of the residuals from your part 2 linear model for UVA and Duke kidney transplants. Interpret the results of the ACF and PACF plots. Do you recommend modeling the residuals? If so, what kind of model should you try based on these plots?  
acf(uva.duke.lm$residuals)
pacf(uva.duke.lm$residuals)
#The plots show that there are correlated residuals, indicating that there is still structure to be modeled
#I would recommend modeling the residuals with an AR(2) model because the PACF plot cuts off after lag 2
#This is in contrast to UVA vs MCV in which we modeled the residuals with an AR(1) model
uva.duke.diff2<-(uva$Kidney[3:29]-duke$Kidney[3:29])
uvaduke.kidney.lm2<- lm(uva.duke.diff2~r11donor$Kidney[3:29]+uva.duke.lm$residuals[1:27]+uva.duke.lm$residuals[2:28])
summary(uvaduke.kidney.lm2)
#The new model has an Adjusted R^2 of .41, indicating a great fit
#A model utility test resulted in a extremely significant p-value (very close to zero), indicating that the model is very significant
#Lastly, the p-values for both the r11donor variable and the residuals variable are significant, implying that the added residuals model is useful in predicting the difference in UVA and MCV kidney transplants

#Generate diagnostic plots for uvamcv.kidney.lm2. What are your observations?
par(mfrow=c(2,2))
plot(uvaduke.kidney.lm2)
par(mfrow=c(1,1))
#the residuals vs fitted plot indicates that while there is some heteroscendasticity, there is much less than in the previous model
#the normal QQ plot shows that the distribution is relatively normal
#there are also a few outliers, with one having a Leverage of close to .5
#overall, this model is a better fit than the previous model without the added residuals
#    Get the fitted values from the regression model
uvaduke.lm2.fit<-fitted(uvaduke.kidney.lm2)
#    Get the residuals from the regression model
uvaduke.lm2.resid<-residuals(uvaduke.kidney.lm2)
#    Get the regression model
uvaduke.lm2.mod<-model.matrix(uvaduke.kidney.lm2)
#     Use the RTSB function to obtain the bootstrap
uvaduke.lm2.boot<-RTSB(uva.duke.diff2,r11donor$Kidney[3:29],uvaduke.lm2.fit,uvaduke.lm2.resid,uvaduke.lm2.mod,5000)
#     The estimates
uvaduke.lm2.boot
summary(uvaduke.kidney.lm2)
boot.ci(uvaduke.lm2.boot, .95, index=2) #for r11donor variable
#According to the confidence intervals, the r11donor variable is significant because Percentile and BCa don't bracket zero
#    Plot the results for the coeffiecient for region 11 donors
par(mfrow = c(1,2))
hist(uvaduke.lm2.boot$t[,2], main = "Region 11 Donors",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uvaduke.lm2.boot$t[,2])
par(mfrow = c(1,1))
#    Plot the results for the coeffiecient for time series components
par(mfrow = c(1,2))
hist(uvaduke.lm2.boot$t[,3], main = "Time Series Component",xlab ="Coefficient Values",   col = "steelblue", breaks = 50)
qqnorm(uvaduke.lm2.boot$t[,3])
par(mfrow = c(1,1))

#    Are the coefficients significant?
#The coefficient for region 11 donors is  significant, due to the histogram of its coefficient values; there is a frequency of about zero for a coefficient value of 0
#The confidence intervals for the region 11 donors coefficient also do not contain zero in them
#The coefficient for the time series components is also  significant due to the histogram being centered around a coefficient value of .8
#The confidence intervals for the time sereies component coefficeint also do not contain zero in them.
#***************************************************************
#
# Part 4: Predicting Differences in Kidney Transplants Part 1
#
#***************************************************************

#Step 4.1 Build an AR model to predict the difference in 2017
uvamcv.diff<-ts(uva$Kidney-mcv$Kidney,1988,2017)
uvamcv.ar<-ar(uvamcv.diff,method="yule-walker")

#Step 4.2 Use the predict function with ar model to forecast 2017 differences between UVA and MCV
uvamcv.pred<-predict(uvamcv.ar,newdata=uvamcv.diff)
uvamcv.pred<-predict(uvamcv.ar,se.fit=T,interval="predict")
uvamcv.pred
# Calculate the CI and plot the time series and prediction and CIs on a graph
# Plot the historical time series, the new prediction, and CIs on a graph
uvamcv.pred$pred+1.96*uvamcv.pred$se #upper CI
uvamcv.pred$pred-1.96*uvamcv.pred$se #lower CI
#The CI is (-45, 21)
plot(uvamcv.diff,ylim=c(-90,30),type='o') #historical time series
segments(2016,uvamcv.pred$pred,col = "red",2017) # prediction
segments(2016,uvamcv.pred$pred+1.96*uvamcv.pred$se,col = "red",2017,lty = "dashed") # upper CI
segments(2016,uvamcv.pred$pred-1.96*uvamcv.pred$se,col = "red",2017,lty = "dashed") # lower CI
# What do you observe?
#Observe that the predicting for 2017 is very close to the actual 2017 value, which is much greater than for 2016.
#The confidence interval is also very wide, with the lower bound being much around -45 and the upper bound being around 21

#Step 4.3 Bootstrapping the difference of UVa and MCV in 2017
#    To obtain a bootstrap estimate of the prediction for 2017
#    use the TSB function in the source file.

#    It takes three arguments:
#    tsint - the time series
#    oth.arg - the data for the new estimate
#    boot.number- number of replications (default=1000)

uvamcv.boot<-TSB(uvamcv.diff,uvamcv.pred$pred,5000)
uvamcv.boot
tsboot.ci(uvamcv.boot)
#    Interpret the results. Is there a significant difference between UVA and MCV in 2017?
#There is a significant difference between UVA and MCV in 2017
#the mean difference was estimated as -13, and the 95% CI was (-17.14,-9.3) which doesn't include zero
#The actual 2017 difference falls with in the bootstrapping interval

#Step 4.4* (bonus) What about Duke? Repeat the above steps and compare for Duke.
uvaduke.kid.diff<-ts(uva$Kidney-duke$Kidney,1988,2017)
uvaduke.ar.kid<-ar(uvaduke.kid.diff, method = "yule-walker")

uvaduke.kidney.pred<-predict(uvaduke.ar.kid,newdata=uvaduke.kid.diff)

plot(uvaduke.kid.diff,type='o')
segments(2016,uvaduke.kidney.pred$pred,col = "red",2017) # Prediction
segments(2016,uvaduke.kidney.pred$pred+1.96*uvaduke.kidney.pred$se,col = "red",2017,lty = "dashed") # upper CI
segments(2016,uvaduke.kidney.pred$pred-1.96*uvaduke.kidney.pred$se,col = "red",2017,lty = "dashed") # lower CI
#The prediction is very close to the actual 2017 difference. 
#The 95% CI is very wide and contains the actual 2017 difference 

diff.boot<-TSB(uvaduke.kid.diff,uvaduke.kidney.pred$pred,5000)
diff.boot
tsboot.ci(diff.boot)

#There is a significant difference between UVA and Duke in 2017
#the mean difference was estimated as -20.96, and the 95% CI was (-24.22,-17.65) which doesn't bracket zero


#***************************************************************
#
# Part 5: Predicting Differences in Kidney Transplants Part 2
#
#***************************************************************


#Step 5.1 Develop an AR model of region 11 kidney donors

#Plot an acf/pacf to estimate the model order for region 11 kidney donors - what model order do you suggest?
acf(r11donor$Kidney[-30])
pacf(r11donor$Kidney[-30])
#an AR(1) model is recommended because the pacf cuts off after lag 1

#Use ar() to fit an ar model to region 11 kidney donors from 1988-2016
r11donor.ar<-ar(r11donor$Kidney[-30],method="yule-walker")
r11donor.ar

#Step 5.2 Forecast the R11 donors and standard errors for 2017 using your ar model from step 5.1. Use forecast from the library(forecast).
library(forecast)
r11donor.pred <- forecast(r11donor.ar,h=1)
r11donor.pred

#Step 5.3 Use the linear model from part 3.2 combined with the forecast of region 11 kidney donors
#to forecast the differences in number of kidney transplants between UVa and MCV for 2017.
#Use the predict() function

#   Creating the new data frame
uvamcv.nd<-data.frame(r11k=r11donor.pred$mean,uvamcv.lm.e1=uva.kidney.lm$residuals[29])


#   Predict the linear model with the time series
uvamcv.new<-predict(uvamcv.kidney.lm2,newdata=uvamcv.nd,num=5000)

#Step 5.4 Bootstrap the Forecast from the linear model combined with the forecast of region 11 kidney donors to forecast the differences in number of kidney transplants between UVa and MCV for 2017.

#   Bootstrap prediction
r11donor.kidney<-r11donor$Kidney[-30]
r11k <- r11donor.kidney[2:29]
uvamcv.dm <- data.frame(uvamcv.diff[2:29], r11k, uva.kidney.lm$residuals[1:28])
r11donor.boot <- RFB(uvamcv.dm, model = uvamcv.kidney.lm2, ndata = uvamcv.nd, num=2000)

#   Bootstrap plot 
plot(r11donor.boot, index = 1)

#   Bootstrap confidence intervals
r11donor.boot.ci<-boot.ci(r11donor.boot, index=1)
r11donor.boot.ci
#   Interpret the results
#The bootstrapping the r11donor$Kidney is not significant because the confidence intervals include zero

#Step 5.5 Plot the current and predictions for each value along with the confidence intervals. Describe your observations.
uvamcv<-uva$Kidney-mcv$Kidney
ticks<-c(1988:2017)
plot(uvamcv~uva$Year,xlim=c(1988,2017), type="o",xlab="Years",ylab="No. of kidney")

segments(2017,r11donor.boot$t0[1],2017, r11donor.boot.ci$percent[4],lty="dashed")
segments(2017,r11donor.boot$t0[1],2017, r11donor.boot.ci$percent[5],lty="dashed")
#Observations?
#The actual difference for 2017 falls in the bootstrapping confidence interval



#Step 5.6* (bonus) What about Duke? Repeat the above steps for Duke.

#Plot an acf/pacf to estimate the model order for region 11 kidney donors - what model order do you suggest?
acf(r11donor$Kidney[-30])
pacf(r11donor$Kidney[-30])
#an AR(1) model is recommended because the pacf cuts off after lag 1

#Use ar() to fit an ar model to region 11 kidney donors from 1988-2016
r11donor.ar<-ar(r11donor$Kidney[-30],method="yule-walker")
r11donor.ar

#Step 5.2 Forecast the R11 donors and standard errors for 2017 using your ar model from step 5.1. Use forecast from the library(forecast).
library(forecast)
r11donor.pred <- forecast(r11donor.ar,h=1)
r11donor.pred

#Step 5.3 Use the linear model from part 3.2 combined with the forecast of region 11 kidney donors
#to forecast the differences in number of kidney transplants between UVa and MCV for 2017.
#Use the predict() function

#   Creating the new data frame
uvamcv.nd<-data.frame(r11k=r11donor.pred$mean,uvamcv.lm.e1=uva.kidney.lm$residuals[29])


#   Predict the linear model with the time series
uvamcv.new<-predict(uvamcv.kidney.lm2,newdata=uvamcv.nd,num=5000)

#Step 5.4 Bootstrap the Forecast from the linear model combined with the forecast of region 11 kidney donors to forecast the differences in number of kidney transplants between UVa and MCV for 2017.

uvaduke.diff<-ts(uva$Kidney-duke$Kidney,1988,2017)
uvaduke.dm <- data.frame(uvaduke.diff[2:29], r11k[1:28], uva.duke.lm$residuals[1:28])
r11donor.boot <- RFB(uvaduke.dm, model = uvaduke.kidney.lm2, ndata = uvaduke.nd, num=2000)

#   Bootstrap plot 
plot(r11donor.boot, index = 1)

#   Bootstrap confidence intervals
r11donor.boot.ci<-boot.ci(r11donor.boot, index=1)
r11donor.boot.ci
#   Interpret the results
#The bootstrapping the r11donor$Kidney is not significant because the confidence intervals include zero

#Step 5.5 Plot the current and predictions for each value along with the confidence intervals. Describe your observations.
uvamcv<-uva$Kidney-mcv$Kidney
ticks<-c(1988:2017)
plot(uvamcv~uva$Year,xlim=c(1988,2017), type="o",xlab="Years",ylab="No. of kidney")

segments(2017,r11donor.boot$t0[1],2017, r11donor.boot.ci$percent[4],lty="dashed")
segments(2017,r11donor.boot$t0[1],2017, r11donor.boot.ci$percent[5],lty="dashed")
#Observations?
#The actual difference for 2017 falls in the bootstrapping confidence interval