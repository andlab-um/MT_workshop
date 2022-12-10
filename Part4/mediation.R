library(mediation) #Mediation package
library(rockchalk) #Graphing simple slopes; moderation
library(multilevel) #Sobel Test
library(bda) #Another Sobel Test option
library(gvlma) #Testing Model Assumptions 
library(stargazer) #Handy regression tables


#lie:AUC
fitM <- lm(AUC ~ drate_m,data=MT) #IV on M; Hours since dawn predicting coffee consumption
fitY <- lm(rewant_z10 ~ AUC + drate_m,data=MT) #IV and M on DV; Hours since dawn and coffee predicting wakefulness
summary(fitM)
summary(fitY)
gvlma(fitM) #data is positively skewed; could log transform (see Chap. 10 on assumptions)
gvlma(fitY)
fitMed <- mediate(fitM, fitY, treat="drate_m", mediator="AUC")
summary(fitMed)
plot(fitMed)
#Bootstrap
fitMedBoot <- mediate(fitM, fitY, boot=T,  treat="drate_m", mediator="AUC")
summary(fitMedBoot)

