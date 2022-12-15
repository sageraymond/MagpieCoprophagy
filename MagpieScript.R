setwd("C:/Users/sager/OneDrive/Desktop/school/MSc/Survey Info/Recon_Maps/Scat/2022Analysis_and_MSPrep")

library(lme4)
library(survival)
library(Hmisc)
library(ggplot2)
library(AICcmodavg)
library(rafalib)
library(MuMIn)
library(MASS)
library(pROC)
library(car)
library(caret)
library(psych)
library(dplyr)
library(tidyr)
library(broom)
library(tidyverse)
library(magrittr)
library(irr)
library(splitstackshape)
library(jtools)
library(ggstance)
library(ggh4x)
library(grid)
library(gridExtra)
library(jtools)
library(interactions)
library(BAMMtools)
library(DescTools)
library(AMR)
library(patchwork)
library(reshape2)
library(AER)
library(glmmTMB)
library(performance)
library(rsq)
library(glmmTMB)
library(ggpubr)
library(corrplot)
library(magick)
library(sjPlot)

##Begn with contextual predictors##
#Read in dataframe containing contextual information
MagCont <- read.csv("AllScats_April52022.csv")

#remove unecessary columns
MagCont <- MagCont %>%
  dplyr::select(-c(ID, Timestamp, DateTime, Year, Day, LocationDeet,
                   Lat, Long, Content, TOSS, Comments.)) %>%
  dplyr::filter(Pee != "Y") %>%
  dplyr::filter(Magpie %in% c("Yes", "No"))

#Group into seasons based on desired divisions
MagCont$BiolSeason<- MagCont$Season
MagCont$Season <- ifelse(MagCont$Month == "January", "Winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "February", "Winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "March", "Winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "April", "Not winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "May", "Not winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "September", "Not winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "October", "Not winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "November", "Not winter", MagCont$Season)
MagCont$Season <- ifelse(MagCont$Month == "December", "Winter", MagCont$Season)

#Convert season to factor
MagCont$Season <- factor(MagCont$Season, levels = c("Not winter", "Winter"))

#Convert freshness to factor
MagCont$Freshness[is.na(MagCont$Freshness)] <- "Unk"
MagCont <- MagCont %>%
  dplyr::filter(Freshness != "Unk")
MagCont$Freshness <- ifelse(MagCont$Freshness == "Old", "Medium", MagCont$Freshness)
MagCont$Freshness <- ifelse(MagCont$Freshness == "Fresh", "Medium", MagCont$Freshness)
MagCont$Freshness <- factor(MagCont$Freshness, levels = c("Medium", "Vold", "Vfresh"))

#Convert size to factor
MagCont$Size <- factor(MagCont$Size, levels = c("Med", "Small", "Large", "Unk"))

#Remove unecessary columns
MagCont <- MagCont %>%
  dplyr::select(-(c(Month, Location, Unk, Pee, BiolSeason)))

#Refine data related to content
MagCont$Fruit <- ifelse(MagCont$Fruit == "Y", "Fruit", "No fruit")
MagCont$Anthro <- ifelse(MagCont$Anthro == "Y", "Anthropogenic", "No anthropogenic")
MagCont$Garbage <- ifelse(MagCont$Garbage == "Y", "Garbage", "No garbage")
MagCont$Veg <- ifelse(MagCont$Veg == "Y", "Vegetation", "No vegetation")
MagCont$Seed <- ifelse(MagCont$Seed == "Y", "Birdseed", "No birdseed")
MagCont$NatPrey <- ifelse(MagCont$NatPrey == "Y", "Natural prey", "No natural prey")

#Convert data related to content to factors
MagCont$Fruit <- factor(MagCont$Fruit, levels = c("No fruit", "Fruit"))
MagCont$Anthro <- factor(MagCont$Anthro, levels = c("No anthropogenic", "Anthropogenic"))
MagCont$Garbage <- factor(MagCont$Garbage, levels = c("No garbage", "Garbage"))
MagCont$Veg <- factor(MagCont$Veg, levels = c("No vegetation", "Vegetation"))
MagCont$Seed <- factor(MagCont$Seed, levels = c("No birdseed", "Birdseed"))
MagCont$NatPrey <- factor(MagCont$NatPrey, levels = c("No natural prey", "Natural prey"))

#Create outcome variable as factor
MagCont$Magpie2 <- ifelse(MagCont$Magpie == "Yes", "1", "0")
MagCont$Magpie2 <- factor(MagCont$Magpie2, levels = c(0,1))

#Create Univariate GLMs for each predictor
Mmod1 <- glm(Magpie2 ~ Season, family = binomial, data = MagCont)
Mmod2 <- glm(Magpie2 ~ Size, family = binomial, data = MagCont)
Mmod3 <- glm(Magpie2 ~ Freshness, family = binomial, data = MagCont)
Mmod4 <- glm(Magpie2 ~ Fruit, family = binomial, data = MagCont)
Mmod5 <- glm(Magpie2 ~ Anthro, family = binomial, data = MagCont)
Mmod6 <- glm(Magpie2 ~ Garbage, family = binomial, data = MagCont)
Mmod7 <- glm(Magpie2 ~ NatPrey, family = binomial, data = MagCont)
Mmod8 <- glm(Magpie2 ~ Seed, family = binomial, data = MagCont)
Mmod9 <- glm(Magpie2 ~ Veg, family = binomial, data = MagCont)

summ(Mmod1, digits = 4, confint = TRUE) #Winter increases coprophagy; retain
summ(Mmod2, digits = 4, confint = TRUE) #avoid size, select unk d; retain
summ(Mmod3, digits = 4, confint = TRUE) #Select older scats, avoid very fresh ones
summ(Mmod4, digits = 4, confint = TRUE) #Avoid fruit, but not sig
summ(Mmod5, digits = 4, confint = TRUE) #select anthro, but not sig
summ(Mmod6, digits = 4, confint = TRUE) #Avoid garbage, retain
summ(Mmod7, digits = 4, confint = TRUE) #avoid meat, retain
summ(Mmod8, digits = 4, confint = TRUE) #Avoid seed, not sig
summ(Mmod9, digits = 4, confint = TRUE) #Avoid veg, retain

#Summarise all info in a table
#export_summs(Mmod1, Mmod2, Mmod3, Mmod4, Mmod5, Mmod6, Mmod7, Mmod8, Mmod9,
             scale = TRUE,
             ci_level = 0.95,
             stars = NULL,
             error_format = '(p = {p.value}), CI [{conf.low}, {conf.high}]',
             model.names = c("Season", "Size", "Freshness", "Fruit", "Anthropogenic",
                             "Garbage", "Natural prey", "Birdseed", "Vegetation"),
             digits = 4,
             to.file = "xlsx",
             file.name = "ContextUnivGLMSOct27.xlsx")

#Get OR
summ(Mmod1, digits = 4, confint = TRUE, exp = TRUE) #Winter increases coprophagy; retain
summ(Mmod2, digits = 4, confint = TRUE, exp = TRUE) #avoid size, select unk d; retain
summ(Mmod3, digits = 4, confint = TRUE, exp = TRUE) #Select older scats, avoid very fresh ones
summ(Mmod4, digits = 4, confint = TRUE, exp = TRUE) #Avoid fruit, but not sig
summ(Mmod5, digits = 4, confint = TRUE, exp = TRUE) #select anthro, but not sig
summ(Mmod6, digits = 4, confint = TRUE, exp = TRUE) #Avoid garbage, retain
summ(Mmod7, digits = 4, confint = TRUE, exp = TRUE) #avoid meat, retain
summ(Mmod8, digits = 4, confint = TRUE, exp = TRUE) #Avoid seed, not sig
summ(Mmod9, digits = 4, confint = TRUE, exp = TRUE) #Avoid veg, retain

#Get mcfadden R2
PseudoR2(Mmod1, which = "McFadden") #0.06732318 
PseudoR2(Mmod2, which = "McFadden") #0.1499049 
PseudoR2(Mmod3, which = "McFadden") #0.06163845 
PseudoR2(Mmod4, which = "McFadden") #0.000450988 
PseudoR2(Mmod5, which = "McFadden") #0.002072263 
PseudoR2(Mmod6, which = "McFadden") #0.004850859 
PseudoR2(Mmod7, which = "McFadden") #0.009498226 
PseudoR2(Mmod8, which = "McFadden") #0.00027383 
PseudoR2(Mmod9, which = "McFadden") #0.008316569 

#Remove predictors associated with P > 0.25
MagCont <- MagCont %>%
  dplyr::select(-(c(Magpie, Fruit, Seed)))

#Dredge (all possible combinations)  
ContextualDredge <- glm(Magpie2 ~ ., family = binomial, data = MagCont)

options(na.action = "na.fail")
ContextRes <- dredge(ContextualDredge, trace = 2)

#Determine models within delta 2 AIC
subset(ContextRes, delta < 2)
Contextmod1a <- get.models(ContextRes, 1)[[1]]
Contextmod2a <- get.models(ContextRes, 2)[[1]]
Contextmod3a <- get.models(ContextRes, 3)[[1]]

#Get model formula
Contextmod1a$call
Contextmod2a$call
Contextmod3a$call

#Create models
Contextmod1 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                     Size + Veg + 1, family = binomial, data = MagCont)
Contextmod2 <- glm(formula = Magpie2 ~ Freshness + Garbage + Season + Size + 
                     Veg + 1, family = binomial, data = MagCont)
Contextmod3 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                     Size + 1, family = binomial, data = MagCont)

#Summarise models
summ(Contextmod1, confint = TRUE, digits = 4)
summ(Contextmod2, confint = TRUE, digits = 4)
summ(Contextmod3, confint = TRUE, digits = 4)

#Create Plot showing beta coefficients
PlotContext <- plot_summs(Contextmod1, Contextmod2, Contextmod3,
                        colors = c("black", "black", "black"),
                        point.size = 30,
                        point.shape = FALSE,
                        coefs = c("Winter" = "SeasonWinter",
                                  "Unknown diameter" = "SizeUnk",
                                  "Fresh (< 24 hours)" = "FreshnessVfresh",
                                  "Garbage (Content)" = "GarbageGarbage",
                                  "Vegetation (Content)" = "VegVegetation",
                                  "Anthropogenic (Content)" = "AnthroAnthropogenic",
                                  "Large diameter" = "SizeLarge",
                                  "Small diameter" = "SizeSmall",
                                  "Old (> 3 weeks)" = "FreshnessVold"))

PlotContext2 <- PlotContext + theme_classic()

PlotContext3 <- PlotContext2 +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Predictor estimate", title = "(A) Contextual predictors") +
  scale_shape_manual(values = c(24,24,24,24,24))

#Conduct model averaging
ContextAvg <- model.avg(Contextmod1, Contextmod2, Contextmod3)
summary(ContextAvg)
confint(ContextAvg, subset = TRUE)
exp(ContextAvg$coefficients)
exp(confint(ContextAvg, subset = TRUE))

#Try adding biologically plausible, a priori interaction terms to final model
Intmod1 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                 Size + Veg + Size*Season + 1, family = binomial, data = MagCont)

Intmod2 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                 Size + Veg + Anthro*Season + 1, family = binomial, data = MagCont)

Intmod3 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                 Size + Veg + Garbage*Season + 1, family = binomial, data = MagCont)

Intmod4 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                 Size + Veg + Veg*Season + 1, family = binomial, data = MagCont)

Intmod5 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                 Size + Veg + Freshness*Size + 1, family = binomial, data = MagCont)

Intmod6 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                 Size + Veg + Freshness*Anthro + 1, family = binomial, data = MagCont)

Intmod7 <- glm(formula = Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                 Size + Veg + Freshness*Garbage + 1, family = binomial, data = MagCont)

#Do the interaction terms improve model AIC?
AICc(Contextmod1) - AICc(Intmod1) #-3.067484
AICc(Contextmod1) - AICc(Intmod2) #-1.421161
AICc(Contextmod1) - AICc(Intmod3) #-1.757168
AICc(Contextmod1) - AICc(Intmod4) #-1.804529
AICc(Contextmod1) - AICc(Intmod5) #-9.766089
AICc(Contextmod1) - AICc(Intmod6) #-2.698476
AICc(Contextmod1) - AICc(Intmod7) #-1.091786
#No improvements

#Do CI for interaction terms overlap zero?
confint(Intmod1)
confint(Intmod2)
confint(Intmod3)
confint(Intmod4)
confint(Intmod5)
confint(Intmod6)
confint(Intmod7)
#Yes, they all do
#Conclusion = interactions do NOT improve model performance!

#Get summary metrics for top models
mean(c(PseudoR2(Contextmod1, which = "Nagelkerke"), PseudoR2(Contextmod2, 
                                                             which = "Nagelkerke"),
       PseudoR2(Contextmod3, which = "Nagelkerke"))) #0.3820695


PseudoR2(Contextmod1, which = "Nagelkerke") #0.3851283 
PseudoR2(Contextmod2, which = "Nagelkerke") #0.3818808 
PseudoR2(Contextmod3, which = "Nagelkerke") #0.3791996 

mean(c(PseudoR2(Contextmod1, which = "McFadden"), PseudoR2(Contextmod2, 
                                                           which = "McFadden"), 
       PseudoR2(Contextmod3, which = "McFadden"))) #0.2494719
PseudoR2(Contextmod1, which = "McFadden") #0.2518363 
PseudoR2(Contextmod2, which = "McFadden") #0.2493238 
PseudoR2(Contextmod3, which = "McFadden") #0.2472556 

#Do cross validation
data_ctrl <- trainControl(method = "cv", number = 5)
model_caret <- train(Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                       Size + Veg + 1, family = binomial, data = MagCont,
                     trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret) #Acc = 0.7636545 kappa = 0.4429998

model_caret2 <- train(Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                        Size + Veg + 1, family = binomial, data = MagCont,
                     trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret2) #Acc = 0.7606343   kappa = 0.4395773


model_caret3 <- train(Magpie2 ~ Anthro + Freshness + Garbage + Season + 
                        Size + 1, family = binomial, data = MagCont,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret3) #Acc = 0.7667307     kappa = 0.4508064

#mean accuracy
mean(c(0.7636545, 0.7606343, 0.7667307)) #0.7636732
#mean cohens kappa
mean(c(0.4429998, 0.4395773, 0.4508064)) #0.4444612
#Note that CV results will be different every time unless set.seed() is used

#ROC 
auc(roc(MagCont$Magpie2~fitted(Contextmod1))) #AUC = 0.8118
auc(roc(MagCont$Magpie2~fitted(Contextmod2))) #AUC = 0.8082
auc(roc(MagCont$Magpie2~fitted(Contextmod3))) #AUC = 0.8068

mean(c(auc(roc(MagCont$Magpie2~fitted(Contextmod1))), 
       auc(roc(MagCont$Magpie2~fitted(Contextmod2))), 
       auc(roc(MagCont$Magpie2~fitted(Contextmod3))))) #0.8089339
#End of contextual predictors#

###Move on to environmental predictors###
#Read in Data
Env <- read.csv("MasterScatOct28.csv")

#Remove all entries except those containing Magpie Info
Env <- Env %>%
  dplyr::filter(Magpie != "DNR") %>%
  dplyr::filter(Magpie!= "Unk")

#select relevant columns
Env <- Env %>%
  dplyr::select(-c(X, FID, DistEdge, EdDens, LC, Curve2, Curve3, Curve4, Curve5,
                   Curve6, Curve7, Aspect, AG, WATER, ScatCount, Ed_Len, DistSWater,
                   DistRiv, DistGarden, Use, AspRad, North))

#Create outcome variable column as factor
Env$Magpie2 <- ifelse(Env$Magpie == "Yes", "1", "0")
Env$Magpie2 <- factor(Env$Magpie2, levels = c(0,1))

#Improve microhabitat type location data
Env$Location <- ifelse(Env$Location == "Bed", "Bed or Den", Env$Location)
Env$Location <- ifelse(Env$Location == "camp", "Other", Env$Location)
Env$Location <- ifelse(Env$Location == "Camp", "Other", Env$Location)
Env$Location <- ifelse(Env$Location == "Den", "Bed or Den", Env$Location)
Env$Location <- ifelse(Env$Location == "Railroad", "Other", Env$Location)
Env$Location <- ifelse(Env$Location == "Yard", "Other", Env$Location)

#Code as factor
Env$Location <- factor(Env$Location, levels = c("Other", "Bed or Den", "Forage",
                                                "Junction", "OffTrail", "RiverCreek",
                                                "Trail"))

#Multiply densities by 1000 (so they are in XX/ km2, not m2)
Env$NatEdDens <- Env$NatEdDens*1000
Env$RoadDens <- Env$RoadDens*1000
Env$BldgDens <- Env$BldgDens*1000

#remove unecessary columns
Env <- Env %>%
  dplyr::select(-(c(Magpie, NatEdLen)))

#Ready to begin modelling approach
#Start by saving a copy of everything in its original form (which will facilitate calculating OR later)
Env$DistWaterOG <- Env$DistWater
Env$DistBldgOG <- Env$DistBldg
Env$DistRoadOG <- Env$DistRoad
Env$DistNatEdgeOG <- Env$DistNatEdge
Env$DistScatOG <- Env$DistScat
Env$NatEdDensOG <- Env$NatEdDens
Env$RoadDensOG <- Env$RoadDens
Env$BldgDensOG <- Env$BldgDens
Env$Curve1OG <- Env$Curve1
Env$SlopeOG <- Env$Slope
Env$ANTHOG <- Env$ANTH
Env$GRASSOG <- Env$GRASS
Env$NATOG <- Env$NAT
Env$DistJuncOG <- Env$DistJunc
Env$DistTrailOG <- Env$DistTrail
Env$DistCampOG <- Env$DistCamp
Env$DistPGOG <- Env$DistPG
Env$EastOG <- Env$East

#Scale and center all variables
Env$DistWater<- scale(Env$DistWater, scale = TRUE, center = TRUE)
Env$DistBldg<- scale(Env$DistBldg, scale = TRUE, center = TRUE)
Env$DistRoad<- scale(Env$DistRoad, scale = TRUE, center = TRUE)
Env$DistNatEdge<- scale(Env$DistNatEdge, scale = TRUE, center = TRUE)
Env$DistScat<- scale(Env$DistScat, scale = TRUE, center = TRUE)
Env$NatEdDens<- scale(Env$NatEdDens, scale = TRUE, center = TRUE)
Env$RoadDens<- scale(Env$RoadDens, scale = TRUE, center = TRUE)
Env$BldgDens<- scale(Env$BldgDens, scale = TRUE, center = TRUE)
Env$Curve1<- scale(Env$Curve1, scale = TRUE, center = TRUE)
Env$Slope<- scale(Env$Slope, scale = TRUE, center = TRUE)
Env$ANTH<- scale(Env$ANTH, scale = TRUE, center = TRUE)
Env$GRASS<- scale(Env$GRASS, scale = TRUE, center = TRUE)
Env$NAT<- scale(Env$NAT, scale = TRUE, center = TRUE)
Env$DistJunc<- scale(Env$DistJunc, scale = TRUE, center = TRUE)
Env$DistTrail<- scale(Env$DistTrail, scale = TRUE, center = TRUE)
Env$DistCamp<- scale(Env$DistCamp, scale = TRUE, center = TRUE)
Env$DistPG<- scale(Env$DistPG, scale = TRUE, center = TRUE)
Env$East<- scale(Env$East, scale = TRUE, center = TRUE)

Env$DistWater<- as.numeric(Env$DistWater)
Env$DistBldg<- as.numeric(Env$DistBldg)
Env$DistRoad<- as.numeric(Env$DistRoad)
Env$DistNatEdge<- as.numeric(Env$DistNatEdge)
Env$DistScat<- as.numeric(Env$DistScat)
Env$NatEdDens<- as.numeric(Env$NatEdDens)
Env$RoadDens<- as.numeric(Env$RoadDens)
Env$BldgDens<- as.numeric(Env$BldgDens)
Env$Curve1<- as.numeric(Env$Curve1)
Env$Slope<- as.numeric(Env$Slope)
Env$ANTH<- as.numeric(Env$ANTH)
Env$GRASS<- as.numeric(Env$GRASS)
Env$NAT<- as.numeric(Env$NAT)
Env$DistJunc<- as.numeric(Env$DistJunc)
Env$DistTrail<- as.numeric(Env$DistTrail)
Env$DistCamp<- as.numeric(Env$DistCamp)
Env$DistPG<- as.numeric(Env$DistPG)
Env$East<- as.numeric(Env$East)

#Create univariate GLMs
EnvMod1 <- glm(Magpie2 ~ DistWater, family = binomial, data = Env)
EnvMod2 <- glm(Magpie2 ~ DistBldg, family = binomial, data = Env)
EnvMod3 <- glm(Magpie2 ~ DistRoad, family = binomial, data = Env)
EnvMod4 <- glm(Magpie2 ~ DistNatEdge, family = binomial, data = Env)
EnvMod5 <- glm(Magpie2 ~ DistScat, family = binomial, data = Env)

EnvMod6 <- glm(Magpie2 ~ NatEdDens, family = binomial, data = Env)
EnvMod7 <- glm(Magpie2 ~ RoadDens, family = binomial, data = Env)
EnvMod8 <- glm(Magpie2 ~ BldgDens, family = binomial, data = Env)

EnvMod9 <- glm(Magpie2 ~ Curve1, family = binomial, data = Env)

EnvMod10 <- glm(Magpie2 ~ Slope, family = binomial, data = Env)

EnvMod11 <- glm(Magpie2 ~ ANTH, family = binomial, data = Env)
EnvMod12 <- glm(Magpie2 ~ GRASS, family = binomial, data = Env)
EnvMod13 <- glm(Magpie2 ~ NAT, family = binomial, data = Env)

EnvMod14 <- glm(Magpie2 ~ East, family = binomial, data = Env)

EnvMod15 <- glm(Magpie2 ~ DistJunc, family = binomial, data = Env)
EnvMod16 <- glm(Magpie2 ~ DistTrail, family = binomial, data = Env)

EnvMod17 <- glm(Magpie2 ~ DistCamp, family = binomial, data = Env)
EnvMod18 <- glm(Magpie2 ~ DistPG, family = binomial, data = Env)

EnvMod19 <- glm(Magpie2 ~ Location, family = binomial, data = Env)


plot_model(EnvMod1, type = "pred", terms = "DistWater")
plot_model(EnvMod2, type = "pred", terms = "DistBldg")
plot_model(EnvMod3, type = "pred", terms = "DistRoad")
plot_model(EnvMod4, type = "pred", terms = "DistNatEdge")
plot_model(EnvMod5, type = "pred", terms = "DistScat")

plot_model(EnvMod6, type = "pred", terms = "NatEdDens")
plot_model(EnvMod7, type = "pred", terms = "RoadDens")
plot_model(EnvMod8, type = "pred", terms = "BldgDens")

plot_model(EnvMod9, type = "pred", terms = "Curve1")

plot_model(EnvMod10, type = "pred", terms = "Slope")

plot_model(EnvMod11, type = "pred", terms = "ANTH")
plot_model(EnvMod12, type = "pred", terms = "GRASS")
plot_model(EnvMod13, type = "pred", terms = "NAT")

plot_model(EnvMod14, type = "pred", terms = "East")

plot_model(EnvMod15, type = "pred", terms = "DistJunc")
plot_model(EnvMod16, type = "pred", terms = "DistTrail")

plot_model(EnvMod17, type = "pred", terms = "DistCamp")
plot_model(EnvMod18, type = "pred", terms = "DistPG")

plot_model(EnvMod19, type = "pred", terms = "Location")


summary(EnvMod1) #MORE cop further from water
summary(EnvMod2) #more cop closer to buildings
summary(EnvMod3) #more cop closer to roads
summary(EnvMod4) #no influence of nat edge distance; do not retain
summary(EnvMod5) #No influence of dist scat; do not retain
summary(EnvMod6) #no influence of nat edge density; do not retain
summary(EnvMod7) #no influence of road density; do not retain
summary(EnvMod8) # marginal influence of building density
summary(EnvMod9) #no influence of curvature; do not retain
summary(EnvMod10) # more soprophagy on flatter slopes
summary(EnvMod11) # less cop in anthropogenic land cover
summary(EnvMod12) # less cop where there's more grass
summary(EnvMod13) # more cop where theres more nat
summary(EnvMod14) # more cop on west aspects
summary(EnvMod15) # more cop further from junctions
summary(EnvMod16) # more cop further from trails
summary(EnvMod17) # more cop further form camps
summary(EnvMod18) # more cop closer to playgrounds
summary(EnvMod19) # Not much influence of anything, but retain because some things are sig

#Develop Decay Terms
#Make decay terms, MANY, for each distance metric. Then choose the best one
Env <- Env %>%
  mutate(Dec_PG0.002 = 1-(exp(-0.002*DistPGOG))) %>%
  mutate(Dec_PG0.004 = 1-(exp(-0.004*DistPGOG))) %>%
  mutate(Dec_PG0.006 = 1-(exp(-0.006*DistPGOG))) %>%
  mutate(Dec_PG0.012 = 1-(exp(-0.012*DistPGOG))) %>%
  mutate(Dec_PG0.03 = 1-(exp(-0.03*DistPGOG))) %>%
  mutate(Dec_PG0.06 = 1-(exp(-0.06*DistPGOG))) %>%
  mutate(Dec_PG0.2 = 1-(exp(-0.2*DistPGOG))) %>%
  mutate(Dec_Water0.002 = 1-(exp(-0.002*DistWaterOG))) %>%
  mutate(Dec_Water0.004 = 1-(exp(-0.004*DistWaterOG))) %>%
  mutate(Dec_Water0.006 = 1-(exp(-0.006*DistWaterOG))) %>%
  mutate(Dec_Water0.012 = 1-(exp(-0.012*DistWaterOG))) %>%
  mutate(Dec_Water0.03 = 1-(exp(-0.03*DistWaterOG))) %>%
  mutate(Dec_Water0.06 = 1-(exp(-0.06*DistWaterOG))) %>%
  mutate(Dec_Water0.2 = 1-(exp(-0.2*DistWaterOG))) %>%
  mutate(Dec_Bldg0.002 = 1-(exp(-0.002*DistBldgOG))) %>%
  mutate(Dec_Bldg0.004 = 1-(exp(-0.004*DistBldgOG))) %>%
  mutate(Dec_Bldg0.006 = 1-(exp(-0.006*DistBldgOG))) %>%
  mutate(Dec_Bldg0.012 = 1-(exp(-0.012*DistBldgOG))) %>%
  mutate(Dec_Bldg0.03 = 1-(exp(-0.03*DistBldgOG))) %>%
  mutate(Dec_Bldg0.06 = 1-(exp(-0.06*DistBldgOG))) %>%
  mutate(Dec_Bldg0.2 = 1-(exp(-0.2*DistBldgOG))) %>%
  mutate(Dec_Road0.002 = 1-(exp(-0.002*DistRoadOG))) %>%
  mutate(Dec_Road0.004 = 1-(exp(-0.004*DistRoadOG))) %>%
  mutate(Dec_Road0.006 = 1-(exp(-0.006*DistRoadOG))) %>%
  mutate(Dec_Road0.012 = 1-(exp(-0.012*DistRoadOG))) %>%
  mutate(Dec_Road0.03 = 1-(exp(-0.03*DistRoadOG))) %>%
  mutate(Dec_Road0.06 = 1-(exp(-0.06*DistRoadOG))) %>%
  mutate(Dec_Road0.2 = 1-(exp(-0.2*DistRoadOG))) %>%
  mutate(Dec_NatEdge0.002 = 1-(exp(-0.002*DistNatEdgeOG))) %>%
  mutate(Dec_NatEdge0.004 = 1-(exp(-0.004*DistNatEdgeOG))) %>%
  mutate(Dec_NatEdge0.006 = 1-(exp(-0.006*DistNatEdgeOG))) %>%
  mutate(Dec_NatEdge0.012 = 1-(exp(-0.012*DistNatEdgeOG))) %>%
  mutate(Dec_NatEdge0.03 = 1-(exp(-0.03*DistNatEdgeOG))) %>%
  mutate(Dec_NatEdge0.06 = 1-(exp(-0.06*DistNatEdgeOG))) %>%
  mutate(Dec_NatEdge0.2 = 1-(exp(-0.2*DistNatEdgeOG))) %>%
  mutate(Dec_Scat0.002 = 1-(exp(-0.002*DistScatOG))) %>%
  mutate(Dec_Scat0.004 = 1-(exp(-0.004*DistScatOG))) %>%
  mutate(Dec_Scat0.006 = 1-(exp(-0.006*DistScatOG))) %>%
  mutate(Dec_Scat0.012 = 1-(exp(-0.012*DistScatOG))) %>%
  mutate(Dec_Scat0.03 = 1-(exp(-0.03*DistScatOG))) %>%
  mutate(Dec_Scat0.06 = 1-(exp(-0.06*DistScatOG))) %>%
  mutate(Dec_Scat0.2 = 1-(exp(-0.2*DistScatOG))) %>%
  mutate(Dec_Junc0.002 = 1-(exp(-0.002*DistJuncOG))) %>%
  mutate(Dec_Junc0.004 = 1-(exp(-0.004*DistJuncOG))) %>%
  mutate(Dec_Junc0.006 = 1-(exp(-0.006*DistJuncOG))) %>%
  mutate(Dec_Junc0.012 = 1-(exp(-0.012*DistJuncOG))) %>%
  mutate(Dec_Junc0.03 = 1-(exp(-0.03*DistJuncOG))) %>%
  mutate(Dec_Junc0.06 = 1-(exp(-0.06*DistJuncOG))) %>%
  mutate(Dec_Junc0.2 = 1-(exp(-0.2*DistJuncOG))) %>%
  mutate(Dec_Trail0.002 = 1-(exp(-0.002*DistTrailOG))) %>%
  mutate(Dec_Trail0.004 = 1-(exp(-0.004*DistTrailOG))) %>%
  mutate(Dec_Trail0.006 = 1-(exp(-0.006*DistTrailOG))) %>%
  mutate(Dec_Trail0.012 = 1-(exp(-0.012*DistTrailOG))) %>%
  mutate(Dec_Trail0.03 = 1-(exp(-0.03*DistTrailOG))) %>%
  mutate(Dec_Trail0.06 = 1-(exp(-0.06*DistTrailOG))) %>%
  mutate(Dec_Trail0.2 = 1-(exp(-0.2*DistTrailOG))) %>%
  mutate(Dec_Camp0.002 = 1-(exp(-0.002*DistCampOG))) %>%
  mutate(Dec_Camp0.004 = 1-(exp(-0.004*DistCampOG))) %>%
  mutate(Dec_Camp0.006 = 1-(exp(-0.006*DistCampOG))) %>%
  mutate(Dec_Camp0.012 = 1-(exp(-0.012*DistCampOG))) %>%
  mutate(Dec_Camp0.03 = 1-(exp(-0.03*DistCampOG))) %>%
  mutate(Dec_Camp0.06 = 1-(exp(-0.06*DistCampOG))) %>%
  mutate(Dec_Camp0.2 = 1-(exp(-0.2*DistCampOG)))

Decmods <- list()
for(i in c(1:9)) {
  if(i==1){decs <- c("DistWater", "Dec_Water0.002", "Dec_Water0.004", "Dec_Water0.006", "Dec_Water0.012", "Dec_Water0.03", "Dec_Water0.06", "Dec_Water0.2")}
  if(i==2){decs <- c("DistBldg", "Dec_Bldg0.002", "Dec_Bldg0.004", "Dec_Bldg0.006", "Dec_Bldg0.012", "Dec_Bldg0.03", "Dec_Bldg0.06", "Dec_Bldg0.2")}
  if(i==3){decs <- c("DistRoad", "Dec_Road0.002", "Dec_Road0.004", "Dec_Road0.006", "Dec_Road0.012", "Dec_Road0.03", "Dec_Road0.06", "Dec_Road0.2")}
  if(i==4){decs <- c("DistNatEdge", "Dec_NatEdge0.002", "Dec_NatEdge0.004", "Dec_NatEdge0.006", "Dec_NatEdge0.012", "Dec_NatEdge0.03", "Dec_NatEdge0.06", "Dec_NatEdge0.2")}
  if(i==5){decs <- c("DistScat", "Dec_Scat0.002", "Dec_Scat0.004", "Dec_Scat0.006", "Dec_Scat0.012", "Dec_Scat0.03", "Dec_Scat0.06", "Dec_Scat0.2")}
  if(i==6){decs <- c("DistJunc", "Dec_Junc0.002", "Dec_Junc0.004", "Dec_Junc0.006", "Dec_Junc0.012", "Dec_Junc0.03", "Dec_Junc0.06", "Dec_Junc0.2")}
  if(i==7){decs <- c("DistTrail", "Dec_Trail0.002", "Dec_Trail0.004", "Dec_Trail0.006", "Dec_Trail0.012", "Dec_Trail0.03", "Dec_Trail0.06", "Dec_Trail0.2")}
  if(i==8){decs <- c("DistCamp", "Dec_Camp0.002", "Dec_Camp0.004", "Dec_Camp0.006", "Dec_Camp0.012", "Dec_Camp0.03", "Dec_Camp0.06", "Dec_Camp0.2")}
  if(i==9){decs <- c("DistPG", "Dec_PG0.002", "Dec_PG0.004", "Dec_PG0.006", "Dec_PG0.012", "Dec_PG0.03", "Dec_PG0.06", "Dec_PG0.2")}
  
  # Create a data frame to store the results from all the models.
  Decmods[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(decs)){
    
    model <- glm(Magpie2 ~ Env[,j], # Create the predictive model
                 data=Env, family = binomial)
    
    null <- glm(Magpie2 ~1, # Create the null model
                data=Env, family = binomial)
    
    
    Decmods[[i]][j,1] <- j # First column is the predictor
    Decmods[[i]][j,2] <- AIC(model) # Second column is AIC
  }
  
  colnames(Decmods[[i]]) <- c("predictor", "AIC")
  # rm(model, decs) # Keep workspace clean
}
DecModels <- dplyr::bind_rows(Decmods)
DecModels
#write.csv(DecModels, file = "DecayMagpieModelOct29.csv")

#Water is best with 0.002, AIC = 836.8111
#Bldg is best with no decay term, AIC = 872.9328
#Road is best with no decay term, AIC = 873.5134
#Dist nat edge is best with 0.2; AIC = 880.6045
#DistScat best with 0.03; AIC = 879.0253
#DistJunc is best with 0.002; AIC = 874.5161
#DistTrail is best with 0.004; AIC = 865.6384
#DistCamp is best with no decay; AIC = 861.4486
#DistPG is best with 0.002; AIC = 867.6228

#Create decay terms
Env <- Env %>%
  dplyr::select(DistWater, DistBldg, DistRoad, DistNatEdge, DistScat, NatEdDens,
                RoadDens, BldgDens, Curve1, Slope, ANTH, GRASS, NAT, DistJunc,
                DistTrail, DistCamp, DistPG, Location, East, Magpie2,
                DistWaterOG, DistBldgOG, DistRoadOG, DistNatEdgeOG, DistScatOG, 
                NatEdDensOG, RoadDensOG, BldgDensOG, Curve1OG, SlopeOG, ANTHOG,
                GRASSOG, NATOG, DistJuncOG, DistTrailOG, DistCampOG, DistPGOG, 
                EastOG, Dec_PG0.002, Dec_Water0.002, Dec_NatEdge0.2, Dec_Scat0.03,
                Dec_Junc0.002, Dec_Trail0.004) %>%
  dplyr::rename("DecPG" = Dec_PG0.002,
                "DecWater" = Dec_Water0.002,
                "DecNatEdge" = Dec_NatEdge0.2,
                "DecScat" = Dec_Scat0.03,
                "DecJunc" = Dec_Junc0.002,
                "DecTrail" = Dec_Trail0.004)

#Create original versions of decay metrics (to facilitate calculating OR later)
Env$DecPGOG <- Env$DecPG
Env$DecJuncOG <- Env$DecJunc
Env$DecNatEdgeOG <- Env$DecNatEdge
Env$DecScatOG <- Env$DecScat
Env$DecTrailOG <- Env$DecTrail
Env$DecWaterOG <- Env$DecWater

#Scale and center the decay terms
Env$DecJunc <- scale(Env$DecJunc, scale = TRUE, center = TRUE)
Env$DecNatEdge <- scale(Env$DecNatEdge, scale = TRUE, center = TRUE)
Env$DecPG <- scale(Env$DecPG, scale = TRUE, center = TRUE)
Env$DecScat <- scale(Env$DecScat, scale = TRUE, center = TRUE)
Env$DecTrail <- scale(Env$DecTrail, scale = TRUE, center = TRUE)
Env$DecWater <- scale(Env$DecWater, scale = TRUE, center = TRUE)

Env$DecJunc <- as.numeric(Env$DecJunc)
Env$DecNatEdge <- as.numeric(Env$DecNatEdge)
Env$DecPG <- as.numeric(Env$DecPG)
Env$DecScat <- as.numeric(Env$DecScat)
Env$DecTrail <- as.numeric(Env$DecTrail)
Env$DecWater <- as.numeric(Env$DecWater)

#Create a table containing summary information for all univariate GLMs
#Table should include info for linear and quadratic
#Table will include decay results

#Create a table containing information about linear univariate GLMs
Envmods_lin <- list()
for(i in c(1:25)) {
  if(i==1){Environs <- ("DistWater")}
  if(i==2){Environs <- ("DistBldg")}
  if(i==3){Environs <- ("DistRoad")}
  if(i==4){Environs <- ("DistNatEdge")}
  if(i==5){Environs <- ("DistScat")}
  if(i==6){Environs <- ("NatEdDens")}
  if(i==7){Environs <- ("RoadDens")}
  if(i==8){Environs <- ("BldgDens")}
  if(i==9){Environs <- ("Curve1")}
  if(i==10){Environs <- ("Slope")}
  if(i==11){Environs <- ("ANTH")}
  if(i==12){Environs <- ("GRASS")}
  if(i==13){Environs <- ("NAT")}
  if(i==14){Environs <- ("DistJunc")}
  if(i==15){Environs <- ("DistTrail")}
  if(i==16){Environs <- ("DistCamp")}
  if(i==17){Environs <- ("GRASS")}
  if(i==18){Environs <- ("DistPG")}
  if(i==19){Environs <- ("East")}
  if(i==20){Environs <- ("DecWater")}
  if(i==21){Environs <- ("DecNatEdge")}
  if(i==22){Environs <- ("DecTrail")}
  if(i==23){Environs <- ("DecScat")}
  if(i==24){Environs <- ("DecJunc")}
  if(i==25){Environs <- ("DecPG")}
  
  # Create a data frame to store the results from all the models.
  Envmods_lin[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(Environs)){
    
    model <- glm(Magpie2 ~ Env[,j], # Create the predictive model
                 data=Env, family = binomial)
    
    Envmods_lin[[i]][j,1] <- j # First column is the predictor
    Envmods_lin[[i]][j,2] <- AIC(model) # Second column is AIC
    Envmods_lin[[i]][j,3:4] <- model$coefficients # thirs and fourth column are beta for intercept and predictor
    Envmods_lin[[i]][j,5:6] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Envmods_lin[[i]][j,7:10] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Envmods_lin[[i]]) <- c("predictor", "AIC", "Intercept_Beta", "Predictor_beta", "Intercept_pval", "Predictor_pval", "Intercept_LowerCI", "Predictor_LowerCI", "Intercept_UpperCI", "Predictor_UpperCI")
  # rm(model, decs) # Keep workspace clean
}
EnvUnivariates_lin <- dplyr::bind_rows(Envmods_lin)

EnvUnivariates_lin <- EnvUnivariates_lin %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictor_pval, Predictor_LowerCI, Predictor_UpperCI) %>%
  dplyr::rename("Beta" = Predictor_beta,
                "pvalue" = Predictor_pval,
                "LowCI" = Predictor_LowerCI,
                "UppCI" = Predictor_UpperCI)

#Repeat this process for quadratic terms
Envmods_quad <- list()
for(i in c(1:25)) {
  if(i==1){Environs <- ("DistWater")}
  if(i==2){Environs <- ("DistBldg")}
  if(i==3){Environs <- ("DistRoad")}
  if(i==4){Environs <- ("DistNatEdge")}
  if(i==5){Environs <- ("DistScat")}
  if(i==6){Environs <- ("NatEdDens")}
  if(i==7){Environs <- ("RoadDens")}
  if(i==8){Environs <- ("BldgDens")}
  if(i==9){Environs <- ("Curve1")}
  if(i==10){Environs <- ("Slope")}
  if(i==11){Environs <- ("ANTH")}
  if(i==12){Environs <- ("GRASS")}
  if(i==13){Environs <- ("NAT")}
  if(i==14){Environs <- ("DistJunc")}
  if(i==15){Environs <- ("DistTrail")}
  if(i==16){Environs <- ("DistCamp")}
  if(i==17){Environs <- ("GRASS")}
  if(i==18){Environs <- ("DistPG")}
  if(i==19){Environs <- ("East")}
  if(i==20){Environs <- ("DecWater")}
  if(i==21){Environs <- ("DecNatEdge")}
  if(i==22){Environs <- ("DecTrail")}
  if(i==23){Environs <- ("DecScat")}
  if(i==24){Environs <- ("DecJunc")}
  if(i==25){Environs <- ("DecPG")}
  
  # Create a data frame to store the results from all the models.
  Envmods_quad[[i]] <- data.frame()
  
  # For every variable you want to test:
  for(j in c(Environs)){
    
    model <- glm(Magpie2 ~ Env[,j] + I(Env[,j]^2), # Create the predictive model
                 data=Env, family = binomial)
    
    Envmods_quad[[i]][j,1] <- j # First column is the predictor
    Envmods_quad[[i]][j,2] <- AIC(model) # Second column is AIC
    Envmods_quad[[i]][j,3:5] <- model$coefficients # thirs, fourth,fifth column are beta for intercept and predictor, predictor^2
    Envmods_quad[[i]][j,6:8] <- coef(summary(model))[,4] # fifth and sixth column are p values for intercept and predictor
    Envmods_quad[[i]][j,9:14] <- confint(model) # fifth and sixth column are p values for intercept and predictor
    
    
  }
  
  colnames(Envmods_quad[[i]]) <- c("predictor", "AIC", "Intercept_Beta", 
                                     "Predictor_beta", "Predictore_beta2", 
                                     "Intercept_pval", "Predictor_pval", 
                                     "Predictor_pval2", "Intercept_LowerCI", 
                                     "Predictor_LowerCI", "Predictor_LowerCI2",
                                     "Intercept_UpperCI", "Predictor_UpperCI",
                                     "Predictor_UpperCI2")
  # rm(model, decs) # Keep workspace clean
}
EnvUnivariates_quad <- dplyr::bind_rows(Envmods_quad)

EnvUnivariates_quad <- EnvUnivariates_quad %>%
  dplyr::select(predictor, AIC, Predictor_beta, Predictore_beta2,
                Predictor_pval, Predictor_pval2, Predictor_LowerCI, Predictor_UpperCI,
                Predictor_LowerCI2, Predictor_UpperCI2) %>%
  dplyr::rename("Q_Beta(lin)" = Predictor_beta,
                "Q_Beta(quad)" = Predictore_beta2,
                "Q_pvalue(lin)" = Predictor_pval,
                "Q_pvalue(quad)" = Predictor_pval2,
                "Q_LowCI(lin)" = Predictor_LowerCI,
                "Q_UppCI(lin)" = Predictor_UpperCI,
                "Q_LowCI(quad)" = Predictor_LowerCI2,
                "Q_UppCI(quad)" = Predictor_UpperCI2,
                "AIC_quad" = AIC)

EnvUnivariates <- cbind(EnvUnivariates_lin, EnvUnivariates_quad)
#Write to csv
#write.csv(EnvUnivariates, file = "MagpieEnvironUnivariates_Oct29.csv")

#Remove predictors that are better as decay terms and predictors that are associated with P > 0.250
Env2 <- Env %>%
  dplyr::select(DistBldg, DistRoad, BldgDens, Curve1, Slope, ANTH, GRASS, NAT,
                DistCamp, East, DecWater, DecNatEdge, DecTrail, DecScat, DecJunc,
                DecPG)

#Determine correlated pairs of predictors
cor(Env2[sapply(Env2, is.numeric)], method = c("pearson"))
#Slope and dec trail are correlated --> Trail performs better
#Slope and dec junc are correlated --> Slope performs better
#ANTH and Nat are correlated --> NAT performs better
#GRASS and NAT Are correlated --> NAT performs better
#Dec Junc and De Trail are correlated --> Dec trail performs better 

#REMOVE GRASS, ANTH, Dec Junc, Slope
Env3 <- Env %>%
  dplyr::select(DistBldg, DistRoad, BldgDens, Curve1, NAT,
                DistCamp, East, DecWater, DecNatEdge, DecTrail, DecScat,
                DecPG, Location, Magpie2)

#Complete all possible combinations
MagpieDredge <- glm(Magpie2 ~ DistBldg + DistRoad + I(DistRoad^2) + BldgDens + 
                      Curve1 + I(Curve1^2) + NAT + DistCamp + East +
                      DecWater + DecNatEdge + DecTrail + DecScat + DecPG + 
                      I(DistCamp^2),
                    family = binomial, data = Env3)

#specify that, when quadratics are included, linear term is also included
subsetMP <- expression(dc(`DistRoad`, `I(DistRoad^2)`) &
                         dc(`DistCamp`, `I(DistCamp^2)`) &
                         dc(`Curve1`, `I(Curve1^2)`))

options(na.action = "na.fail")
MagpieRes <- dredge(MagpieDredge, subset = subsetMP, trace = 2)

#Get models within 2 AIC
subset(MagpieRes, delta < 2)

#view models
EnvMagmod1a <- get.models(MagpieRes, 1)[[1]]
EnvMagmod2a <- get.models(MagpieRes, 2)[[1]]
EnvMagmod3a <- get.models(MagpieRes, 3)[[1]]
EnvMagmod4a <- get.models(MagpieRes, 4)[[1]]
EnvMagmod5a <- get.models(MagpieRes, 5)[[1]]
EnvMagmod6a <- get.models(MagpieRes, 6)[[1]]
EnvMagmod7a <- get.models(MagpieRes, 7)[[1]]
EnvMagmod8a <- get.models(MagpieRes, 8)[[1]]
EnvMagmod9a <- get.models(MagpieRes, 9)[[1]]
EnvMagmod10a <- get.models(MagpieRes, 10)[[1]]
EnvMagmod11a <- get.models(MagpieRes, 11)[[1]]

#Create models
EnvMagmod1a$call
EnvMagmod1 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 1, family = binomial, 
                  data = Env3)

EnvMagmod2a$call
EnvMagmod2 <- glm(formula = Magpie2 ~ BldgDens + DecPG + DecScat + DecTrail + 
                    DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                    1, family = binomial, data = Env3)

EnvMagmod3a$call
EnvMagmod3 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + NAT + 1, family = binomial, 
                  data = Env3)

EnvMagmod4a$call
EnvMagmod4 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + East + NAT + 1, family = binomial, 
                  data = Env3)

EnvMagmod5a$call
EnvMagmod5 <- glm(formula = Magpie2 ~ BldgDens + DecPG + DecScat + DecTrail + 
                    DecWater + DistCamp + I(DistCamp^2) + DistRoad + NAT + 1, 
                  family = binomial, data = Env3)

EnvMagmod6a$call
EnvMagmod6 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + NAT + 1, family = binomial, data = Env3)

EnvMagmod7a$call
EnvMagmod7 <- glm(formula = Magpie2 ~ DecNatEdge + DecPG + DecScat + DecTrail + 
                    DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                    1, family = binomial, data = Env3)

EnvMagmod8a$call
EnvMagmod8 <- glm(formula = Magpie2 ~ BldgDens + DecPG + DecScat + DecTrail + 
                    DecWater + DistCamp + I(DistCamp^2) + East + NAT + 1, family = binomial, 
                  data = Env3)

EnvMagmod9a$call
EnvMagmod9 <- glm(formula = Magpie2 ~ Curve1 + DecPG + DecScat + DecTrail + 
                    DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                    1, family = binomial, data = Env3)

EnvMagmod10a$call
EnvMagmod10 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                     DistCamp + I(DistCamp^2) + DistRoad + I(DistRoad^2) + East + 
                     NAT + 1, family = binomial, data = Env3)

EnvMagmod11a$call
EnvMagmod11 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                     DistBldg + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                     1, family = binomial, data = Env3)


#Add biologically plausible interaction terms
EnvIntmod1 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + Curve1 + 
                    DecWater:Curve1 + 1, family = binomial, 
                  data = Env3)


EnvIntmod2 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                    DistCamp:DecTrail + 1, family = binomial, 
                  data = Env3)

EnvIntmod3 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                    DistCamp:NAT + 1, family = binomial, 
                  data = Env3)

EnvIntmod4 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                    DistCamp:DistRoad + 1, family = binomial, 
                  data = Env3)

EnvIntmod5 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + BldgDens +
                    DistCamp:BldgDens + 1, family = binomial, 
                  data = Env3)

EnvIntmod6 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    DecScat:NAT + 1, family = binomial, 
                  data = Env3)

EnvIntmod7 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    DecScat:NAT + 1, family = binomial, 
                  data = Env3)

EnvIntmod8 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    DecScat:DistRoad + 1, family = binomial, 
                  data = Env3)

EnvIntmod9 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + BldgDens +
                    DecScat:BldgDens + 1, family = binomial, 
                  data = Env3)

EnvIntmod10 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    DecTrail:NAT + 1, family = binomial, 
                  data = Env3)

EnvIntmod11 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    DecTrail:DistRoad + 1, family = binomial, 
                  data = Env3)

EnvIntmod12 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    DecPG:NAT + 1, family = binomial, 
                  data = Env3)

EnvIntmod13 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    DecPG:DistRoad + 1, family = binomial, 
                  data = Env3)

EnvIntmod14 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT +
                    NAT:DistRoad + 1, family = binomial, 
                  data = Env3)

EnvIntmod15 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + BldgDens +
                    NAT:BldgDens + 1, family = binomial, 
                  data = Env3)

EnvIntmod16 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                    DistCamp + I(DistCamp^2) + DistRoad + East + NAT + DecNatEdge +
                    NAT:DecNatEdge + 1, family = binomial, 
                  data = Env3)

EnvIntmod17 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                     DistCamp + I(DistCamp^2) + DistRoad + East + NAT + DistBldg +
                     DistBldg:DistRoad + 1, family = binomial, 
                   data = Env3)

EnvIntmod18 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                     DistCamp + I(DistCamp^2) + DistRoad + East + NAT + BldgDens +
                     BldgDens:DistRoad + 1, family = binomial, 
                   data = Env3)

EnvIntmod19 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                     DistCamp + I(DistCamp^2) + DistRoad + East + NAT + Curve1 +
                     East:Curve1 + 1, family = binomial, 
                   data = Env3)

EnvIntmod20 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                     DistCamp + I(DistCamp^2) + DistRoad + East + NAT + DistBldg + BldgDens +
                     BldgDens:DistBldg + 1, family = binomial, 
                   data = Env3)

#Do interaction terms improve model AIC?
AICc(EnvMagmod1) - AICc(EnvIntmod1) #-3.776037
AICc(EnvMagmod1) - AICc(EnvIntmod2) #-1.655664
AICc(EnvMagmod1) - AICc(EnvIntmod3) #-1.884955
AICc(EnvMagmod1) - AICc(EnvIntmod4) #-1.957818
AICc(EnvMagmod1) - AICc(EnvIntmod5) #1.429384; improvement, but not enough for inclusion
AICc(EnvMagmod1) - AICc(EnvIntmod6) #-1.827956
AICc(EnvMagmod1) - AICc(EnvIntmod7) #-1.827956
AICc(EnvMagmod1) - AICc(EnvIntmod8) #-1.668372
AICc(EnvMagmod1) - AICc(EnvIntmod9) #-2.437316
AICc(EnvMagmod1) - AICc(EnvIntmod10) #2.467051; improvement sufficient for inclusion
AICc(EnvMagmod1) - AICc(EnvIntmod11) #0.8253529; improvement, but not enough for inclusion
AICc(EnvMagmod1) - AICc(EnvIntmod12) #-0.0335327
AICc(EnvMagmod1) - AICc(EnvIntmod13) #-0.4638404
AICc(EnvMagmod1) - AICc(EnvIntmod14) #-2.051653
AICc(EnvMagmod1) - AICc(EnvIntmod15) #-1.886294
AICc(EnvMagmod1) - AICc(EnvIntmod16) #-3.607521
AICc(EnvMagmod1) - AICc(EnvIntmod17) #-3.861703
AICc(EnvMagmod1) - AICc(EnvIntmod18) #-2.546245
AICc(EnvMagmod1) - AICc(EnvIntmod19) #-3.824595
AICc(EnvMagmod1) - AICc(EnvIntmod20) #-3.568773

#Repeat dredge
MagpieDredge2 <- glm(Magpie2 ~ DistBldg + DistRoad + I(DistRoad^2) + BldgDens + 
                      Curve1 + I(Curve1^2) + NAT + DistCamp + East +
                      DecWater + DecNatEdge + DecTrail + DecScat + DecPG + 
                      I(DistCamp^2) + DecTrail:NAT,
                    family = binomial, data = Env3)

#specify that, when quadratics are included, linear term is also included
subsetMP <- expression(dc(`DistRoad`, `I(DistRoad^2)`) &
                         dc(`DistCamp`, `I(DistCamp^2)`) &
                         dc(`Curve1`, `I(Curve1^2)`))

options(na.action = "na.fail")
MagpieRes2 <- dredge(MagpieDredge2, subset = subsetMP, trace = 2)

#Get models within 2 AIC
subset(MagpieRes2, delta < 2)

#Make the models
I_EnvMagmod1a <- get.models(MagpieRes2, 1)[[1]]
I_EnvMagmod2a <- get.models(MagpieRes2, 2)[[1]]
I_EnvMagmod3a <- get.models(MagpieRes2, 3)[[1]]
I_EnvMagmod4a <- get.models(MagpieRes2, 4)[[1]]
I_EnvMagmod5a <- get.models(MagpieRes2, 5)[[1]]
I_EnvMagmod6a <- get.models(MagpieRes2, 6)[[1]]
I_EnvMagmod7a <- get.models(MagpieRes2, 7)[[1]]
I_EnvMagmod8a <- get.models(MagpieRes2, 8)[[1]]
I_EnvMagmod9a <- get.models(MagpieRes2, 9)[[1]]
I_EnvMagmod10a <- get.models(MagpieRes2, 10)[[1]]

I_EnvMagmod1a$call
I_EnvMagmod2a$call
I_EnvMagmod3a$call
I_EnvMagmod4a$call
I_EnvMagmod5a$call
I_EnvMagmod6a$call
I_EnvMagmod7a$call
I_EnvMagmod8a$call
I_EnvMagmod9a$call
I_EnvMagmod10a$call


I_EnvMagmod1 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                      DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 1, 
                    family = binomial, data = Env3)

I_EnvMagmod2 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                      DistCamp + I(DistCamp^2) + DistRoad + East + NAT + DecTrail:NAT + 
                      1, family = binomial, data = Env3)

I_EnvMagmod3 <- glm(formula = Magpie2 ~ BldgDens + DecPG + DecScat + DecTrail + 
      DecWater + DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 
      1, family = binomial, data = Env3)

I_EnvMagmod4 <- glm(formula = Magpie2 ~ BldgDens + DecPG + DecScat + DecTrail + 
                      DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                      DecTrail:NAT + 1, family = binomial, data = Env3)

I_EnvMagmod5 <- glm(formula = Magpie2 ~ DecNatEdge + DecPG + DecScat + DecTrail + 
                      DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                      DecTrail:NAT + 1, family = binomial, data = Env3)

I_EnvMagmod6 <- glm(formula = Magpie2 ~ DecNatEdge + DecPG + DecScat + DecTrail + 
                      DecWater + DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 
                      1, family = binomial, data = Env3)

I_EnvMagmod7 <- glm(formula = Magpie2 ~ Curve1 + DecPG + DecScat + DecTrail + 
                      DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                      DecTrail:NAT + 1, family = binomial, data = Env3)

I_EnvMagmod8 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                      DistCamp + I(DistCamp^2) + NAT + DecTrail:NAT + 1, family = binomial, 
                    data = Env3)

I_EnvMagmod9 <- glm(formula = Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                      DistCamp + I(DistCamp^2) + DistRoad + I(DistRoad^2) + East + 
                      NAT + DecTrail:NAT + 1, family = binomial, data = Env3)

I_EnvMagmod10 <- glm(formula = Magpie2 ~ Curve1 + DecPG + DecScat + DecTrail + 
                       DecWater + DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 
                       1, family = binomial, data = Env3)


#Summarise  
summary(I_EnvMagmod1)
summary(I_EnvMagmod2)
summary(I_EnvMagmod3)
summary(I_EnvMagmod4)
summary(I_EnvMagmod5)
summary(I_EnvMagmod6)
summary(I_EnvMagmod7)
summary(I_EnvMagmod8)
summary(I_EnvMagmod9)
summary(I_EnvMagmod10)


#Conduct model averaging
EnvAvg <- model.avg(I_EnvMagmod1, I_EnvMagmod2, I_EnvMagmod3, I_EnvMagmod4,
                    I_EnvMagmod5, I_EnvMagmod6, I_EnvMagmod7, I_EnvMagmod8,
                    I_EnvMagmod9, I_EnvMagmod10)

summary(EnvAvg)
confint(EnvAvg, subset = TRUE)

#Nagelkerke pseudo R2
PseudoR2(I_EnvMagmod1, which = "Nagelkerke") 
PseudoR2(I_EnvMagmod2, which = "Nagelkerke")  
PseudoR2(I_EnvMagmod3, which = "Nagelkerke") 
PseudoR2(I_EnvMagmod4, which = "Nagelkerke") 
PseudoR2(I_EnvMagmod5, which = "Nagelkerke") 
PseudoR2(I_EnvMagmod6, which = "Nagelkerke")  
PseudoR2(I_EnvMagmod7, which = "Nagelkerke") 
PseudoR2(I_EnvMagmod8, which = "Nagelkerke") 
PseudoR2(I_EnvMagmod9, which = "Nagelkerke") 
PseudoR2(I_EnvMagmod10, which = "Nagelkerke")  

mean(c(PseudoR2(I_EnvMagmod1, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod2, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod3, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod4, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod5, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod6, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod7, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod8, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod9, which = "Nagelkerke"),
       PseudoR2(I_EnvMagmod10, which = "Nagelkerke"))) #0.2262512

#McFadden pseudo R2
PseudoR2(I_EnvMagmod1, which = "McFadden") 
PseudoR2(I_EnvMagmod2, which = "McFadden")  
PseudoR2(I_EnvMagmod3, which = "McFadden") 
PseudoR2(I_EnvMagmod4, which = "McFadden") 
PseudoR2(I_EnvMagmod5, which = "McFadden") 
PseudoR2(I_EnvMagmod6, which = "McFadden")  
PseudoR2(I_EnvMagmod7, which = "McFadden") 
PseudoR2(I_EnvMagmod8, which = "McFadden") 
PseudoR2(I_EnvMagmod9, which = "McFadden") 
PseudoR2(I_EnvMagmod10, which = "McFadden")  

mean(c(PseudoR2(I_EnvMagmod1, which = "McFadden"),
       PseudoR2(I_EnvMagmod2, which = "McFadden"),
       PseudoR2(I_EnvMagmod3, which = "McFadden"),
       PseudoR2(I_EnvMagmod4, which = "McFadden"),
       PseudoR2(I_EnvMagmod5, which = "McFadden"),
       PseudoR2(I_EnvMagmod6, which = "McFadden"),
       PseudoR2(I_EnvMagmod7, which = "McFadden"),
       PseudoR2(I_EnvMagmod8, which = "McFadden"),
       PseudoR2(I_EnvMagmod9, which = "McFadden"),
       PseudoR2(I_EnvMagmod10, which = "McFadden"))) #0.1373087

#ROC
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod1))) #0.7256
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod2))) #0.727
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod3))) #0.726
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod4))) #0.725
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod5))) #0.7272
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod6))) #0.7247
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod7))) #0.7268
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod8))) #0.7264
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod9))) #0.727
auc(roc(Env3$Magpie2~fitted(I_EnvMagmod10))) #0.7247

mean(c(auc(roc(Env3$Magpie2~fitted(I_EnvMagmod1))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod2))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod3))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod4))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod5))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod6))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod7))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod8))),
       auc(roc(Env3$Magpie2~fitted(I_EnvMagmod9))),
        auc(roc(Env3$Magpie2~fitted(I_EnvMagmod10))))) #0.7259396

#cross validation
data_ctrl <- trainControl(method = "cv", number = 5)
model_caret1 <- train(Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                        DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 1,
                      family = binomial, data = Env3,
                          trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret1) 

model_caret2 <- train(Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                        DistCamp + I(DistCamp^2) + DistRoad + East + NAT + DecTrail:NAT + 
                        1, family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret2) #



model_caret3 <- train(Magpie2 ~ BldgDens + DecPG + DecScat + DecTrail + 
                        DecWater + DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 
                        1, family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret3) 

model_caret4 <- train(Magpie2 ~ BldgDens + DecPG + DecScat + DecTrail + 
                        DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                        DecTrail:NAT + 1, family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret4) 

model_caret5 <- train(Magpie2 ~ DecNatEdge + DecPG + DecScat + DecTrail + 
                        DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                        DecTrail:NAT + 1, family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret5) 

model_caret6 <- train(Magpie2 ~ DecNatEdge + DecPG + DecScat + DecTrail + 
                        DecWater + DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 
                        1, 
                      family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret6) #

model_caret7 <- train(Magpie2 ~ Curve1 + DecPG + DecScat + DecTrail + 
                        DecWater + DistCamp + I(DistCamp^2) + DistRoad + East + NAT + 
                        DecTrail:NAT + 1, 
                      family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret7) #

model_caret8 <- train(Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                        DistCamp + I(DistCamp^2) + NAT + DecTrail:NAT + 1, 
                      family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret8) #

model_caret9 <- train(Magpie2 ~ DecPG + DecScat + DecTrail + DecWater + 
                        DistCamp + I(DistCamp^2) + DistRoad + I(DistRoad^2) + East + 
                        NAT + DecTrail:NAT + 1, 
                      family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret9) 

model_caret10 <- train(Magpie2 ~ Curve1 + DecPG + DecScat + DecTrail + 
                         DecWater + DistCamp + I(DistCamp^2) + East + NAT + DecTrail:NAT + 
                         1, 
                      family = binomial, data = Env3,
                      trControl = data_ctrl,  method = "glm", na.action = na.pass)
print(model_caret10) 

#get mean accuracy
model_caret1$results[2]
model_caret2$results[2]
model_caret3$results[2]
model_caret4$results[2]
model_caret5$results[2]
model_caret6$results[2]
model_caret7$results[2]
model_caret8$results[2]
model_caret9$results[2]
model_caret10$results[2]


mean(c(0.68564, 0.6886994, 0.6691842, 0.6781394, 0.682617, 0.6870722, 0.6916283,
       0.6871325, 0.6900123, 0.6886477))
#0.6848773

model_caret1$results[3]
model_caret2$results[3]
model_caret3$results[3]
model_caret4$results[3]
model_caret5$results[3]
model_caret6$results[3]
model_caret7$results[3]
model_caret8$results[3]
model_caret9$results[3]
model_caret10$results[3]

mean(c(0.2687731, 0.2789592, 0.236671, 0.2527001, 0.2620211, 0.2687161, 0.2847894,
       0.2745051, 0.2847144, 0.2774048))

#To get OR for model averaged results, repeat with predictors that aren't scaled and centered

S_I_EnvMagmod1 <- glm(formula = Magpie2 ~ DecPGOG + DecScatOG + DecTrailOG + DecWaterOG + 
                      DistCampOG + I(DistCampOG^2) + EastOG + NATOG + DecTrailOG:NATOG + 1, 
                    family = binomial, data = Env)

S_I_EnvMagmod2 <- glm(formula = Magpie2 ~ DecPGOG + DecScatOG + DecTrailOG + DecWaterOG + 
                      DistCampOG + I(DistCampOG^2) + DistRoadOG + EastOG + NATOG + DecTrailOG:NATOG + 
                      1, family = binomial, data = Env)

S_I_EnvMagmod3 <- glm(formula = Magpie2 ~ BldgDensOG + DecPGOG + DecScatOG + DecTrailOG + 
                      DecWaterOG + DistCampOG + I(DistCampOG^2) + EastOG + NATOG + DecTrailOG:NATOG + 
                      1, family = binomial, data = Env)

S_I_EnvMagmod4 <- glm(formula = Magpie2 ~ BldgDensOG + DecPGOG + DecScatOG + DecTrailOG + 
                      DecWaterOG + DistCampOG + I(DistCampOG^2) + DistRoadOG + EastOG + NATOG + 
                      DecTrailOG:NATOG + 1, family = binomial, data = Env)

S_I_EnvMagmod5 <- glm(formula = Magpie2 ~ DecNatEdgeOG + DecPGOG + DecScatOG + DecTrailOG + 
                      DecWaterOG + DistCampOG + I(DistCampOG^2) + DistRoadOG + EastOG + NATOG + 
                      DecTrailOG:NATOG + 1, family = binomial, data = Env)

S_I_EnvMagmod6 <- glm(formula = Magpie2 ~ DecNatEdgeOG + DecPGOG + DecScatOG + DecTrailOG + 
                      DecWaterOG + DistCampOG + I(DistCampOG^2) + EastOG + NATOG + DecTrailOG:NATOG + 
                      1, family = binomial, data = Env)

S_I_EnvMagmod7 <- glm(formula = Magpie2 ~ Curve1OG + DecPGOG + DecScatOG + DecTrailOG + 
                      DecWaterOG + DistCampOG + I(DistCampOG^2) + DistRoadOG + EastOG + NATOG + 
                      DecTrailOG:NATOG + 1, family = binomial, data = Env)

S_I_EnvMagmod8 <- glm(formula = Magpie2 ~ DecPGOG + DecScatOG + DecTrailOG + DecWaterOG + 
                      DistCampOG + I(DistCampOG^2) + NATOG + DecTrailOG:NATOG + 1, family = binomial, 
                    data = Env)

S_I_EnvMagmod9 <- glm(formula = Magpie2 ~ DecPGOG + DecScatOG + DecTrailOG + DecWaterOG + 
                      DistCampOG + I(DistCampOG^2) + DistRoadOG + I(DistRoadOG^2) + EastOG + 
                      NATOG + DecTrailOG:NATOG + 1, family = binomial, data = Env)

S_I_EnvMagmod10 <- glm(formula = Magpie2 ~ Curve1OG + DecPGOG + DecScatOG + DecTrailOG + 
                       DecWaterOG + DistCampOG + I(DistCampOG^2) + EastOG + NATOG + DecTrailOG:NATOG + 
                       1, family = binomial, data = Env)


EnvAvg_Notscaled <- model.avg(S_I_EnvMagmod1, S_I_EnvMagmod2, S_I_EnvMagmod3,
                              S_I_EnvMagmod4, S_I_EnvMagmod5, S_I_EnvMagmod6,
                              S_I_EnvMagmod7, S_I_EnvMagmod8, S_I_EnvMagmod9,
                              S_I_EnvMagmod10)
summary(EnvAvg_Notscaled)
exp(EnvAvg_Notscaled$coefficients)
options(scipen=999)
exp(confint(EnvAvg_Notscaled))

#Make a figure similar to PlotContext3 showing beta coefficients
EnvContext <- plot_summs(I_EnvMagmod1, I_EnvMagmod2, I_EnvMagmod3, I_EnvMagmod4, I_EnvMagmod5,
                         I_EnvMagmod6, I_EnvMagmod7, I_EnvMagmod8, I_EnvMagmod9, I_EnvMagmod10,
                          colors = c("black", "black", "black",
                                     "black", "black", "black",
                                     "black", "black", "black",
                                     "black"),
                          point.size = 30,
                          point.shape = FALSE,
                          coefs = c("Distance to camp" = "DistCamp",
                                    "Decay distance to water" = "DecWater",
                                    "Distance to camp (Squared)" = "I(DistCamp^2)",
                                    "Natural land cover" = "NAT",
                                    "Decay distance to maintained trail" = "DecTrail",
                                    "Decay distance to trail x\nNatural land cover" = "DecTrail:NAT",
                                    "Decay distance to playground" = "DecPG",
                                    "Decay distance to scat" = "DecScat",
                                    "East index" = "East",
                                    "Distance to road" = "DistRoad",
                                    "Building density" = "BldgDens",
                                    "Distance to building" = "DistBldg",
                                    "Decay distance to natural edge" = "DecNatEdge",
                                    "Topogrpahic concavity" = "Curve1",
                                    "Distance to road (Squared)" = "I(DistRoad^2)"))

EnvContext2 <- EnvContext + theme_classic()

EnvContext3 <- EnvContext2 +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Predictor estimate", title = "(B) Environmental predictors") +
  scale_shape_manual(values = c(24,24,24,24,24,24,24,24,24,24,24))

#Create figure of the two coefficient estimate plots beside each other
Plot1 <- ggarrange(PlotContext3, EnvContext3, widths = c(2,3))
ggsave("Fig2.jpg", Plot1, width = 10, height = 5.5, dpi = 700) #arranged horizontally
Plot2 <- ggarrange(PlotContext3, EnvContext3, ncol = 1, heights = c(2,4))
ggsave("CoefficientEstimatesMagpie2.png", Plot2, width = 5, height = 8, dpi = 700) #arranged vertically

#Create photo figure
Pic1 <- image_read("ScatMP.jpg")
Pic2 <- image_read("ScatMP2.jpg")
Pic3 <- image_read("Magpie3.jpeg")
Pic4 <- image_read("Magpie4.jpg")

image_browse(Pic1)
image_browse(Pic2)
image_browse(Pic3)
image_browse(Pic4)

image_info(Pic1)
image_info(Pic2)
image_info(Pic3)
image_info(Pic4)

#Target is images that are 2442 wide, 1831 height, but start by zooming in + cropping
#crop to desired extent
Pic1a <- image_rotate(Pic1, 90)
image_info(Pic1a)
image_browse(Pic1a)
#Rescale to height 1831
Pic1b <- image_scale(Pic1a, "x1831")
image_info(Pic1b)#ALL GOOD 


#Rescale to height 1831
Pic2a <- image_scale(Pic2, "x1831")
image_info(Pic2a) #ALL GOOD

#crop to desired extent
Pic3a <- image_crop(Pic3, "700x500+250+150")
image_browse(Pic3a)
#Rescale to height 1831
Pic3b <- image_scale(Pic3a, "x1831")
image_info(Pic3b)
#Now crop again to desired width as well
Pic3c <- image_crop(Pic3b, "2442x1831+30")
image_browse(Pic3c)
image_info(Pic3c) #ALL GOOD

#crop to desired extent
Pic4a <- image_crop(Pic4, "1500x1000+700+200")
image_browse(Pic4a)
#Rescale to height 1831
Pic4b <- image_scale(Pic4a, "x1831")
image_info(Pic4b)
#Now crop again to desired width as well
Pic4c <- image_crop(Pic4b, "2442x1831+30")
image_browse(Pic4c)
image_info(Pic4c) #ALL GOOD

#Add white borders
Pic1c <- image_border(Pic1b, color = "white", "10x10")
Pic2b <- image_border(Pic2a, color = "white", "10x10")
Pic3d <- image_border(Pic3c, color = "white", "10x10")
Pic4d <- image_border(Pic4c, color = "white", "10x10")

#Pic1b, Pic2a, Pic3c, Pic4c
img1 <- c(Pic1c, Pic2b)
img2 <- c(Pic3d, Pic4d)

PhotoFig1 <- image_append(img1, stack = FALSE)
PhotoFig2 <- image_append(img2, stack = FALSE)

image_browse(PhotoFig1)
image_browse(PhotoFig2)

#NOw add all four photos together
img3 <- c(PhotoFig1, PhotoFig2)

PhotoFig3 <- image_append(img3, stack = TRUE)
image_browse(PhotoFig3)

PhotoFig4 <- image_rotate(PhotoFig3, degrees = 90)
image_browse(PhotoFig4)

#Save two-paneled photo figure
image_write(PhotoFig4, path = "FIG1.jpg", format = "jpg")


#Make figure showing interaction plot + key marginal plots
IntPlot1 <- interact_plot(model = S_I_EnvMagmod1, pred = DecTrailOG,
              modx = NATOG, legend.main = "Natural land cover",
              line.thickness = 2,
              colors = c("darkgrey", "darkgrey", "black"))

IntPlot2 <- IntPlot1 + theme_classic() + ggtitle("(A)")


IntPlot3 <- IntPlot2 +
  theme(legend.position = c(0.30, 0.8),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1),
        legend.title = element_text(colour = "black", face = "plain", size = 12),
        legend.text = element_text(colour = "black", face = "plain", size = 12),
        axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Distance to maintained trail (m)", y = "Likelihood of coprophagy by magpies")
IntPlot3 ###NEED TO CALCULATE X AXIS VALUES AND LIMIT TO REASONABLE LIMITS

#Create marginal plots for variables that are extra important
CampPlot1 <- plot_model(S_I_EnvMagmod1, type = "pred", terms = "DistCampOG")

Campplot2 <- CampPlot1 + theme_classic() + ggtitle("(C)")

CampPlot3 <- Campplot2 +
  theme(axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Distance to camp (m)", y = "Likelihood of coprophagy by magpies") +
  ylim(0,0.5) +xlim(0,7000)
Camplot4 <- CampPlot3 + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))

CampPlot5 <- Camplot4 + geom_hline(yintercept=0.42, linetype="dashed", 
           color = "black", size=1)

ggpredict(EnvMagmod1_OG, terms = "DistCampOG [0, 3200, 6200]")


PGPlot1 <- plot_model(S_I_EnvMagmod1, type = "pred", terms = "DecPGOG")

PGPlot2 <- PGPlot1 + theme_classic() + ggtitle("(B)")

PGPlot3 <- PGPlot2 +
  theme(axis.title.y = element_text(colour = "black", face = "plain", size = 12),
        axis.title.x = element_text(colour = "black", face = "plain", size = 12),
        axis.text.y = element_text(colour = "black", face = "plain", size = 12),
        axis.text.x = element_text(colour = "black", face = "plain", size = 12)) +
  labs(x = "Distance to playground (m)", y = "Likelihood of coprophagy by magpies")
PGPlot3 
PGPlot4 <- PGPlot3 + scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
PGPlot4

ggpredict(EnvMagmod1_OG, terms = "DecPGOG[0, 0.25]")


#Arrange into a composute figure with annotations
CompositeMarginals <- ggarrange(IntPlot3, PGPlot4, CampPlot5, ncol = 1, nrow = 3)
ggsave("CompositeMarginalsDec7.png", CompositeMarginals, width = 4, height = 10, dpi = 700)
