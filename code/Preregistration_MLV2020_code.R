## ----setup,include=FALSE,cache=FALSE,echo=FALSE--------------------------

library(MASS)
library(knitr)
library(xtable)
library(tidyverse)
library(gridExtra)
library(lme4)
library(reshape2)


## ----echo=FALSE----------------------------------------------------------
critt<-qt((0.05/6)/2,df=39)
critz<-qnorm((0.05/6)/2)


## ----dillonloaddata, include=TRUE, echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
# Bayesian analysis of original Experiment 1 from Dillon et al 2013
orig <- read.table("../data/originaldataExp1/data_experiment1_jml.txt",header=TRUE)
orig$subj <- factor(orig$subj)
orig$item <- factor(orig$item)
orig$cond <- factor(orig$cond)
# rename conditions to match our data
levels(orig$cond) <- c('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h')
# cast data to wide format to have the same data structure as in the analysis of our data above
orig <- dcast(orig, subj+item+cond+region ~ fixationtype, value.var = "value")

orig <- orig[, c('subj', 'item', 'cond', 'region', 'ff','fp', 'pr',  'rp', 'rr', 'tt')]
# rename columns:
colnames(orig) <- c('subj', 'item', 'cond', 'roi', 'FFD','FPRT', 'FPR', 'RPD','RRT','TFT')


## ----dillonvars, include=TRUE, echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
Nsubj_orig <- length(unique(factor(orig$subj)))
Nitem_orig <- length(unique(factor(subset(orig, cond!='filler')$item)))


## ----dilloncontrasts, include=TRUE,  echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
# Condition Labels (using our labels; same contrasts as above in analysis of our data); for documentation of contrasts see above. 

orig$Dep <- ifelse(orig$cond %in% c('a', 'b', 'c', 'd'), .5, -.5) # main effect of dependency type: agr=0.5, refl=-0.5
orig$Gram <- ifelse(orig$cond %in% c('a', 'b', 'e', 'f'), -.5, .5) # main effect of grammaticality: gram=-.5, ungram=.5
orig$Int_gram <- ifelse(orig$cond %in% c('a','e'), .5, ifelse(orig$cond %in% c('b', 'f'), -.5, 0) ) # interference in grammatical sentences: distr-match=0.5, distr-mismatch=-0.5
orig$Int_ungram <- ifelse(orig$cond %in% c('d', 'h'), .5, ifelse(orig$cond %in% c('c', 'g'), -.5, 0)) # interference in ungrammatical sentences: distr-match=0.5, distr-mismatch=-0.5
orig$DepxInt_gram <-ifelse(orig$cond %in% c('a', 'f'), .5, ifelse(orig$cond %in% c('b', 'e'), -.5, 0))
orig$DepxInt_ungram <- ifelse(orig$cond %in% c('d', 'g'), .5, ifelse(orig$cond %in% c('c', 'h'), -.5, 0))
orig$DepxGram <- ifelse(orig$cond %in% c('c', 'd', 'e', 'f'), .5, -.5)

#orig$DepxInt <- ifelse(orig$cond %in% c('a', 'd', 'f', 'g'), 0.5, -0.5)
#orig$Int <- ifelse(orig$cond %in% c('a', 'd', 'e', 'h'), 0.5, -0.5) # main effect of interference: int=0.5, no int=-0.5
#orig$GramxInt <- ifelse(orig$cond %in% c('b', 'd', 'f', 'h'),0.5,-0.5)
#orig$DepxGramxInt <- ifelse(orig$cond %in% c('b', 'd', 'e', 'g'), 0.5, -0.5)

orig$Int_gram_refl <- ifelse(orig$cond %in% c('e'), .5, ifelse(orig$cond %in% c('f'), -.5, 0))
orig$Int_gram_agr <- ifelse(orig$cond %in% c('a'), .5, ifelse(orig$cond %in% c('b'), -.5, 0))
orig$Int_ungram_refl <- ifelse(orig$cond %in% c('h'), .5, ifelse(orig$cond %in% c('g'), -.5, 0))
orig$Int_ungram_agr <- ifelse(orig$cond %in% c('d'), .5, ifelse(orig$cond %in% c('c'), -.5, 0))


## ----dilloncrit, include=TRUE,echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
# in original data, critical region is 
crit_orig <- subset(orig, roi==5 & cond!='filler') 


## ----analysesDillon,eval=FALSE,echo=FALSE,cache=TRUE---------------------
## library(lme4)
## mFFD <- lmer(FFD~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit_orig,control=lmerControl(calc.derivs=FALSE))
## resFFD<-summary(mFFD)$coefficients[8,]
## 
## mFPRT<-lmer(FPRT~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit_orig,control=lmerControl(calc.derivs=FALSE))
## resFPRT<-summary(mFPRT)$coefficients[8,]
## 
## mFPR<-glmer(FPR~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit_orig,family=binomial())
## resFPR<-summary(mFPR)$coefficients[8,1:3]
## 
## mRPD<-lmer(RPD~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit_orig,control=lmerControl(calc.derivs=FALSE))
## resRPD<-summary(mRPD)$coefficients[8,]
## 
## mRRT<-lmer(RRT~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit_orig,control=lmerControl(calc.derivs=FALSE))
## resRRT<-summary(mRRT)$coefficients[8,]
## 
## mTFT<-lmer(TFT~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit_orig,control=lmerControl(calc.derivs=FALSE))
## 
## mTFTmax<-lmer(log(TFT+1)~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram ||subj)+(1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram ||item),crit_orig,control=lmerControl(calc.derivs=FALSE))
## resTFT<-summary(mTFT)$coefficients[8,]
## resTFTmax<-summary(mTFTmax)
## 
## res<-rbind(resFFD,resFPRT,resFPR,resRPD,resRRT,resTFT)
## res<-data.frame(DepVar=c("FFD","FPRT","FPR","RPD","RRT","TFT"),res)


## ----xtableres,echo=FALSE,eval=FALSE-------------------------------------
## xtable(res,include.rownames=FALSE)


## ----loaddata, include=TRUE, cache=TRUE, echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
d1 <- read.table("../data/dataDillonRepV1.txt",header=TRUE)
# delete subject 201t (same subject as 201t_2)
d1 <- subset(d1, d1$subject!='201t')
# rename subject 168bry and 201t_2
d1$subject <- ifelse(d1$subject=='168bry', '168b', as.character(d1$subject))
d1$subject <- ifelse(d1$subject=='201t_2', '201b', as.character(d1$subject))
d1$subject <- factor(d1$subject)
# edit subject ids such that they match with the working memory data
d1$subject <- unlist(strsplit(as.character(d1$subject), split=c('b')))
# remove 0s at id onset
d1$subject<-as.factor(as.integer(d1$subject))

d2 <- read.table("../data/dataDillonRepV2.txt",header=TRUE)
#summary(factor(d1$RESPONSE_ACCURACY)) 
#summary(factor(d2$RESPONSE_ACCURACY)) # 1=correct; 0=incorrect; -1=trial without a question
#xtabs(~RESPONSE_ACCURACY+item, d2)
# items followed by comprehension question
q <- sort(unique(subset(d2, RESPONSE_ACCURACY!=-1)$item))
# code absence of question as -1 in RESPONSE_ACCURACY of d1 (analoguously to coding in d2)
d1$RESPONSE_ACCURACY <- ifelse(!d1$item %in% q, -1, d1$RESPONSE_ACCURACY)

# Correct some wrong subject id labels in dataset V2 (=d2)
## 321_dr2 is the correct file for participant 321; 321_dr3 was incorrectly labeled and should have bene 386_dr3.
d2$subject <- factor(ifelse(d2$subject=='321_dr3', '386_dr3', as.character(d2$subject)))
##337_dr1 is the correct file for participant 337; 337_dr2 was incorrectly labeled and should be 338_dr1.
d2$subject <- factor(ifelse(d2$subject=='337_dr2', '338_dr1', as.character(d2$subject)))
# The data for participant 334 (from V2) is labeled as 223_dr6 => should be 334_dr6.
d2$subject <- factor(ifelse(d2$subject=='223_dr6', '334_dr6', as.character(d2$subject)))
# Data for participant 331 is labeled as 362_dr3 => should be 331_dr3.
d2$subject <- factor(ifelse(d2$subject=='2362_dr3', '331_dr3', as.character(d2$subject)))


d2$subject <- unlist(strsplit(as.character(d2$subject), split='[_dlr].*'))
d2$subject <- factor(d2$subject)


d <- rbind(d1,d2)
#colnames(d)
d$FPR <- as.numeric(as.logical(d$RBRC)) # first-pass regression
d$REG <- as.numeric(as.logical(d$TRC)) # whether or not there was an regression
d <- d[, c('subject', 'item', 'condition', 'RESPONSE_ACCURACY', 'roi', 'FFD','FPRT',  'FPR', 'RPD', 'RRT','TFT')]
colnames(d) <- c('subj', 'item', 'cond', 'acc', 'roi', 'FFD','FPRT',  'FPR', 'RPD', 'RRT','TFT')
d$item <- factor(d$item)
d$subj <- factor(d$subj)

#write.table(d, file='data/dataJMVV.txt', sep = '\t')
#d <- read.table(file='data/dataJMVV.txt', sep = '\t')


## ----vars, include=TRUE, cache=TRUE, echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
Nsubj <- length(unique(factor(d$subj)))
#Nsubj_noWM <- length(which(is.na(op$span))) # number of subjects from which no wm data recorded
Nitem <- length(unique(factor(subset(d, cond!='filler')$item)))
Nfiller <- length(unique(factor(subset(d, cond=='filler')$item)))


## ----accs, include=TRUE, cache=TRUE, echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
dq <- subset(d, item %in% q & roi==1) # comprehension accuracies
acc <- round(tapply(dq$acc, dq$cond, mean), 2)
# remove fillers
d <- subset(d, cond!='filler')


## ----contrasts, include=TRUE, cache=TRUE,  echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
# Condition Labels (our labels are different from the ones used by Dillon et al 2013)
# a - d Agreement conditions.  
# a. grammatical, interference 
# b. grammatical, no interference
# c. ungrammatical, no interference
# d. ungrammatical, interference

# e - h Reflexives conditions. 
# e. grammatical, interference 
# f. grammatical, no interference 
# g. ungrammatical, no interference
# h. ungrammatical, interference


# Model 1: 
# main effect of dependency type (dep); positive effect => agr > refl
# main effect of grammaticality (gram);  positive effect => ungram > gram
# dep x gram
# interference within gram (int_gram); positive effect => inhib interference
# interference within ungram (int_ungram); => positive effect => inhib interference
# int_gram x dep; 
# int_ungram x dep; 
#### Only int_gram, int_ungram, int_gram x dep, int_ungram x dep are of theoretical interest


# Model 2: resolving interaction between dep and int_gram and between dep and int_ungram 
# main effect of dependency type (dep);  
# main effect of grammaticality (gram); 
# dep x gram
# interference within gram  in reflexives (int_gram_refl); positive effect => inhib interference
# interference within gram  in reflexives (int_gram_agr); positive effect => inhib interference
# interference within ungram in reflexives (int_ungram_refl); positive effect => inhib interference
# interference within ungram in agreement (int_ungram_agr); positive effect => inhib interference
#### Only int_gram_refl, int_gram_agr,  int_ungram_refl, int_ungram_agr are of theoretical interest

# Model 1wm: same as Model 1, with the following additional effects:
# main effect of working memory
# interactions between working memory and all other fixed effects (including interactions)

# Model 2wm: same as Model 2, with the following additional effects:
# main effect of working memory
# interactions between working memory and all other fixed effects (including interaction)


d$Dep <- ifelse(d$cond %in% c('a', 'b', 'c', 'd'), .5, -.5) # main effect of dependency type: agr=0.5, refl=-0.5
d$Gram <- ifelse(d$cond %in% c('a', 'b', 'e', 'f'), -.5, .5) # main effect of grammaticality: gram=-.5, ungram=.5
d$Int_gram <- ifelse(d$cond %in% c('a','e'), .5, ifelse(d$cond %in% c('b', 'f'), -.5, 0) ) # interference in grammatical sentences: distr-match=0.5, distr-mismatch=-0.5
d$Int_ungram <- ifelse(d$cond %in% c('d', 'h'), .5, ifelse(d$cond %in% c('c', 'g'), -.5, 0)) # interference in ungrammatical sentences: distr-match=0.5, distr-mismatch=-0.5
d$DepxInt_gram <-ifelse(d$cond %in% c('a', 'f'), .5, ifelse(d$cond %in% c('b', 'e'), -.5, 0))
d$DepxInt_ungram <- ifelse(d$cond %in% c('d', 'g'), .5, ifelse(d$cond %in% c('c', 'h'), -.5, 0))
d$DepxGram <- ifelse(d$cond %in% c('c', 'd', 'e', 'f'), .5, -.5)

d$DepxInt <- ifelse(d$cond %in% c('a', 'd', 'f', 'g'), 0.5, -0.5)
d$Int <- ifelse(d$cond %in% c('a', 'd', 'e', 'h'), 0.5, -0.5) # main effect of interference: int=0.5, no int=-0.5
d$GramxInt <- ifelse(d$cond %in% c('b', 'd', 'f', 'h'), 0.5, -0.5)
d$DepxGramxInt <- ifelse(d$cond %in% c('b', 'd', 'e', 'g'), 0.5, -0.5)

d$Int_gram_refl <- ifelse(d$cond %in% c('e'), .5, ifelse(d$cond %in% c('f'), -.5, 0))
d$Int_gram_agr <- ifelse(d$cond %in% c('a'), .5, ifelse(d$cond %in% c('b'), -.5, 0))
d$Int_ungram_refl <- ifelse(d$cond %in% c('h'), .5, ifelse(d$cond %in% c('g'), -.5, 0))
d$Int_ungram_agr <- ifelse(d$cond %in% c('d'), .5, ifelse(d$cond %in% c('c'), -.5, 0))


## ----crit, include=TRUE, cache=TRUE, echo=FALSE, warning=TRUE, error=TRUE, message=TRUE----
crit <- subset(d, roi==12 & cond!='filler') 


## ----analysesDillonrep,eval=FALSE,echo=FALSE,cache=TRUE------------------
## mFFD <- lmer(FFD~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit,control=lmerControl(calc.derivs=FALSE))
## resFFD<-summary(mFFD)$coefficients[8,]
## 
## mFPRT<-lmer(FPRT~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit,control=lmerControl(calc.derivs=FALSE))
## resFPRT<-summary(mFPRT)$coefficients[8,]
## 
## 
## mFPR<-glmer(FPR~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit,family=binomial())
## resFPR<-summary(mFPR)$coefficients[8,1:3]
## 
## 
## 
## mRPD<-lmer(RPD~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit,control=lmerControl(calc.derivs=FALSE))
## resRPD<-summary(mRPD)$coefficients[8,]
## 
## mRRT<-lmer(RRT~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit,control=lmerControl(calc.derivs=FALSE))
## resRRT<-summary(mRRT)$coefficients[8,]
## 
## mTFT<-lmer(TFT~ 1+Dep+Gram+DepxGram+Int_gram+Int_ungram+DepxInt_gram+DepxInt_ungram+ (1+DepxInt_ungram ||subj)+(1+DepxInt_ungram ||item),crit,control=lmerControl(calc.derivs=FALSE))
## 
## resTFT<-summary(mTFT)$coefficients[8,]
## 
## res<-rbind(resFFD,resFPRT,resFPR,resRPD,resRRT,resTFT)
## res<-data.frame(DepVar=c("FFD","FPRT","FPR","RPD","RRT","TFT"),res)


## ----xtableresrep,echo=FALSE,eval=FALSE----------------------------------
## xtable(res,include.rownames=FALSE)

