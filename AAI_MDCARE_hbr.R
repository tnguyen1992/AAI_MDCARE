library(readr)
library(glmmTMB)
library(car)
library(emmeans)
library(effects)
library(sjmisc)
setwd("~/ucloud/Projects /MCARE")

MCARE_wtc_hbr <- read_csv("MCARE_wtc_hbr.csv", 
                          col_names = FALSE, na=c("NaN", "NA", "Inf"))
DCARE_wtc_hbr <- read_csv("DCARE_wtc_hbr.csv", 
                          col_names = FALSE, na=c("NaN", "NA", "Inf"))

colnames(MCARE_wtc_hbr) <- list("ID", "cond", "ch", "wtc")
colnames(DCARE_wtc_hbr) <- list("wtc", "ID", "cond","skip", "ch","skip" )

MCARE_wtc_hbr$caregiver <- 1
DCARE_wtc_hbr$caregiver <- 2

DCARE_wtc_hbr$ID <- DCARE_wtc_hbr$ID+90

DCARE_wtc_hbr <- data.frame(DCARE_wtc_hbr$ID, DCARE_wtc_hbr$cond, DCARE_wtc_hbr$ch, DCARE_wtc_hbr$wtc, DCARE_wtc_hbr$caregiver)
colnames(DCARE_wtc_hbr) <- list("ID", "cond", "ch", "wtc", "caregiver")


xdata.wtc.hbr <- rbind(MCARE_wtc_hbr,DCARE_wtc_hbr)



xdata.wtc.hbr$roi1 <- ifelse(xdata.wtc.hbr$ch=="1"|xdata.wtc.hbr$ch=="2"|xdata.wtc.hbr$ch=="3"|xdata.wtc.hbr$ch=="4","ldlpfc",NA)
xdata.wtc.hbr$roi2 <- ifelse(xdata.wtc.hbr$ch=="5"|xdata.wtc.hbr$ch=="6"|xdata.wtc.hbr$ch=="7"|xdata.wtc.hbr$ch=="8","rdlpfc",NA)
xdata.wtc.hbr$roi3 <- ifelse(xdata.wtc.hbr$ch=="9"|xdata.wtc.hbr$ch=="10"|xdata.wtc.hbr$ch=="11"|xdata.wtc.hbr$ch=="12","ltpj",NA)
xdata.wtc.hbr$roi4 <- ifelse(xdata.wtc.hbr$ch=="13"|xdata.wtc.hbr$ch=="14"|xdata.wtc.hbr$ch=="15"|xdata.wtc.hbr$ch=="16","rtpj",NA)

paste5 <- function(..., sep = " ", collapse = NULL, na.rm = F) {
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

xdata.wtc.hbr$roi <- paste5(xdata.wtc.hbr$roi1,xdata.wtc.hbr$roi2,xdata.wtc.hbr$roi3,xdata.wtc.hbr$roi4, na.rm=T)
xdata.wtc.hbr <- xdata.wtc.hbr[,c(1:5,10)]

xdata.wtc.hbr<-xdata.wtc.hbr[!(xdata.wtc.hbr$cond=="4"),]

xdata.wtc.hbr$cond <- as.factor(xdata.wtc.hbr$cond)
xdata.wtc.hbr$cond <- factor(xdata.wtc.hbr$cond, levels=c(1,2,3), labels=c("coop", "ind", "rest"))
xdata.wtc.hbr$roi <- as.factor(xdata.wtc.hbr$roi)
xdata.wtc.hbr$caregiver <- factor(xdata.wtc.hbr$caregiver, levels=c(1,2), labels=c("mother", "father"))
##############################################################################
# #rpa analysis
rpa_wtc_hbo <- read_csv("rpa_wtc_hbo.csv",
                        col_names = FALSE)
colnames(rpa_wtc_hbo) <- list("wtc", "ID", "cond", "skip", "ch", "roi")
rpa_wtc_hbo$caregiver <- 1
rpa_wtc_hbo <- rpa_wtc_hbo[,c(2,3,5,1,7,6)]
rpa_wtc_hbo$roi <- factor(rpa_wtc_hbo$roi, levels=c(1,2,3,4),labels=c("ldlpfc", "rdlpfc",
                                                    "ltpj", "rtpj"))

# 
# 
rpa_df <- read_csv("~/ucloud/Projects /MCARE/DCARE_wtc_hbo_rpa.csv",
                              col_names = FALSE)
colnames(rpa_df) <- c("ID", "cond", "ch", "wtc")
rpa_df$caregiver <- 2
# 
rpa_df$ID <- rpa_df$ID+90
rpa_df<-rpa_df[!(rpa_df$cond=="4"),]
rpa_df$roi1 <- ifelse(rpa_df$ch=="1"|rpa_df$ch=="2"|rpa_df$ch=="3"|rpa_df$ch=="4","ldlpfc",NA)
rpa_df$roi2 <- ifelse(rpa_df$ch=="5"|rpa_df$ch=="6"|rpa_df$ch=="7"|rpa_df$ch=="8","rdlpfc",NA)
rpa_df$roi3 <- ifelse(rpa_df$ch=="9"|rpa_df$ch=="10"|rpa_df$ch=="11"|rpa_df$ch=="12","ltpj",NA)
rpa_df$roi4 <- ifelse(rpa_df$ch=="13"|rpa_df$ch=="14"|rpa_df$ch=="15"|rpa_df$ch=="16","rtpj",NA)
rpa_df$roi <- paste5(rpa_df$roi1,rpa_df$roi2,rpa_df$roi3,rpa_df$roi4, na.rm=T)
rpa_df <- rpa_df[,c(1:5,10)]
# 
DCARE_ID <- data.frame(DCARE_wtc_hbo$ID)
colnames(DCARE_ID) <- list("ID")
DCARE_ID <- unique(DCARE_ID)
rpa_df <- merge(DCARE_ID,rpa_df, by="ID", all.x=T )

MDCARE_rpa <- rbind(rpa_wtc_hbo,rpa_df)
MDCARE_rpa$c.sex <- factor(MDCARE_rpa$caregiver, levels=c(1,2), labels=c("mother-child dyads", "father-child dyads"))
MDCARE_rpa$random <- 1
MDCARE_rpa$caregiver <- factor(MDCARE_rpa$caregiver, levels=c(1,2), labels=c("mother", "father"))
MDCARE_rpa$cond <- factor(MDCARE_rpa$cond, levels=c(1,2,3), labels=c("coop", "ind", "rest"))

# 
xdata.wtc.hbr$random <- 0
# 
finaldata <- rbind(xdata.wtc.hbr,MDCARE_rpa)
# 
# finaldata <- xdata.wtc.hbr
# 
finaldata$random <- factor(finaldata$random, levels=c(0,1), labels=c("true", "random"))
# 
finaldata$wtc[finaldata$wtc<=0]<-NA
# 
rpa.df <- merge(MDCARE.sample, finaldata, by="ID", all.x=T)
# 
rpa.df <- rpa.df[complete.cases(rpa.df), ]
rpa.df<-rpa.df[rpa.df$cond != 'rest', ] 

# model
#rpa.df$wtc <- rpa.df$wtc+0.0001
m1 <- lmer(wtc ~ cond*random*roi+(1|ID),
              data=rpa.df,
              #family=beta_family(link="logit"),
              na.action=na.omit)

Anova(m1)
summary(m1)
plot(effect("random:roi", m1))
emmeans(m1, pairwise ~ random|cond, adjust="tukey", type="response")



##########################################################
# condition comparison
# exclude resting phase
xdata.wtc.hbr<-xdata.wtc.hbr[xdata.wtc.hbr$cond != 'rest', ] 
library(sjPlot)
xdata.wtc.hbr$c.sex <- factor(xdata.wtc.hbr$caregiver, labels=c("mother-child dyads", "father-child dyads"))

xdata.wtc.hbr$wtc <- xdata.wtc.hbr$wtc+0.0001
cond.m <- glmmTMB(wtc ~ cond  * roi * c.sex +(1|ID), 
                  data=xdata.wtc.hbr, 
                  family=beta_family(link="logit"), 
                  na.action=na.exclude)
Anova(cond.m)
summary(cond.m)
plot_model(cond.m, type="eff", terms=c("cond"))#,# colors="gs",
#            show.data=F, jitter=20, title="")
library(ggplot2)
ggplot(data=xdata.wtc.hbr,aes(x=cond,y=wtc, color=cond))+
  geom_jitter(width=0.1)+
  geom_boxplot(aes(alpha=0.3))+
  #geom_smooth(method='lm',formula=y~x,color="black")+
  xlab("conditions")+ 
  ylab("wtc")+
  facet_grid(c.sex~roi)+
  theme_bw(base_size=14)+
  theme(legend.position='none')
#ylim(0,1)
emmeans(cond.m, pairwise ~ cond|roi|c.sex, type="logit")
confint(emmeans(cond.m, pairwise ~ cond|roi|caregiver, type="logit"))

# plot_model(cond.m, type="eff", terms=c("cond"), colors="gs", 
#            show.data=F, jitter=1, title="")
# emmeans(cond.m, pairwise ~ cond|c.sex, type="log")


# additional contrasts (coop - ind)
# extract coop
xdata.wtc.hbr.coop <- subset(xdata.wtc.hbr, cond == "coop")
xdata.wtc.hbr.ind <- subset(xdata.wtc.hbr, cond == "ind")
xdata.wtc.hbr.ind$wtc.ind  <- xdata.wtc.hbr.ind$wtc 
xdata.wtc.hbr.ind <- xdata.wtc.hbr.ind[,c(1:3,9)]
xdata.wtc.c <- merge(xdata.wtc.hbr.coop,xdata.wtc.hbr.ind, by=c("ID", "ch"), all.x=T)
xdata.wtc.c$wtc.cont <- xdata.wtc.c$wtc - xdata.wtc.c$wtc.ind

library(lme4)
cond.m2 <- lmer(wtc.cont ~ roi * c.sex +(1|ID), 
                data=xdata.wtc.c, 
                #family=beta_family(link="logit"), 
                na.action=na.exclude)
Anova(cond.m2)
summary(cond.m2)
plot(allEffects(cond.m2))
emmeans(cond.m2, pairwise ~ c.sex|roi, type="log", correct='Tukey')
##############################################################################
# hypothesis 2
library(readxl)
AAI_data <- read_excel("AAI_data.xlsx", 
                       skip = 1)
AAI_data$ID <- ifelse(AAI_data$caregiver==2,AAI_data$ID+90,AAI_data$ID) 
AAI_data$caregiver <- factor(AAI_data$caregiver, levels=c(1,2), labels=c("mother", "father"))
xdata.AAI <- merge(AAI_data, xdata.wtc.hbr, by=c("ID", "caregiver"), all=T)

xdata.AAI$secure <- factor(xdata.AAI$secure, levels=c(0,1), labels=c("insecure", "secure"))
xdata.AAI$savx <- factor(xdata.AAI$savx, levels=c(1,2,3), labels=c("secure", "AV", "AX"))
xdata.AAI<-xdata.AAI[xdata.AAI$cond != "rest", ]
xdata.AAI<-xdata.AAI[xdata.AAI$cond != "ind", ]
xdata.AAI <- xdata.AAI[complete.cases(c(xdata.AAI$secure)), ]

library(DescTools)


ID<-unique(xdata.AAI$ID)
ID<-data.frame(ID)
desc.AAI <- merge(AAI_data, ID, by=c("ID"), all.y=T)
desc.AAI <- desc.AAI[complete.cases(c(desc.AAI$secure)), ]

MDCARE.sample<-as.data.frame(ID)

xdata.AAI$attachment <- factor(xdata.AAI$secure)
Desc(xdata.AAI$attachment, na.rm=T)

aai.m <- glmmTMB(wtc ~  roi  * attachment * c.sex +(1|ID), 
                 data=xdata.AAI, 
                 family=beta_family(link="logit"), 
                 na.action=na.exclude)
Anova(aai.m)
summary(aai.m)
plot_model(aai.m, type="eff", terms=c("roi","attachment", "c.sex"), colors="gs", title="")
plot(allEffects(aai.m))
emmeans(aai.m, pairwise ~attachment*c.sex|roi,type='response')
emmeans(aai.m, pairwise ~c.sex,type='response')

plot(effect("roi:secure:caregiver",aai.m))


ggplot(data=xdata.AAI,aes(x=attachment,y=wtc, color=attachment))+
  geom_jitter(width=0.1)+
  geom_boxplot(aes(alpha=0.3))+
  #geom_smooth(method='lm',formula=y~x,color="black")+
  xlab("caregiver attachment representation")+ 
  ylab("wtc")+
  facet_grid(c.sex~roi)+
  theme_bw(base_size=14)+
  theme(legend.position='none')



# additional contrasts (coop - ind)
# extract coop
xdata.AAI.coop <- subset(xdata.AAI, cond == "coop")
xdata.AAI.ind <- subset(xdata.AAI, cond == "ind")
xdata.AAI.ind$wtc.ind  <- xdata.AAI.ind$wtc 
xdata.AAI.ind <- xdata.AAI.ind[,c(1,9,10)]
xdata.AAI.c <- merge(xdata.AAI.coop,xdata.AAI.ind, by=c("ID", "ch"), all.x=T)
xdata.AAI.c$wtc.cont <- xdata.AAI.c$wtc.x - xdata.AAI.c$wtc.y

library(lme4)
aai.m1 <- glmmTMB(wtc.cont ~  roi  * attachment * c.sex +(1|ID), 
                  data=xdata.AAI.c, 
                  family=beta_family(link="logit"), 
                  na.action=na.exclude)
Anova(aai.m1)
summary(aai.m1)
emmeans(aai.m1,pairwise ~c.sex|roi|attachment)




aai.m <- glmmTMB(wtc ~  roi  * savx * caregiver+(1|ID), 
                 data=xdata.AAI, 
                 family=beta_family(link="logit"), 
                 na.action=na.exclude)
Anova(aai.m)
emmeans(aai.m, pairwise ~ savx|caregiver|roi)
plot(effect("roi:savx:caregiver",aai.m))


aai.m <- glmmTMB(wtc ~  roi  * coherence * caregiver+(1|ID), 
                 data=xdata.AAI, 
                 family=beta_family(link="logit"), 
                 na.action=na.exclude)
Anova(aai.m)
plot_model(aai.m, type="eff", terms=c("coherence","roi", "c.sex"), colors="simply")
emtrends(aai.m,  pairwise ~ roi|caregiver, var="coherence")

##############
# sensitivity
library(haven)
int.xdata <- read_sav('M-CARE_D-CARE_InteraktionsdatenERLANGEN.sav')
int.xdata$ID<- int.xdata$dyad
int.xdata$caregiver<- int.xdata$parent

int.xdata$ID <- ifelse(int.xdata$caregiver==2,int.xdata$ID+90,int.xdata$ID) 
int.xdata$caregiver <- factor(int.xdata$caregiver, levels=c(1,2), labels=c("mother", "father"))

int.wtc.xdata <- merge(xdata.wtc.hbr, int.xdata, by=c("ID", "caregiver"))
int.wtc.xdata.coop<-int.wtc.xdata[int.wtc.xdata$cond == "coop", ]

int.wtc.xdata.coop <- merge(MDCARE.sample, int.wtc.xdata.coop, by="ID", all.x=T)


int.m <- glmmTMB(wtc ~  roi * caregiver * supportive_presence_puzzle +(1|ID), 
                 data=int.wtc.xdata.coop, 
                 family=beta_family(link="logit"), 
                 na.action=na.exclude, 
                 control=glmmTMBControl(optCtrl = list(iter.max = 300, eval.max = 400)))
Anova(int.m)
summary(int.m)
plot_model(int.m, type="eff", terms=c("supportive_presence_puzzle", "caregiver", "roi"), colors="simply")
plot(effect("caregiver:supportive_presence_puzzle",int.m))
emtrends(int.m, pairwise ~ caregiver|roi, var="supportive_presence_puzzle")



########## reciprocity
library(pscl)

rec.xdata <- read_delim("MDCARE_rec.csv", 
                        delim = ";", escape_double = FALSE, col_types = cols(...3 = col_skip(), 
                                                                             ...4 = col_skip(), ...5 = col_skip(), 
                                                                             ...6 = col_skip(), ...7 = col_skip()), 
                        locale = locale(decimal_mark = ",", 
                                        grouping_mark = "."), trim_ws = TRUE)
rec.xdata$rec <- rec.xdata$reciprocity
rec.xdata <- merge(rec.xdata,ID, by="ID", all.y=T)
pdata1 <- merge(int.xdata,rec.xdata, by="ID" , all.x=T)



rec.m1 <- glmmTMB(rec ~ supportive_presence_puzzle * caregiver, 
                  data=pdata1)
summary(rec.m1)
Anova(rec.m1)
plot(allEffects(rec.m1))
emtrends(rec.m1, ~caregiver, var="supportive_presence_puzzle")

pdata2 <- merge(AAI_data,rec.xdata, by="ID" , all.x=T)
pdata2.m<-pdata2[pdata2$caregiver == "mother", ]
pdata2.f<-pdata2[pdata2$caregiver == "father", ]

pdata2$attachment <- factor(pdata2$secure, levels=c(0,1), labels=c("secure", "insecure")) 
pdata2$attachment <- factor(pdata2$savx, levels=c(1,2,3), labels=c("secure", "preoccupied", "avoidant")) 

hist(pdata2.m$rec)
hist(pdata2.f$rec)

rec.m2 <- glmmTMB(rec ~ attachment * caregiver , 
                  data=pdata2,
                  family=gaussian())
summary(rec.m2)
Anova(rec.m2)
emmeans(rec.m2, pairwise ~ caregiver)
plot(allEffects(rec.m2))
plot_model(rec.m2, type="pred", terms=c( "caregiver"))

pdata3 <- merge(xdata.wtc.hbr,rec.xdata, by="ID" , all.x=T)
pdata3$rec.z <- scale(log(pdata3$reciprocity), center=T, scale=T)
pdata3.c<-pdata3[pdata3$cond == "coop", ]
pdata3.m<-pdata3.c[pdata3.c$caregiver == "father", ]

recm <- glmmTMB(wtc ~  roi * reciprocity+(1|ID), 
                data=pdata3.m, 
                family=beta_family(link="logit"), 
                na.action=na.exclude, 
                control=glmmTMBControl(optCtrl = list(iter.max = 300, eval.max = 400)))
Anova(recm)
summary(recm)
emtrends(recm,  ~roi|caregiver , var="reciprocity")
plot_model(recm, type="eff", terms=c( "reciprocity", "roi"))

hist(pdata3.c$rec.z)
hist(pdata3.c$reciprocity)


pdata4 <- merge(storystem,rec.xdata, by=c("ID") , all.x=T)
pdata4<- merge(pdata4,demo, by="ID")


rec.m4 <- glmmTMB(rec ~ M_COH4 * sex.child , 
                  data=pdata4)
summary(rec.m4)
plot(allEffects(rec.m4))

##########
# Parental reflective functioning
library(haven)
ment.xdata <- read_sav("~/ucloud/Projects /MCARE/M-CARE_D-CARE_PRFQ_MK021121.sav")

ment.xdata$ID<- ment.xdata$dyad
ment.xdata$caregiver<- ment.xdata$parent

ment.xdata$ID <- ifelse(ment.xdata$caregiver==2,ment.xdata$ID+90,ment.xdata$ID) 
ment.xdata$caregiver <- factor(ment.xdata$caregiver, levels=c(1,2), labels=c("mother", "father"))

ment.wtc.xdata <- merge(xdata.wtc.hbr, ment.xdata, by=c("ID", "caregiver"))
ment.wtc.xdata.coop<-ment.wtc.xdata[ment.wtc.xdata$cond == "coop", ]

ment.wtc.xdata.coop <- merge(MDCARE.sample, ment.wtc.xdata.coop, by="ID", all.x=T)


ment.m <- glmmTMB(wtc ~  roi * caregiver * (PRFQ_CM+PRFQ_IC+PRFQ_PM) +(1|ID), 
                  data=ment.wtc.xdata.coop, 
                  family=beta_family(link="logit"), 
                  na.action=na.exclude, 
                  control=glmmTMBControl(optCtrl = list(iter.max = 300, eval.max = 400)))
Anova(ment.m)
summary(ment.m)
plot_model(ment.m, type="eff", terms=c("PRFQ_IC", "caregiver"), colors="simply")
plot(effect("PRFQ_PM",ment.m))
emtrends(ment.m, pairwise ~ caregiver, var="PRFQ_IC")


# story stems
library(readxl)
storystem <- read_excel("ECR_STEMS_and_AAI_FINAL_2022.xlsx", 
                        sheet = "ECR,_AAI_and_STEMS")

#library(haven)
#storystem <- read_sav("MD_CARE_DATA_AAI_STORY_STEMS_28-05-22.sav")

#storystem$ID<- storystem$ID_CHECK
#storystem$caregiver<- storystem$D_or_MCARE

storystem$ID <- ifelse(storystem$caregiver==2,storystem$ID+90,storystem$ID) 
storystem$caregiver <- factor(storystem$caregiver, levels=c(1,2), labels=c("mother", "father"))

storystem <- merge(ID, storystem, by=c("ID"), all.x=T)


stem.wtc.xdata <- merge(xdata.wtc.hbr, storystem, by=c("ID", "caregiver"))
stem.wtc.xdata.coop<-stem.wtc.xdata[stem.wtc.xdata$cond == "coop", ]

stem.wtc.xdata.coop$M_COH_DICH <- as.factor(stem.wtc.xdata.coop$M_COH_DICH)

library(readxl)
D_CARE_demographic <- read_excel("cs-transfer/D_CARE_demographic.xlsx")
dcare.demo <- D_CARE_demographic[,c(1,4)]
colnames(dcare.demo)<- list("ID","sex.child" )
dcare.demo$ID <- dcare.demo$ID+90

M_CARE_demographic_TN <- read_excel("cs-transfer/M_CARE_demographic_TN.xlsx")
mcare.demo <- M_CARE_demographic_TN[,c(1,4)]
colnames(mcare.demo)<- list("ID","sex.child" )

demo <- rbind(mcare.demo, dcare.demo)

stem.wtc.xdata.coop<- merge(stem.wtc.xdata.coop,demo, by="ID")
storystem<- merge(storystem,demo, by="ID")

library(doBy)
summary_df <- summaryBy(M_COH4 ~ sex.child, data = stem.wtc.xdata.coop, FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
summary(lm(M_COH4 ~ sex.child, data = stem.wtc.xdata.coop))

stem.wtc.xdata.coop$attachment <- factor(stem.wtc.xdata.coop$SvsISclass, levels=c(0,1), labels=c("secure", "insecure"))
stem.wtc.xdata.coop.m<-stem.wtc.xdata.coop[stem.wtc.xdata.coop$caregiver == "mother", ]

stem.m <- glmmTMB(wtc ~  roi*M_COH4* attachment  +(1|ID), 
                  data=stem.wtc.xdata.coop, 
                  family=beta_family(link="logit"), 
                  na.action=na.exclude, 
                  control=glmmTMBControl(optCtrl = list(iter.max = 300, eval.max = 400)))
Anova(stem.m)
summary(stem.m)
plot(allEffects(stem.m))
emtrends(stem.m, pairwise ~attachment, var="M_COH4")
emtrends(stem.m, pairwise ~sex.child|roi, var="M_AVD4")

stem.wtc.xdata.coop$sex.child <- factor(stem.wtc.xdata.coop$sex.child, labels=c("parent-daughter dyads", "parent-son dyads"))

ggplot(data=stem.wtc.xdata.coop,aes(x=M_COH4,y=wtc))+
  geom_jitter(width=0.1)+
  #geom_boxplot(aes(alpha=0.3))+
  geom_smooth(method='lm',formula=y~x)+
  xlab("child attachment coherence")+ 
  ylab("wtc")+
  facet_grid(sex.child~roi)+
  theme_bw(base_size=14)+
  theme(legend.position='none')






