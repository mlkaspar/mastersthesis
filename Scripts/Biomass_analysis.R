require(readr)
require(reshape2)
require(ggplot2)
require(scales)
require(dplyr)
require(gridExtra)
##########################################################
##########################################################
# Initial -----------------------------------------------------------------
#functions
stderror = function(x){
  sd(x)/sqrt(length(x))
}
scientific_e10 <- function(x) {
  parse(text=gsub("e", "%*%10^", scientific_format()(x)))
}

#change number of decimal places for more correct calculations:
options(digits=14)

colorsstr1 = c("gold", "gold3", "darkorange", "darkorange3", "firebrick1", "firebrick3")
colorsstr2 = c("cadetblue1", "cadetblue3", "steelblue1", "steelblue3", "royalblue1", "royalblue3")
clrs = c(D="azure4", L = "azure")
clrs2 = c("darkgoldenrod1","darkorange", "firebrick")
clrs2_2 = c("cadetblue1", "steelblue1", "royalblue2")

reds = c("darkorange","darkorange", "#FF3333", "#FF3333", "#9F3620", "#9F3620")
blues = c("cadetblue2", "cadetblue2", "#33BBFF", "#33BBFF", "#000FFF", "#000FFF")
colrs = append(reds, blues)
##########################################################
##########################################################
#load tables with biomass data
strain1_d3 = read.table(choose.files(), header=TRUE, sep=";")
strain2_d3 = read.table(choose.files(), header=TRUE, sep=";")

strain1_d11 = read.table(choose.files(), header=TRUE, sep=";")
strain2_d11 = read.table(choose.files(), header=TRUE, sep=";")

strain1_d18 = read.table(choose.files(), header=TRUE, sep=";")
strain2_d18 = read.table(choose.files(), header=TRUE, sep=";")

#load table with cell numbers from flowcytometer
Alltripl_table = read.table(choose.files(), header=TRUE, sep=";")

#backups
strain1_d3_BU = strain1_d3
strain2_d3_BU = strain2_d3

strain1_d11_BU = strain1_d11
strain2_d11_BU = strain2_d11

strain1_d18_BU = strain1_d18
strain2_d18_BU = strain2_d18
##########################################################
##########################################################
# Calculations and statistics for each of the three sampling days -----------------------------------------------------------------

day = "d3"

#full join: join 2 tables by whatever column (here they are named the same, but you can also do by = c("variable1" = "column1"))
# %>% is like a pipe. Take the result from full_join and filter them directly, and save in variable "output"
output = full_join(strain1_d3, Alltripl_table, by = c("Name" = "Name")) %>% filter(Samplingday==day)
output = output[!grepl(2, output$Name),] #delete data for strain2_d3

strain1_d3 = output

strain1_d3$biomass = as.double(as.character(strain1_d3$Carbon)) * (strain1_d3$Cellnumber/1e09) #get ug/ul = g/l

#95% interval
strain1_d3_95 = subset(strain1_d3, strain1_d3$biomass > quantile(strain1_d3$biomass, 0.025))
strain1_d3_95 = subset(strain1_d3, strain1_d3$biomass < quantile(strain1_d3$biomass, 0.975))

#keep order of names and ID variables
strain1_d3_95$Name_f = factor(strain1_d3_95$Name, levels = c("DL1A","DL1B", "DL1C","LL1A", "LL1B", "LL1C", "DM1A","DM1B","DM1C", "LM1A", "LM1B", "LM1C", "DH1A","DH1B","DH1C",  "LH1A", "LH1B", "LH1C"))
strain1_d3_95$ID_f = factor(strain1_d3_95$ID, levels = c("DL1", "LL1", "DM1", "LM1", "DH1", "LH1"))

strain1_d3_95$biomasslog = log(strain1_d3_95$biomass)

#aggregate is awesome! makes a "summarized" data frame based on whatever function you want, in this case "mean".

strain1_d3_means = aggregate(strain1_d3_95[, 16], list(strain1_d3_95$Name), mean) #mean of log data
strain1_d3_STEs = aggregate(strain1_d3_95[, 16], list(strain1_d3_95$Name), stderror) #stderror is a function defined at the to

#renaming and adding columns
strain1_d3_means = as.data.frame(strain1_d3_means)
colnames(strain1_d3_means) = c("Name", "Biomass")
strain1_d3_means$ID = factor(gsub('.{1}$', '', strain1_d3_means$Name)) # replace last 1 character with empty
strain1_d3_means$Light = factor(gsub('.{3}$', '', strain1_d3_means$Name)) # replace last 3 character with empty
strain1_d3_means$Medium = factor(gsub('^{1}.|.{2}$', '', strain1_d3_means$Name)) #replace first and two last

#subsetting
strain1_d3_means$Medium_f = factor(strain1_d3_means$Medium, levels=c("L", "M", "H"))
strain1_d3.lowc = subset(strain1_d3_means, Medium=="L")
strain1_d3.mediumc = subset(strain1_d3_means, Medium=="M")
strain1_d3.highc = subset(strain1_d3_means, Medium=="H")

strain1_d3.light = subset(strain1_d3_means, Light == "L")
strain1_d3.dark = subset(strain1_d3_means, Light == "D")
##########################################################
#save R output with sink(). commands are not exported! you can use print() for comments inbetween
sink(file="strain1_d3_log.txt")
strain1_d3.Lo1 = lm(data=strain1_d3.lowc, Biomass ~ Light)
summary(strain1_d3.Lo1)
anova(strain1_d3.Lo1)
confint(strain1_d3.Lo1)
a1 = aov(strain1_d3.Lo1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d3.lowc, Biomass ~ Light)

strain1_d3.M1 = lm(data=strain1_d3.mediumc, Biomass ~ Light)
summary(strain1_d3.M1)
anova(strain1_d3.M1)
confint(strain1_d3.M1)
a1 = aov(strain1_d3.M1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d3.mediumc, Biomass ~ Light)

strain1_d3.H1 = lm(data=strain1_d3.highc, Biomass ~ Light)
summary(strain1_d3.H1)
anova(strain1_d3.H1)
confint(strain1_d3.H1)
a1 = aov(strain1_d3.H1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d3.highc, Biomass ~ Light)
############
strain1_d3.Li1 = lm(data=strain1_d3.light, Biomass ~ Medium)
summary(strain1_d3.Li1)
anova(strain1_d3.Li1)
confint(strain1_d3.Li1)
aL1 = aov(strain1_d3.Li1)
TukeyHSD(x=aL1, "Medium", conf.level=0.95)

strain1_d3.D1 = lm(data=strain1_d3.dark, Biomass ~ Medium)
summary(strain1_d3.D1)
anova(strain1_d3.D1)
confint(strain1_d3.D1)
aD1 = aov(strain1_d3.D1)
TukeyHSD(x=aD1, "Medium", conf.level=0.95)
sink()
####################################################################################################################
#strain 2
strain2_d3 = as.data.frame(strain2_d3)

output = full_join(strain2_d3, Alltripl_table, by = c("Name" = "Name")) %>% filter(Samplingday==day)
output = output[!grepl(1, output$Name),] 

strain2_d3 = output

strain2_d3$biomass = as.double(as.character(strain2_d3$Carbon)) * (strain2_d3$Cellnumber/1e09) 

strain2_d3_95 = subset(strain2_d3, strain2_d3$biomass > quantile(strain2_d3$biomass, 0.025))
strain2_d3_95 = subset(strain2_d3, strain2_d3$biomass < quantile(strain2_d3$biomass, 0.975))

strain2_d3_95$Name_f = factor(strain2_d3_95$Name, levels = c("DL2A","DL2B", "DL2C","LL2A", "LL2B", "LL2C", "DM2A","DM2B","DM2C", "LM2A", "LM2B", "LM2C", "DH2A","DH2B","DH2C",  "LH2A", "LH2B", "LH2C"))
strain2_d3_95$ID_f = factor(strain2_d3_95$ID, levels = c("DL2", "LL2", "DM2", "LM2", "DH2", "LH2"))

strain2_d3_95$biomasslog = log(strain2_d3_95$biomass)
strain2_d3_means = aggregate(strain2_d3_95[, 16], list(strain2_d3_95$Name), mean) 
strain2_d3_STEs = aggregate(strain2_d3_95[, 16], list(strain2_d3_95$Name), stderror)

strain2_d3_means = as.data.frame(strain2_d3_means)
colnames(strain2_d3_means) = c("Name", "Biomass")
strain2_d3_means$ID = factor(gsub('.{1}$', '', strain2_d3_means$Name)) 
strain2_d3_means$Light = factor(gsub('.{3}$', '', strain2_d3_means$Name))
strain2_d3_means$Medium = factor(gsub('^{1}.|.{2}$', '', strain2_d3_means$Name))

strain2_d3_means$Medium_f = factor(strain2_d3_means$Medium, levels=c("L", "M", "H"))
strain2_d3.lowc = subset(strain2_d3_means, Medium=="L")
strain2_d3.mediumc = subset(strain2_d3_means, Medium=="M")
strain2_d3.highc = subset(strain2_d3_means, Medium=="H")

strain2_d3.light = subset(strain2_d3_means, Light == "L")
strain2_d3.dark = subset(strain2_d3_means, Light == "D")
##########################################################
sink(file="strain2_d3_log.txt")

strain2_d3.Lo1 = lm(data=strain2_d3.lowc, Biomass ~ Light)
summary(strain2_d3.Lo1)
anova(strain2_d3.Lo1)
confint(strain2_d3.Lo1)
a1 = aov(strain2_d3.Lo1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d3.lowc, Biomass ~ Light)

strain2_d3.M1 = lm(data=strain2_d3.mediumc, Biomass ~ Light)
summary(strain2_d3.M1)
anova(strain2_d3.M1)
confint(strain2_d3.M1)
a1 = aov(strain2_d3.M1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d3.mediumc, Biomass ~ Light)

strain2_d3.H1 = lm(data=strain2_d3.highc, Biomass ~ Light)
summary(strain2_d3.H1)
anova(strain2_d3.H1)
confint(strain2_d3.H1)
a1 = aov(strain2_d3.H1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d3.highc, Biomass ~ Light)
############
strain2_d3.Li1 = lm(data=strain2_d3.light, Biomass ~ Medium)
summary(strain2_d3.Li1)
anova(strain2_d3.Li1)
confint(strain2_d3.Li1)
aL1 = aov(strain2_d3.Li1)
TukeyHSD(x=aL1, "Medium", conf.level=0.95)

strain2_d3.D1 = lm(data=strain2_d3.dark, Biomass ~ Medium)
summary(strain2_d3.D1)
anova(strain2_d3.D1)
confint(strain2_d3.D1)
aD1 = aov(strain2_d3.D1)
TukeyHSD(x=aD1, "Medium", conf.level=0.95)
sink()
##########################################################
##########################################################
##########################################################
##########################################################

day = "d11"

strain1_d11 = lapply(strain1_d11, function(x){
  gsub("-two", "", x)
})
strain1_d11 = as.data.frame(strain1_d11)

output = full_join(strain1_d11, Alltripl_table, by = c("Name" = "Name")) %>% filter(Samplingday==day)
output = output[!grepl(2, output$Name),] #delete data for strain2_d11

strain1_d11 = output

strain1_d11$biomass = as.double(as.character(strain1_d11$Carbon)) * (strain1_d11$Cellnumber/1e09) #get ug/ul = g/l

strain1_d11_95 = subset(strain1_d11, strain1_d11$biomass > quantile(strain1_d11$biomass, 0.025))
strain1_d11_95 = subset(strain1_d11, strain1_d11$biomass < quantile(strain1_d11$biomass, 0.975))

strain1_d11_95$Name_f = factor(strain1_d11_95$Name, levels = c("DL1A","DL1B", "DL1C","LL1A", "LL1B", "LL1C", "DM1A","DM1B","DM1C", "LM1A", "LM1B", "LM1C", "DH1A","DH1B","DH1C",  "LH1A", "LH1B", "LH1C"))
strain1_d11_95$ID_f = factor(strain1_d11_95$ID, levels = c("DL1", "LL1", "DM1", "LM1", "DH1", "LH1"))

strain1_d11_95$biomasslog = log(strain1_d11_95$biomass)
strain1_d11_means = aggregate(strain1_d11_95[, 16], list(strain1_d11_95$Name), mean) 
strain1_d11_STEs = aggregate(strain1_d11_95[, 16], list(strain1_d11_95$Name), stderror)

strain1_d11_means = as.data.frame(strain1_d11_means)
colnames(strain1_d11_means) = c("Name", "Biomass")
strain1_d11_means$ID = factor(gsub('.{1}$', '', strain1_d11_means$Name)) 
strain1_d11_means$Light = factor(gsub('.{3}$', '', strain1_d11_means$Name))
strain1_d11_means$Medium = factor(gsub('^{1}.|.{2}$', '', strain1_d11_means$Name))

strain1_d11_means$Medium_f = factor(strain1_d11_means$Medium, levels=c("L", "M", "H"))
strain1_d11.lowc = subset(strain1_d11_means, Medium=="L")
strain1_d11.mediumc = subset(strain1_d11_means, Medium=="M")
strain1_d11.highc = subset(strain1_d11_means, Medium=="H")

strain1_d11.light = subset(strain1_d11_means, Light == "L")
strain1_d11.dark = subset(strain1_d11_means, Light == "D")
##########################################################
sink(file="strain1_d11_log.txt")
strain1_d11.Lo1 = lm(data=strain1_d11.lowc, Biomass ~ Light)
summary(strain1_d11.Lo1)
anova(strain1_d11.Lo1)
confint(strain1_d11.Lo1)
a1 = aov(strain1_d11.Lo1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d11.lowc, Biomass ~ Light)

strain1_d11.M1 = lm(data=strain1_d11.mediumc, Biomass ~ Light)
summary(strain1_d11.M1)
anova(strain1_d11.M1)
confint(strain1_d11.M1)
a1 = aov(strain1_d11.M1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d11.mediumc, Biomass ~ Light)

strain1_d11.H1 = lm(data=strain1_d11.highc, Biomass ~ Light)
summary(strain1_d11.H1)
anova(strain1_d11.H1)
confint(strain1_d11.H1)
a1 = aov(strain1_d11.H1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d11.highc, Biomass ~ Light)
############
strain1_d11.Li1 = lm(data=strain1_d11.light, Biomass ~ Medium)
summary(strain1_d11.Li1)
anova(strain1_d11.Li1)
confint(strain1_d11.Li1)
aL1 = aov(strain1_d11.Li1)
TukeyHSD(x=aL1, "Medium", conf.level=0.95)

strain1_d11.D1 = lm(data=strain1_d11.dark, Biomass ~ Medium)
summary(strain1_d11.D1)
anova(strain1_d11.D1)
confint(strain1_d11.D1)
aD1 = aov(strain1_d11.D1)
TukeyHSD(x=aD1, "Medium", conf.level=0.95)
sink()
####################################################################################################################

strain2_d11 = as.data.frame(strain2_d11)

output = full_join(strain2_d11, Alltripl_table, by = c("Name" = "Name")) %>% filter(Samplingday==day)
output = output[!grepl(1, output$Name),] #delete data for strain1

strain2_d11 = output

strain2_d11$biomass = as.double(as.character(strain2_d11$Carbon)) * (strain2_d11$Cellnumber/1e09) 

strain2_d11 = strain2_d11[complete.cases(strain2_d11),] #remove NAs

strain2_d11_95 = subset(strain2_d11, strain2_d11$biomass > quantile(strain2_d11$biomass, 0.025))
strain2_d11_95 = subset(strain2_d11, strain2_d11$biomass < quantile(strain2_d11$biomass, 0.975))

strain2_d11_95$Name_f = factor(strain2_d11_95$Name, levels = c("DL2A","DL2B", "DL2C","LL2A", "LL2B", "LL2C", "DM2A","DM2B","DM2C", "LM2A", "LM2B", "LM2C", "DH2A","DH2B","DH2C",  "LH2A", "LH2B", "LH2C"))
strain2_d11_95$ID_f = factor(strain2_d11_95$ID, levels = c("DL2", "LL2", "DM2", "LM2", "DH2", "LH2"))

strain2_d11_95$biomasslog = log(strain2_d11_95$biomass)
strain2_d11_means = aggregate(strain2_d11_95[, 16], list(strain2_d11_95$Name), mean) 
strain2_d11_STEs = aggregate(strain2_d11_95[, 16], list(strain2_d11_95$Name), stderror)

strain2_d11_means = as.data.frame(strain2_d11_means)
colnames(strain2_d11_means) = c("Name", "Biomass")
strain2_d11_means$ID = factor(gsub('.{1}$', '', strain2_d11_means$Name)) 
strain2_d11_means$Light = factor(gsub('.{3}$', '', strain2_d11_means$Name))
strain2_d11_means$Medium = factor(gsub('^{1}.|.{2}$', '', strain2_d11_means$Name))

strain2_d11_means$Medium_f = factor(strain2_d11_means$Medium, levels=c("L", "M", "H"))
strain2_d11.lowc = subset(strain2_d11_means, Medium=="L")
strain2_d11.mediumc = subset(strain2_d11_means, Medium=="M")
strain2_d11.highc = subset(strain2_d11_means, Medium=="H")

strain2_d11.light = subset(strain2_d11_means, Light == "L")
strain2_d11.dark = subset(strain2_d11_means, Light == "D")
##########################################################
sink(file="strain2_d11_t.txt")

strain2_d11.Lo1 = lm(data=strain2_d11.lowc, Biomass ~ Light)
summary(strain2_d11.Lo1)
anova(strain2_d11.Lo1)
confint(strain2_d11.Lo1)
a1 = aov(strain2_d11.Lo1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d11.lowc, Biomass ~ Light)

strain2_d11.M1 = lm(data=strain2_d11.mediumc, Biomass ~ Light)
summary(strain2_d11.M1)
anova(strain2_d11.M1)
confint(strain2_d11.M1)
a1 = aov(strain2_d11.M1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d11.mediumc, Biomass ~ Light)

strain2_d11.H1 = lm(data=strain2_d11.highc, Biomass ~ Light)
summary(strain2_d11.H1)
anova(strain2_d11.H1)
confint(strain2_d11.H1)
a1 = aov(strain2_d11.H1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d11.highc, Biomass ~ Light)
############
strain2_d11.Li1 = lm(data=strain2_d11.light, Biomass ~ Medium)
summary(strain2_d11.Li1)
anova(strain2_d11.Li1)
confint(strain2_d11.Li1)
aL1 = aov(strain2_d11.Li1)
TukeyHSD(x=aL1, "Medium", conf.level=0.95)

strain2_d11.D1 = lm(data=strain2_d11.dark, Biomass ~ Medium)
summary(strain2_d11.D1)
anova(strain2_d11.D1)
confint(strain2_d11.D1)
aD1 = aov(strain2_d11.D1)
TukeyHSD(x=aD1, "Medium", conf.level=0.95)
sink()
##########################################################
##########################################################
##########################################################
##########################################################

day = "d18"

strain1_d18 = lapply(strain1_d18, function(x){
  gsub("-two", "", x)
})
strain1_d18 = as.data.frame(strain1_d18)

output = full_join(strain1_d18, Alltripl_table, by = c("Name" = "Name")) %>% filter(Samplingday==day)
output = output[!grepl(2, output$Name),] #delete data for strain2_d18

strain1_d18 = output

strain1_d18$biomass = as.double(as.character(strain1_d18$Carbon)) * (strain1_d18$Cellnumber/1e09) 

strain1_d18_95 = subset(strain1_d18, strain1_d18$biomass > quantile(strain1_d18$biomass, 0.025))
strain1_d18_95 = subset(strain1_d18, strain1_d18$biomass < quantile(strain1_d18$biomass, 0.975))

strain1_d18_95$Name_f = factor(strain1_d18_95$Name, levels = c("DL1A","DL1B", "DL1C","LL1A", "LL1B", "LL1C", "DM1A","DM1B","DM1C", "LM1A", "LM1B", "LM1C", "DH1A","DH1B","DH1C",  "LH1A", "LH1B", "LH1C"))
strain1_d18_95$ID_f = factor(strain1_d18_95$ID, levels = c("DL1", "LL1", "DM1", "LM1", "DH1", "LH1"))

strain1_d18_95$biomasslog = log(strain1_d18_95$biomass)
strain1_d18_means = aggregate(strain1_d18_95[, 16], list(strain1_d18_95$Name), mean)
strain1_d18_STEs = aggregate(strain1_d18_95[, 16], list(strain1_d18_95$Name), stderror)

strain1_d18_means = as.data.frame(strain1_d18_means)
colnames(strain1_d18_means) = c("Name", "Biomass")
strain1_d18_means$ID = factor(gsub('.{1}$', '', strain1_d18_means$Name)) 
strain1_d18_means$Light = factor(gsub('.{3}$', '', strain1_d18_means$Name))
strain1_d18_means$Medium = factor(gsub('^{1}.|.{2}$', '', strain1_d18_means$Name))

strain1_d18_means$Medium_f = factor(strain1_d18_means$Medium, levels=c("L", "M", "H"))
strain1_d18.lowc = subset(strain1_d18_means, Medium=="L")
strain1_d18.mediumc = subset(strain1_d18_means, Medium=="M")
strain1_d18.highc = subset(strain1_d18_means, Medium=="H")

strain1_d18.light = subset(strain1_d18_means, Light == "L")
strain1_d18.dark = subset(strain1_d18_means, Light == "D")
##########################################################
sink(file="strain1_d18_log.txt")
strain1_d18.Lo1 = lm(data=strain1_d18.lowc, Biomass ~ Light)
summary(strain1_d18.Lo1)
anova(strain1_d18.Lo1)
confint(strain1_d18.Lo1)
a1 = aov(strain1_d18.Lo1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d18.lowc, Biomass ~ Light)

strain1_d18.M1 = lm(data=strain1_d18.mediumc, Biomass ~ Light)
summary(strain1_d18.M1)
anova(strain1_d18.M1)
confint(strain1_d18.M1)
a1 = aov(strain1_d18.M1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d18.mediumc, Biomass ~ Light)

strain1_d18.H1 = lm(data=strain1_d18.highc, Biomass ~ Light)
summary(strain1_d18.H1)
anova(strain1_d18.H1)
confint(strain1_d18.H1)
a1 = aov(strain1_d18.H1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain1_d18.highc, Biomass ~ Light)
############
strain1_d18.Li1 = lm(data=strain1_d18.light, Biomass ~ Medium)
summary(strain1_d18.Li1)
anova(strain1_d18.Li1)
confint(strain1_d18.Li1)
aL1 = aov(strain1_d18.Li1)
TukeyHSD(x=aL1, "Medium", conf.level=0.95)

strain1_d18.D1 = lm(data=strain1_d18.dark, Biomass ~ Medium)
summary(strain1_d18.D1)
anova(strain1_d18.D1)
confint(strain1_d18.D1)
aD1 = aov(strain1_d18.D1)
TukeyHSD(x=aD1, "Medium", conf.level=0.95)
sink()
####################################################################################################################
strain2_d18 = lapply(strain2_d18, function(x){
  gsub("-two", "", x)
})
strain2_d18 = as.data.frame(strain2_d18)

output = full_join(strain2_d18, Alltripl_table, by = c("Name" = "Name")) %>% filter(Samplingday==day)
output = output[!grepl(1, output$Name),] #delete data for strain1

strain2_d18 = output

strain2_d18$biomass = as.double(as.character(strain2_d18$Carbon)) * (strain2_d18$Cellnumber/1e09) 

strain2_d18 = strain2_d18[complete.cases(strain2_d18),] #remove NAs

strain2_d18_95 = subset(strain2_d18, strain2_d18$biomass > quantile(strain2_d18$biomass, 0.025))
strain2_d18_95 = subset(strain2_d18, strain2_d18$biomass < quantile(strain2_d18$biomass, 0.975))

strain2_d18_95$Name_f = factor(strain2_d18_95$Name, levels = c("DL2A","DL2B", "DL2C","LL2A", "LL2B", "LL2C", "DM2A","DM2B","DM2C", "LM2A", "LM2B", "LM2C", "DH2A","DH2B","DH2C",  "LH2A", "LH2B", "LH2C"))
strain2_d18_95$ID_f = factor(strain2_d18_95$ID, levels = c("DL2", "LL2", "DM2", "LM2", "DH2", "LH2"))

strain2_d18_95$biomasslog = log(strain2_d18_95$biomass)
strain2_d18_means = aggregate(strain2_d18_95[, 16], list(strain2_d18_95$Name), mean) 
strain2_d18_STEs = aggregate(strain2_d18_95[, 16], list(strain2_d18_95$Name), stderror)

strain2_d18_means = as.data.frame(strain2_d18_means)
colnames(strain2_d18_means) = c("Name", "Biomass")
strain2_d18_means$ID = factor(gsub('.{1}$', '', strain2_d18_means$Name)) 
strain2_d18_means$Light = factor(gsub('.{3}$', '', strain2_d18_means$Name))
strain2_d18_means$Medium = factor(gsub('^{1}.|.{2}$', '', strain2_d18_means$Name))

strain2_d18_means$Medium_f = factor(strain2_d18_means$Medium, levels=c("L", "M", "H"))
strain2_d18.lowc = subset(strain2_d18_means, Medium=="L")
strain2_d18.mediumc = subset(strain2_d18_means, Medium=="M")
strain2_d18.highc = subset(strain2_d18_means, Medium=="H")

strain2_d18.light = subset(strain2_d18_means, Light == "L")
strain2_d18.dark = subset(strain2_d18_means, Light == "D")
##########################################################
sink(file="strain2_d18_log.txt")

strain2_d18.Lo1 = lm(data=strain2_d18.lowc, Biomass ~ Light)
summary(strain2_d18.Lo1)
anova(strain2_d18.Lo1)
confint(strain2_d18.Lo1)
a1 = aov(strain2_d18.Lo1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d18.lowc, Biomass ~ Light)

strain2_d18.M1 = lm(data=strain2_d18.mediumc, Biomass ~ Light)
summary(strain2_d18.M1)
anova(strain2_d18.M1)
confint(strain2_d18.M1)
a1 = aov(strain2_d18.M1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d18.mediumc, Biomass ~ Light)

strain2_d18.H1 = lm(data=strain2_d18.highc, Biomass ~ Light)
summary(strain2_d18.H1)
anova(strain2_d18.H1)
confint(strain2_d18.H1)
a1 = aov(strain2_d18.H1)
TukeyHSD(x=a1, "Light", conf.level=0.95)
t.test(data=strain2_d18.highc, Biomass ~ Light)
############
strain2_d18.Li1 = lm(data=strain2_d18.light, Biomass ~ Medium)
summary(strain2_d18.Li1)
anova(strain2_d18.Li1)
confint(strain2_d18.Li1)
aL1 = aov(strain2_d18.Li1)
TukeyHSD(x=aL1, "Medium", conf.level=0.95)

strain2_d18.D1 = lm(data=strain2_d18.dark, Biomass ~ Medium)
summary(strain2_d18.D1)
anova(strain2_d18.D1)
confint(strain2_d18.D1)
aD1 = aov(strain2_d18.D1)
TukeyHSD(x=aD1, "Medium", conf.level=0.95)
sink()

##########################################################
##########################################################
##########################################################
##########################################################
# Graphs -----------------------------------------------------------------

#make greek letters in R: expression(). paste is used to glue together various elements, in this case text strings and an expression
ugml = expression(paste("Biomass ", mu, "g/ml"))

d3_L1_plot = ggplot(strain1_d3.light, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab(ugml) + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

d3_D1_plot = ggplot(strain1_d3.dark, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab(ugml) + xlab("") +
  scale_y_continuous(limits=c(0,0.7))
#####################################d3_strain2
d3_L2_plot = ggplot(strain2_d3.light, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2_2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab(ugml) + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

d3_D2_plot = ggplot(strain2_d3.dark, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2_2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab(ugml) + xlab("") +
  scale_y_continuous(limits=c(0,0.7))
######################################
######################################d11
d11_L1_plot = ggplot(strain1_d11.light, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

d11_D1_plot = ggplot(strain1_d11.dark, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))
#####################################d11_strain2
d11_L2_plot = ggplot(strain2_d11.light, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2_2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

d11_D2_plot = ggplot(strain2_d11.dark, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2_2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

######################################
######################################d18
d18_L1_plot = ggplot(strain1_d18.light, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

d18_D1_plot = ggplot(strain1_d18.dark, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))
#####################################d18_strain2
d18_L2_plot = ggplot(strain2_d18.light, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2_2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

d18_D2_plot = ggplot(strain2_d18.dark, aes(x=Medium_f, y=(Biomass))) +  
  geom_boxplot(aes(fill=Medium_f)) + 
  scale_fill_manual(values=clrs2_2) + 
  guides(fill=FALSE) +
  ggtitle("") +
  theme_classic() +
  ylab("") + xlab("") +
  scale_y_continuous(limits=c(0,0.7))

#################################################################################
#################################################################################
#arrange several plots together in one:
totalgraph = grid.arrange(d3_D1_plot, d11_D1_plot, d18_D1_plot, d3_L1_plot, d11_L1_plot, d18_L1_plot, d3_D2_plot, d11_D2_plot, d18_D2_plot, d3_L2_plot, d11_L2_plot, d18_L2_plot, ncol = 6, widths=c(1.1, 1, 1, 1.1, 1, 1))
#first one (d3_D1) has y axis label -> not needed for the other ones per row
##########################################################
##########################################################
##########################################################
##########################################################

# Statistics 2 -----------------------------------------------------------------


s1d3_meanbyname = aggregate(strain1_d3_95[,"biomass"], list(strain1_d3_95$Name), mean)
s1d3_meanbyname$ID = factor(gsub('.{1}$', '', s1d3_meanbyname$Group.1)) #replace last character with empty

s1d3_meanbygroup = aggregate(s1d3_meanbyname[,"x"], list(s1d3_meanbyname$ID), function(x) c(mean = mean(x), sd = stderror(x)))

s2d3_meanbyname = aggregate(strain2_d3_95[,"biomass"], list(strain2_d3_95$Name), mean)
s2d3_meanbyname$ID = factor(gsub('.{1}$', '', s2d3_meanbyname$Group.1)) #replace last character with empty

#you can also summarize by two functions (mean and sd) on same group:
s2d3_meanbygroup = aggregate(s2d3_meanbyname[,"x"], list(s2d3_meanbyname$ID), function(x) c(mean = mean(x), sd = stderror(x)))
#remember: stderror is a selfwritten function at the top
#####################

s1d11_meanbyname = aggregate(strain1_d11_95[,"biomass"], list(strain1_d11_95$Name), mean)
s1d11_meanbyname$ID = factor(gsub('.{1}$', '', s1d11_meanbyname$Group.1)) #replace last character with empty

s1d11_meanbygroup = aggregate(s1d11_meanbyname[,"x"], list(s1d11_meanbyname$ID), function(x) c(mean = mean(x), sd = stderror(x)))

s2d11_meanbyname = aggregate(strain2_d11_95[,"biomass"], list(strain2_d11_95$Name), mean)
s2d11_meanbyname$ID = factor(gsub('.{1}$', '', s2d11_meanbyname$Group.1)) #replace last character with empty

s2d11_meanbygroup = aggregate(s2d11_meanbyname[,"x"], list(s2d11_meanbyname$ID), function(x) c(mean = mean(x), sd = stderror(x)))
######################

s1d18_meanbyname = aggregate(strain1_d18_95[,"biomass"], list(strain1_d18_95$Name), mean)
s1d18_meanbyname$ID = factor(gsub('.{1}$', '', s1d18_meanbyname$Group.1)) #replace last character with empty

s1d18_meanbygroup = aggregate(s1d18_meanbyname[,"x"], list(s1d18_meanbyname$ID), function(x) c(mean = mean(x), sd = stderror(x)))

s2d18_meanbyname = aggregate(strain2_d18_95[,"biomass"], list(strain2_d18_95$Name), mean)
s2d18_meanbyname$ID = factor(gsub('.{1}$', '', s2d18_meanbyname$Group.1)) #replace last character with empty

s2d18_meanbygroup = aggregate(s2d18_meanbyname[,"x"], list(s2d18_meanbyname$ID), function(x) c(mean = mean(x), sd = stderror(x)))
#####################
sink(file="mean_vals_biomass.txt")
print("day 3")
s1d3_meanbygroup
s2d3_meanbygroup

print("day 11")
s1d11_meanbygroup
s2d11_meanbygroup

print("day 18")
s1d18_meanbygroup
s2d18_meanbygroup
sink()

# Carbon uptake -----------------------------------------------------------------

#copied values to Excel, calculated uptake by start concentration

carbon_uptake_df = read.table(choose.files(), header=TRUE, sep=";")
carbon_uptake_df$Sday_f = factor(carbon_uptake_df$Samplingday, levels=c("d3", "d11","d18"))
carbon_uptake_df$Name = factor(carbon_uptake_df$Name, levels=c("DL1", "LL1", "DM1", "LM1", "DH1", "LH1", "DL2", "LL2", "DM2", "LM2", "DH2", "LH2"))

colorsstr1 = c("goldenrod2", "gold", "darkorange3", "darkorange", "firebrick3", "firebrick1")
colorsstr2 = c("cadetblue3", "cadetblue1", "dodgerblue3", "deepskyblue1", "steelblue4", "steelblue2")
anothercolorvector = c("darkorange", "darkorange", "firebrick", "firebrick",  "red3", "red3", "dodgerblue", "dodgerblue", "steelblue", "steelblue", "#000FFF", "#000FFF")

#subset without "Low"
no_Low_subset$Sday_f = factor(no_Low_subset$Samplingday, levels=c("d3", "d11", "d18"))

ggplot(no_Low_subset, aes(x=Sday_f, y=Uptake, group = Name)) + 
  geom_point(aes(color=Name, shape=Light), size= 4, position=position_dodge(width=0.4)) +
  scale_color_manual("Names", values = c("darkorange", "darkorange", "firebrick", "firebrick", "dodgerblue", "dodgerblue", "steelblue", "steelblue")) +
  geom_errorbar(aes(ymin=Uptake-SE, ymax=Uptake+SE),color="gray30", width=0.2, size=0.7, position=position_dodge(width=0.4)) +
  facet_grid(Strain ~ .) +
  xlab("Sampling day") + ylab("Carbon uptake (%)") +
  theme_bw()

#subset with  only "Low"
low_only_subset = subset(carbon_uptake_df, Medium == "L")
low_only_subset$Name = factor(low_only_subset$Name, levels=c("DL1", "LL1", "DL2", "LL2"))
low_only_subset$Sday_f = factor(low_only_subset$Samplingday, levels=c("d3", "d11", "d18"))

ggplot(low_only_subset, aes(x=Day, y=Uptake, group = Name, linetype=Light)) + 
  geom_line(aes(colour=Name), size=1) +
  scale_color_manual(values = c("red3", "red3", "#000FFF", "#000FFF")) +
  geom_errorbar(limits, width=0.03, linetype=1) + 
  theme_bw()

ggplot(low_only_subset, aes(x=Sday_f, y=Uptake, group = Name)) + 
  geom_point(aes(color=Name, shape=Light), size= 4, position=position_dodge(width=0.4)) +
  scale_color_manual(values = c("#F5A70A", "#F5A70A", "#0A58F5", "#0A58F5")) +
  geom_errorbar(aes(ymin=Uptake-SE, ymax=Uptake+SE),color="gray30", width=0.2, size=0.7, position=position_dodge(width=0.4)) +
  theme_bw() +
  xlab("Sampling day") + ylab("Carbon uptake (%)") 

