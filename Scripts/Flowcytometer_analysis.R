rm(list=ls())
##############################################################
###Required packages
require(ggplot2)
require(scales)
require(reshape2)
##############################################################
###Load table with data from the Flowycytometer
FlowCytoTable = read.table(choose.files(), header=TRUE, sep=";")
##############################################################

# Initial -----------------------------------------------------------------

#Keep order of samples:
FlowCytoTable$Name = factor(FlowCytoTable$Name, levels=c("DL1", "LL1", "DM1", "LM1", "DH1", "LH1", "DL2", "LL2", "DM2", "LM2", "DH2", "LH2"))
#Calculate limits for errorbars
limits = aes(ymax = Mean + SE, ymin = Mean - SE, colour=FlowCytoTable$Name)
#function for y-scale format:
scientific_10 <- function(x) {
  parse(text=gsub("1e", "10^", scientific_format()(x)))
}
#Define order of colors for graph:
reds = c("darkorange","darkorange", "#FF3333", "#FF3333", "#9F3620", "#9F3620")
blues = c("cadetblue2", "cadetblue2", "#33BBFF", "#33BBFF", "#000FFF", "#000FFF")
colrs = append(reds, blues)
#black and white:
bw = c("#E7E7E7", "#E7E7E7", "#CACACA", "#CACACA", "#909090", "#909090", "#737373", "#737373", "#575757", "#575757", "#1D1D1D", "#1D1D1D")
##############################################################
# Graphs -----------------------------------------------------------------
###basic plot:
basic_plot = ggplot(data=FlowCytoTable, aes(x=TotDays, y=Mean, group=Name, linetype=Light)) + 
  geom_line(aes(colour=Name), size=1) + 
  geom_errorbar(limits, width=0.03, linetype=1) + 
  ggtitle("Growth curves") + 
  xlab("Time (days)") +
  scale_x_continuous(breaks=seq(0,110,2))+
  #xlab("Time (h)") + 
  ylab("Number of cells/ml") +
  theme_bw()
#show plot:
basic_plot

###log-scale, color, no shapes:
log_col = basic_plot +
  scale_y_log10(limits=c(1e03, 1e08), breaks=c(1e04, 1e05, 1e06, 1e07, 1e08), labels = scientific_10) +  
  scale_color_manual(values = colrs)  +
  annotation_logticks(sides="l") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=60, vjust=0.5, size=7))
#show plot:
log_col

###log-scale, black and white, shapes:
log_bw = basic_plot +
  geom_point(aes(shape=ID), size=2) +
  scale_y_log10(limits=c(1e02, 1e08), breaks=c(1e04, 1e05, 1e06, 1e07, 1e08)) +  
  scale_color_manual(values = bw)  +
  annotation_logticks(sides="l") +
  scale_shape(guide=FALSE) +
  guides(shape=FALSE)
#show plot:
log_bw

##facets
#keep order of ID variable
FlowCytoTable$IDf = factor(FlowCytoTable$ID, levels=c("L1", "M1", "H1", "L2", "M2", "H2"))

fac_col_log = ggplot(data=FlowCytoTable, aes(x=TotDays, y=Mean, group=Name, linetype=Light)) + 
  geom_line(aes(colour=Name), size=1) + 
  geom_errorbar(limits, width=0.03, linetype=1) + 
  ggtitle("Growth curves") + 
  xlab("Time (days)") +
  scale_x_continuous(breaks=seq(0,110,5)) +
  ylab("Number of cells/ml") +
  scale_color_manual(values = colrs)  +
  scale_y_log10(limits=c(1e03, 1e08), breaks=c(1e04, 1e05, 1e06, 1e07, 1e08), labels = scientific_10) +   #scale_color_manual(values = colrs)
  annotation_logticks(sides="l") +
  theme_classic() +
  #facet_grid(. ~ IDf) #vertical ~horizontal
  facet_wrap( ~ IDf, ncol=3) +
  theme(axis.text.x = element_text(angle=60, vjust=0.5, size=7))
#show plot:
fac_col_log


# Statistics -----------------------------------------------------------------
sink("growthexperiment_stats.txt")

##########Strain 1, all sampling days
strain1sub = subset(FlowCytoTable, Strain == "Rhod")
strain1sub$Logs = log(strain1sub$Mean)
#####
print("Strain1")
Low1_subset = subset(strain1sub, ID == "L1")
Med1_subset = subset(strain1sub, ID == "M1")
High1_subset = subset(strain1sub, ID == "H1")
#####
print("Strain 1, low")
Low1_subset.lm = lm(data=Low1_subset, Logs ~ Light)
summary(Low1_subset.lm)
t.test(data=Low1_subset, Logs ~ Light)
#####
print("Strain 1, medium")
Med1_subset.lm = lm(data=Med1_subset, Logs ~ Light)
summary(Med1_subset.lm)
t.test(data=Med1_subset, Logs ~ Light)
#####
print("Strain 1, high")
High1_subset.lm = lm(data=High1_subset, Logs ~ Light)
summary(High1_subset.lm)
t.test(data=High1_subset, Logs ~ Light)
##########
print("Strain 2")
strain2sub = subset(FlowCytoTable, Strain == "NoRhod")
strain2sub$Logs = log(strain2sub$Mean)
#####
Low2_subset = subset(strain2sub, ID == "L2")
Med2_subset = subset(strain2sub, ID == "M2")
High2_subset = subset(strain2sub, ID == "H2")
#####
print("Strain 2, low")
Low2_subset.lm = lm(data=Low2_subset, Logs ~ Light)
summary(Low2_subset.lm)
t.test(data=Low2_subset, Logs ~ Light)
#####
print("Strain 2, medium")
Med2_subset.lm = lm(data=Med2_subset, Logs ~ Light)
summary(Med2_subset.lm)
t.test(data=Med2_subset, Logs ~ Light)
#####
print("Strain 2, high")
High2_subset.lm = lm(data=High2_subset, Logs ~ Light)
summary(High2_subset.lm)
t.test(data=High2_subset, Logs ~ Light)
##########
sink()

##################################################################
##################################################################
## with all triplicates instead of mean
#load table with data of all triplicates from the Flowcytometer
Alltripl = read.table(choose.files(), header=TRUE, sep=";")

#subsampling
AlltriplAll = Alltripl
Alltripld22 = subset(Alltripl, Samplingday == "d22")
Alltripld18 = subset(Alltripl, Samplingday == "d18")
Alltripld11 = subset(Alltripl, Samplingday == "d11")
Alltripld3 = subset(Alltripl, Samplingday == "d3")

################All days
sink("growthexperiment_stats_All_log.txt")
##########
#####
print("Strain1")
Low1_subset = subset(AlltriplAll, ID == "L1")
Med1_subset = subset(AlltriplAll, ID == "M1")
High1_subset = subset(AlltriplAll, ID == "H1")
#####
print("Strain 1, low")
Low1_subset.lm = lm(data=Low1_subset, Logs ~ Samplingday + Light)
summary(Low1_subset.lm)
anova(Low1_subset.lm)
t.test(data=Low1_subset, Logs ~ Light)
a1 = aov(data=Low1_subset, Logs ~ Light + Samplingday)
TukeyHSD(x=a1, "Samplingday", conf.level=0.95)

#####
print("Strain 1, medium")
Med1_subset.lm = lm(data=Med1_subset, Logs ~ Light)
summary(Med1_subset.lm)
t.test(data=Med1_subset, Logs ~ Light)
#####
print("Strain 1, high")
High1_subset.lm = lm(data=High1_subset, Logs ~ Light)
summary(High1_subset.lm)
t.test(data=High1_subset, Logs ~ Light)
##########
print("Strain 2")
#####
Low2_subset = subset(AlltriplAll, ID == "L2")
Med2_subset = subset(AlltriplAll, ID == "M2")
High2_subset = subset(AlltriplAll, ID == "H2")
#####
print("Strain 2, low")
Low2_subset.lm = lm(data=Low2_subset, Logs ~ Light)
summary(Low2_subset.lm)
t.test(data=Low2_subset, Logs ~ Light)
#####
print("Strain 2, medium")
Med2_subset.lm = lm(data=Med2_subset, Logs ~ Light)
summary(Med2_subset.lm)
t.test(data=Med2_subset, Logs ~ Light)
#####
print("Strain 2, high")
High2_subset.lm = lm(data=High2_subset, Logs ~ Light)
summary(High2_subset.lm)
t.test(data=High2_subset, Logs ~ Light)
##########
sink()

#######################samplingday 3

sink("growthexperiment_stats_d3_log.txt")
##########
#####
print("Strain1")
Low1_subset = subset(Alltripld3, ID == "L1")
Med1_subset = subset(Alltripld3, ID == "M1")
High1_subset = subset(Alltripld3, ID == "H1")
#####
print("Strain 1, low")
Low1_subset.lm = lm(data=Low1_subset, Logs ~ Samplingday + Light)
summary(Low1_subset.lm)
anova(Low1_subset.lm)
t.test(data=Low1_subset, Logs ~ Light)
a1 = aov(data=Low1_subset, Logs ~ Light + Samplingday)
TukeyHSD(x=a1, "Samplingday", conf.level=0.95)

#####
print("Strain 1, medium")
Med1_subset.lm = lm(data=Med1_subset, Logs ~ Light)
summary(Med1_subset.lm)
t.test(data=Med1_subset, Logs ~ Light)
#####
print("Strain 1, high")
High1_subset.lm = lm(data=High1_subset, Logs ~ Light)
summary(High1_subset.lm)
t.test(data=High1_subset, Logs ~ Light)
##########
print("Strain 2")
#####
Low2_subset = subset(Alltripld3, ID == "L2")
Med2_subset = subset(Alltripld3, ID == "M2")
High2_subset = subset(Alltripld3, ID == "H2")
#####
print("Strain 2, low")
Low2_subset.lm = lm(data=Low2_subset, Logs ~ Light)
summary(Low2_subset.lm)
t.test(data=Low2_subset, Logs ~ Light)
#####
print("Strain 2, medium")
Med2_subset.lm = lm(data=Med2_subset, Logs ~ Light)
summary(Med2_subset.lm)
t.test(data=Med2_subset, Logs ~ Light)
#####
print("Strain 2, high")
High2_subset.lm = lm(data=High2_subset, Logs ~ Light)
summary(High2_subset.lm)
t.test(data=High2_subset, Logs ~ Light)
##########
sink()

#######################samplingday 11

sink("growthexperiment_stats_d11_log.txt")
##########
#####
print("Strain1")
Low1_subset = subset(Alltripld11, ID == "L1")
Med1_subset = subset(Alltripld11, ID == "M1")
High1_subset = subset(Alltripld11, ID == "H1")
#####
print("Strain 1, low")
Low1_subset.lm = lm(data=Low1_subset, Logs ~ Samplingday + Light)
summary(Low1_subset.lm)
anova(Low1_subset.lm)
t.test(data=Low1_subset, Logs ~ Light)
a1 = aov(data=Low1_subset, Logs ~ Light + Samplingday)
TukeyHSD(x=a1, "Samplingday", conf.level=0.95)

#####
print("Strain 1, medium")
Med1_subset.lm = lm(data=Med1_subset, Logs ~ Light)
summary(Med1_subset.lm)
t.test(data=Med1_subset, Logs ~ Light)
#####
print("Strain 1, high")
High1_subset.lm = lm(data=High1_subset, Logs ~ Light)
summary(High1_subset.lm)
t.test(data=High1_subset, Logs ~ Light)
##########
print("Strain 2")
#####
Low2_subset = subset(Alltripld11, ID == "L2")
Med2_subset = subset(Alltripld11, ID == "M2")
High2_subset = subset(Alltripld11, ID == "H2")
#####
print("Strain 2, low")
Low2_subset.lm = lm(data=Low2_subset, Logs ~ Light)
summary(Low2_subset.lm)
t.test(data=Low2_subset, Logs ~ Light)
#####
print("Strain 2, medium")
Med2_subset.lm = lm(data=Med2_subset, Logs ~ Light)
summary(Med2_subset.lm)
t.test(data=Med2_subset, Logs ~ Light)
#####
print("Strain 2, high")
High2_subset.lm = lm(data=High2_subset, Logs ~ Light)
summary(High2_subset.lm)
t.test(data=High2_subset, Logs ~ Light)
##########
sink()

#######################samplingday 18

sink("growthexperiment_stats_d18_log.txt")
##########
#####
print("Strain1")
Low1_subset = subset(Alltripld18, ID == "L1")
Med1_subset = subset(Alltripld18, ID == "M1")
High1_subset = subset(Alltripld18, ID == "H1")
#####
print("Strain 1, low")
Low1_subset.lm = lm(data=Low1_subset, Logs ~ Samplingday + Light)
summary(Low1_subset.lm)
anova(Low1_subset.lm)
t.test(data=Low1_subset, Logs ~ Light)
a1 = aov(data=Low1_subset, Logs ~ Light + Samplingday)
TukeyHSD(x=a1, "Samplingday", conf.level=0.95)

#####
print("Strain 1, medium")
Med1_subset.lm = lm(data=Med1_subset, Logs ~ Light)
summary(Med1_subset.lm)
t.test(data=Med1_subset, Logs ~ Light)
#####
print("Strain 1, high")
High1_subset.lm = lm(data=High1_subset, Logs ~ Light)
summary(High1_subset.lm)
t.test(data=High1_subset, Logs ~ Light)
##########
print("Strain 2")
#####
Low2_subset = subset(Alltripld18, ID == "L2")
Med2_subset = subset(Alltripld18, ID == "M2")
High2_subset = subset(Alltripld18, ID == "H2")
#####
print("Strain 2, low")
Low2_subset.lm = lm(data=Low2_subset, Logs ~ Light)
summary(Low2_subset.lm)
t.test(data=Low2_subset, Logs ~ Light)
#####
print("Strain 2, medium")
Med2_subset.lm = lm(data=Med2_subset, Logs ~ Light)
summary(Med2_subset.lm)
t.test(data=Med2_subset, Logs ~ Light)
#####
print("Strain 2, high")
High2_subset.lm = lm(data=High2_subset, Logs ~ Light)
summary(High2_subset.lm)
t.test(data=High2_subset, Logs ~ Light)
##########
sink()

#######################samplingday 22

sink("growthexperiment_stats_d22_log.txt")
##########
#####
print("Strain1")
Low1_subset = subset(Alltripld22, ID == "L1")
Med1_subset = subset(Alltripld22, ID == "M1")
High1_subset = subset(Alltripld22, ID == "H1")
#####
print("Strain 1, low")
Low1_subset.lm = lm(data=Low1_subset, Logs ~ Samplingday + Light)
summary(Low1_subset.lm)
anova(Low1_subset.lm)
t.test(data=Low1_subset, Logs ~ Light)
a1 = aov(data=Low1_subset, Logs ~ Light + Samplingday)
TukeyHSD(x=a1, "Samplingday", conf.level=0.95)

#####
print("Strain 1, medium")
Med1_subset.lm = lm(data=Med1_subset, Logs ~ Light)
summary(Med1_subset.lm)
t.test(data=Med1_subset, Logs ~ Light)
#####
print("Strain 1, high")
High1_subset.lm = lm(data=High1_subset, Logs ~ Light)
summary(High1_subset.lm)
t.test(data=High1_subset, Logs ~ Light)
##########
print("Strain 2")
#####
Low2_subset = subset(Alltripld22, ID == "L2")
Med2_subset = subset(Alltripld22, ID == "M2")
High2_subset = subset(Alltripld22, ID == "H2")
#####
print("Strain 2, low")
Low2_subset.lm = lm(data=Low2_subset, Logs ~ Light)
summary(Low2_subset.lm)
t.test(data=Low2_subset, Logs ~ Light)
#####
print("Strain 2, medium")
Med2_subset.lm = lm(data=Med2_subset, Logs ~ Light)
summary(Med2_subset.lm)
t.test(data=Med2_subset, Logs ~ Light)
#####
print("Strain 2, high")
High2_subset.lm = lm(data=High2_subset, Logs ~ Light)
summary(High2_subset.lm)
t.test(data=High2_subset, Logs ~ Light)
##########
sink()

##################################################################
##################################################################
##################################################################
#ANOVA of all subsets

models = sapply(x, function(my){
  lm(Logs ~ Light, data=y, Samplingday == my)
}, simplify=FALSE)

x = unique(Low1_subset$Samplingday)
y = Low1_subset
models = sapply(x, function(my){
  lm(Logs ~ Light, data=y, Samplingday == my)
}, simplify=FALSE)
ANOVA.tables_low1 = sapply(models, anova, simplify=FALSE)
sink("anova_low1.txt")
ANOVA.tables_low1
sink()

x = unique(Med1_subset$Samplingday)
y = Med1_subset
models = sapply(x, function(my){
  lm(Logs ~ Light, data=y, Samplingday == my)
}, simplify=FALSE)
ANOVA.tables_med1 = sapply(models, anova, simplify=FALSE)
sink("anova_med1.txt")
ANOVA.tables_med1
sink()

x = unique(High1_subset$Samplingday)
y = High1_subset
models = sapply(x, function(my){
  lm(Logs ~ Light, data=y, Samplingday == my)
}, simplify=FALSE)
ANOVA.tables_high1 = sapply(models, anova, simplify=FALSE)
sink("anova_high1.txt")
ANOVA.tables_high1
sink()

x = unique(Low2_subset$Samplingday)
y = Low2_subset
models = sapply(x, function(my){
  lm(Logs ~ Light, data=y, Samplingday == my)
}, simplify=FALSE)
ANOVA.tables_low2 = sapply(models, anova, simplify=FALSE)
sink("anova_low2.txt")
ANOVA.tables_low2
sink()

x = unique(Med2_subset$Samplingday)
y = Med2_subset
models = sapply(x, function(my){
  lm(Logs ~ Light, data=y, Samplingday == my)
}, simplify=FALSE)
ANOVA.tables_med2 = sapply(models, anova, simplify=FALSE)
sink("anova_med2.txt")
ANOVA.tables_med2
sink()

x = unique(High2_subset$Samplingday)
y = High2_subset
models = sapply(x, function(my){
  lm(Logs ~ Light, data=y, Samplingday == my)
}, simplify=FALSE)
ANOVA.tables_high2 = sapply(models, anova, simplify=FALSE)
sink("anova_high2.txt")
ANOVA.tables_high2
sink()

