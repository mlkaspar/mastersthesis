rm(list=ls())
###
require(ggplot2)
require(scales)
require(plyr)
require(reshape2)
library(cowplot)
library("gridExtra")

# Initial --------------------------------------------------------------------
#load tables with morphotype counts
d3 = read.table(choose.files(), header = TRUE, sep=";")
d11 = read.table(choose.files(), header = TRUE, sep=";")
d18 = read.table(choose.files(), header = TRUE, sep=";")

# Calculate and plot morphotypes for each samplingday ------------------------

##############################d3####################################
#rm all bact
d3 = d3[d3$Type != "All_bact",]

######################################################strain 1, d3
d3_S1 = subset(d3, Strain == 1)
#aggregate by percentage
d3_S1_agg = aggregate(Type~Light, d3_S1, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_S1a = data.frame(Light=d3_S1_agg, d3_S1_agg$Type) 

#melt for long format
d3_S1am = melt(d3_S1a, id=1, variable.name="Type", value.name="Percentage")
d3_S1am = rename(d3_S1am, c("Light.Light" = "Light"))
d3_S1am = d3_S1am[-c(1,2),]
d3_S1am$label = paste0(sprintf("%.2f", (d3_S1am$Percentage*100), "%"))

#plot
d3_plot = ggplot(d3_S1am, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 3") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d3_plot
ggsave("d3_morphotypes.png", d3_plot, width = 5, height = 6)


######################################################
############################################################################################################strain 2, d3
d3_S2 = subset(d3, Strain == 2)
d3_S2_agg = aggregate(Type~Light, d3_S2, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_S2a = data.frame(Light=d3_S2_agg, d3_S2_agg$Type) 

d3_S2am = melt(d3_S2a, id=1, variable.name="Type", value.name="Percentage")
d3_S2am = rename(d3_S2am, c("Light.Light" = "Light"))
d3_S2am = d3_S2am[-c(1,2),]
d3_S2am$label = paste0(sprintf("%.2f", (d3_S2am$Percentage*100), "%"))


d3_plot2 = ggplot(d3_S2am, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 3") 
d3_plot2
ggsave("d3_morphotypes2.png", d3_plot2, width = 5, height = 6)


##############################d11####################################
# rm all bact
d11 = d11[d11$Type != "All_bact",]

######################################################strain 1, d11
d11_S1 = subset(d11, Strain == 1)
d11_S1_agg = aggregate(Type~Light, d11_S1, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_S1a = data.frame(Light=d11_S1_agg, d11_S1_agg$Type) 

d11_S1am = melt(d11_S1a, id=1, variable.name="Type", value.name="Percentage")
d11_S1am = rename(d11_S1am, c("Light.Light" = "Light"))
d11_S1am = d11_S1am[-c(1,2),]
d11_S1am$label = paste0(sprintf("%.2f", (d11_S1am$Percentage*100), "%"))


d11_plot = ggplot(d11_S1am, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 11")
d11_plot
ggsave("d11_morphotypes.png", d11_plot, width = 5, height = 6)

######################################################strain 2, d11
d11_S2 = subset(d11, Strain == 2)
d11_S2_agg = aggregate(Type~Light, d11_S2, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_S2a = data.frame(Light=d11_S2_agg, d11_S2_agg$Type) 

d11_S2am = melt(d11_S2a, id=1, variable.name="Type", value.name="Percentage")
d11_S2am = rename(d11_S2am, c("Light.Light" = "Light"))
d11_S2am = d11_S2am[-c(1,2),]
d11_S2am$label = paste0(sprintf("%.2f", (d11_S2am$Percentage*100), "%"))


d11_plot2 = ggplot(d11_S2am, aes(x=Light, y=Percentage, fill=Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 11")
d11_plot
ggsave("d11_morphotypes2.png", d11_plot2, width = 5, height = 6)

##############################d18####################################
#rm all bact
d18 = d18[d18$Type != "All_bact",]

######################################################strain 1, d18
d18_S1 = subset(d18, Strain == 1)
d18_S1_agg = aggregate(Type~Light, d18_S1, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_S1a = data.frame(Light=d18_S1_agg, d18_S1_agg$Type) 

d18_S1am = melt(d18_S1a, id=1, variable.name="Type", value.name="Percentage")
d18_S1am = rename(d18_S1am, c("Light.Light" = "Light"))
d18_S1am = d18_S1am[-c(1,2),]
d18_S1am$label = paste0(sprintf("%.2f", (d18_S1am$Percentage*100), "%"))


d18_plot = ggplot(d18_S1am, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 18")
d18_plot
ggsave("d18_morphotypes.png", d18_plot, width = 5, height = 6)

######################################################strain 2, d18
d18_S2 = subset(d18, Strain == 2)
d18_S2_agg = aggregate(Type~Light, d18_S2, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_S2a = data.frame(Light=d18_S2_agg, d18_S2_agg$Type) 

d18_S2am = melt(d18_S2a, id=1, variable.name="Type", value.name="Percentage")
d18_S2am = rename(d18_S2am, c("Light.Light" = "Light"))
d18_S2am = d18_S2am[-c(1,2),]
d18_S2am$label = paste0(sprintf("%.2f", (d18_S2am$Percentage*100), "%"))

d18_plot2 = ggplot(d18_S2am, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 18")
d18_plot2
ggsave("d18_morphotypes2.png", d18_plot2, width = 5, height = 6)

############################################################################################################
#Arrange plots in a grid

#save legend
legend = get_legend(d3_plot)

d3_str1_nolegend = d3_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d11_str1_nolegend = d11_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d18_str1_nolegend = d18_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))

d3_str2_nolegend = d3_plot2 + theme(legend.position = "none")+ ylab("") + theme(plot.title = element_text(size = 10))
d11_str2_nolegend = d11_plot2 + theme(legend.position = "none")+ ylab("") + theme(plot.title = element_text(size = 10))
d18_str2_nolegend = d18_plot2 + theme(legend.position = "none")+ ylab("") + theme(plot.title = element_text(size = 10))

grid.arrange(d3_str1_nolegend, d11_str1_nolegend, d18_str1_nolegend, legend, ncol=4, widths=c(1.5,1.5,1.5,0.8))
grid.arrange(d3_str2_nolegend, d11_str2_nolegend, d18_str2_nolegend, legend, ncol=4, widths=c(1.5,1.5,1.5,0.8))

grid.arrange(d3_str1_nolegend, d11_str1_nolegend, d18_str1_nolegend, legend, d3_str2_nolegend, d11_str2_nolegend, d18_str2_nolegend, legend, ncol=4, widths=c(1.5,1.5,1.5,0.8))


# Repeat for each subset of Low, Medium, High and for each samplingday -----------------------------------------------------------------

######################################################
######################################################
d3_rbind = rbind(d3_S1, d3_S2)
d3_rbind$Samplingday = c(rep("d3", times= (nrow(d3_rbind))))
d11_rbind = rbind(d11_S1, d11_S2)
d11_rbind$Samplingday = c(rep("d11", times= (nrow(d11_rbind))))
d18_rbind = rbind(d18_S1, d18_S2)
d18_rbind$Samplingday = c(rep("d18", times= (nrow(d18_rbind))))

morphoall = rbind(d3_rbind, d11_rbind, d18_rbind)
morpho_str1 = subset(morphoall, Strain == 1)
morpho_str1$type_f = factor(morpho_str1$Type, levels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci"))
ggplot(morpho_str1, aes(x=Samplingday)) + 
  geom_bar(stat="count", aes(fill=Type, y = (..count..)/sum(..count..))) + 
  facet_grid(Light ~ Medium) +
  scale_fill_brewer(palette="Set1", name="Morphotype") + #labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")
  ylab("Frequency in total dataset") + ggtitle("Morphotype frequencies, strain 1")

morpho_str2 = subset(morphoall, Strain == 2)
ggplot(morpho_str2, aes(x=Samplingday)) + 
  geom_bar(stat="count", aes(fill=Type, y = (..count..)/sum(..count..))) + 
  facet_grid(Light ~ Medium) +
  scale_fill_brewer(palette="Set1", name="Morphotype") + #labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")
  ylab("Frequency in total dataset") + ggtitle("Morphotype frequencies, strain 2")

######################################################
######################################################
#medium subsets
#d3
#low
d3_s1_low = subset(d3_S1, Medium == "L")
d3_s1_low_agg = aggregate(Type~Light, d3_s1_low, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_s1_lowa = data.frame(Light=d3_s1_low_agg, d3_s1_low_agg$Type) 

d3_s1_lowam = melt(d3_s1_lowa, id=1, variable.name="Type", value.name="Percentage")
d3_s1_lowam = rename(d3_s1_lowam, c("Light.Light" = "Light"))
d3_s1_lowam = d3_s1_lowam[-c(1,2),]
d3_s1_lowam$label = paste0(sprintf("%.2f", (d3_s1_lowam$Percentage*100), "%"))


d3_S1_low_plot = ggplot(d3_s1_lowam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 3") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d3_S1_low_plot
#medium

d3_s1_med = subset(d3_S1, Medium == "M")
d3_s1_med_agg = aggregate(Type~Light, d3_s1_med, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_s1_meda = data.frame(Light=d3_s1_med_agg, d3_s1_med_agg$Type) 

d3_s1_medam = melt(d3_s1_meda, id=1, variable.name="Type", value.name="Percentage")
d3_s1_medam = rename(d3_s1_medam, c("Light.Light" = "Light"))
d3_s1_medam = d3_s1_medam[-c(1,2),]
d3_s1_medam$label = paste0(sprintf("%.2f", (d3_s1_medam$Percentage*100), "%"))


d3_S1_med_plot = ggplot(d3_s1_medam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 3") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d3_S1_med_plot

#high
d3_s1_high = subset(d3_S1, Medium == "H")
d3_s1_high_agg = aggregate(Type~Light, d3_s1_high, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_s1_higha = data.frame(Light=d3_s1_high_agg, d3_s1_high_agg$Type) 

d3_s1_higham = melt(d3_s1_higha, id=1, variable.name="Type", value.name="Percentage")
d3_s1_higham = rename(d3_s1_higham, c("Light.Light" = "Light"))
d3_s1_higham = d3_s1_higham[-c(1,2),]
d3_s1_higham$label = paste0(sprintf("%.2f", (d3_s1_higham$Percentage*100), "%"))


d3_S1_high_plot = ggplot(d3_s1_higham, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 3") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d3_S1_high_plot
###############################
d3_s2_low = subset(d3_S2, Medium == "L")
d3_s2_low_agg = aggregate(Type~Light, d3_s2_low, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_s2_lowa = data.frame(Light=d3_s2_low_agg, d3_s2_low_agg$Type) 

d3_s2_lowam = melt(d3_s2_lowa, id=1, variable.name="Type", value.name="Percentage")
d3_s2_lowam = rename(d3_s2_lowam, c("Light.Light" = "Light"))
d3_s2_lowam = d3_s2_lowam[-c(1,2),]
d3_s2_lowam$label = paste0(sprintf("%.2f", (d3_s2_lowam$Percentage*100), "%"))


d3_s2_low_plot = ggplot(d3_s2_lowam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 3") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d3_s2_low_plot
#medium

d3_s2_med = subset(d3_S2, Medium == "M")
d3_s2_med_agg = aggregate(Type~Light, d3_s2_med, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_s2_meda = data.frame(Light=d3_s2_med_agg, d3_s2_med_agg$Type) 

d3_s2_medam = melt(d3_s2_meda, id=1, variable.name="Type", value.name="Percentage")
d3_s2_medam = rename(d3_s2_medam, c("Light.Light" = "Light"))
d3_s2_medam = d3_s2_medam[-c(1,2),]
d3_s2_medam$label = paste0(sprintf("%.2f", (d3_s2_medam$Percentage*100), "%"))


d3_s2_med_plot = ggplot(d3_s2_medam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 3") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d3_s2_med_plot

#high
d3_s2_high = subset(d3_S2, Medium == "H")
d3_s2_high_agg = aggregate(Type~Light, d3_s2_high, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d3_s2_higha = data.frame(Light=d3_s2_high_agg, d3_s2_high_agg$Type) 

d3_s2_higham = melt(d3_s2_higha, id=1, variable.name="Type", value.name="Percentage")
d3_s2_higham = rename(d3_s2_higham, c("Light.Light" = "Light"))
d3_s2_higham = d3_s2_higham[-c(1,2),]
d3_s2_higham$label = paste0(sprintf("%.2f", (d3_s2_higham$Percentage*100), "%"))


d3_s2_high_plot = ggplot(d3_s2_higham, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 3") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d3_s2_high_plot

#####################################################
#####################################################
#d11
#low
d11_s1_low = subset(d11_S1, Medium == "L")
d11_s1_low_agg = aggregate(Type~Light, d11_s1_low, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_s1_lowa = data.frame(Light=d11_s1_low_agg, d11_s1_low_agg$Type) 

d11_s1_lowam = melt(d11_s1_lowa, id=1, variable.name="Type", value.name="Percentage")
d11_s1_lowam = rename(d11_s1_lowam, c("Light.Light" = "Light"))
d11_s1_lowam = d11_s1_lowam[-c(1,2),]
d11_s1_lowam$label = paste0(sprintf("%.2f", (d11_s1_lowam$Percentage*100), "%"))


d11_S1_low_plot = ggplot(d11_s1_lowam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 11") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d11_S1_low_plot
#medium

d11_s1_med = subset(d11_S1, Medium == "M")
d11_s1_med_agg = aggregate(Type~Light, d11_s1_med, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_s1_meda = data.frame(Light=d11_s1_med_agg, d11_s1_med_agg$Type) 

d11_s1_medam = melt(d11_s1_meda, id=1, variable.name="Type", value.name="Percentage")
d11_s1_medam = rename(d11_s1_medam, c("Light.Light" = "Light"))
d11_s1_medam = d11_s1_medam[-c(1,2),]
d11_s1_medam$label = paste0(sprintf("%.2f", (d11_s1_medam$Percentage*100), "%"))


d11_S1_med_plot = ggplot(d11_s1_medam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 11") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d11_S1_med_plot

#high
d11_s1_high = subset(d11_S1, Medium == "H")
d11_s1_high_agg = aggregate(Type~Light, d11_s1_high, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_s1_higha = data.frame(Light=d11_s1_high_agg, d11_s1_high_agg$Type) 

d11_s1_higham = melt(d11_s1_higha, id=1, variable.name="Type", value.name="Percentage")
d11_s1_higham = rename(d11_s1_higham, c("Light.Light" = "Light"))
d11_s1_higham = d11_s1_higham[-c(1,2),]
d11_s1_higham$label = paste0(sprintf("%.2f", (d11_s1_higham$Percentage*100), "%"))


d11_S1_high_plot = ggplot(d11_s1_higham, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 11") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d11_S1_high_plot
###############################
d11_s2_low = subset(d11_S2, Medium == "L")
d11_s2_low_agg = aggregate(Type~Light, d11_s2_low, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_s2_lowa = data.frame(Light=d11_s2_low_agg, d11_s2_low_agg$Type) 

d11_s2_lowam = melt(d11_s2_lowa, id=1, variable.name="Type", value.name="Percentage")
d11_s2_lowam = rename(d11_s2_lowam, c("Light.Light" = "Light"))
d11_s2_lowam = d11_s2_lowam[-c(1,2),]
d11_s2_lowam$label = paste0(sprintf("%.2f", (d11_s2_lowam$Percentage*100), "%"))


d11_s2_low_plot = ggplot(d11_s2_lowam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 11") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d11_s2_low_plot
#medium

d11_s2_med = subset(d11_S2, Medium == "M")
d11_s2_med_agg = aggregate(Type~Light, d11_s2_med, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_s2_meda = data.frame(Light=d11_s2_med_agg, d11_s2_med_agg$Type) 

d11_s2_medam = melt(d11_s2_meda, id=1, variable.name="Type", value.name="Percentage")
d11_s2_medam = rename(d11_s2_medam, c("Light.Light" = "Light"))
d11_s2_medam = d11_s2_medam[-c(1,2),]
d11_s2_medam$label = paste0(sprintf("%.2f", (d11_s2_medam$Percentage*100), "%"))


d11_s2_med_plot = ggplot(d11_s2_medam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 11") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d11_s2_med_plot

#high
d11_s2_high = subset(d11_S2, Medium == "H")
d11_s2_high_agg = aggregate(Type~Light, d11_s2_high, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d11_s2_higha = data.frame(Light=d11_s2_high_agg, d11_s2_high_agg$Type) 

d11_s2_higham = melt(d11_s2_higha, id=1, variable.name="Type", value.name="Percentage")
d11_s2_higham = rename(d11_s2_higham, c("Light.Light" = "Light"))
d11_s2_higham = d11_s2_higham[-c(1,2),]
d11_s2_higham$label = paste0(sprintf("%.2f", (d11_s2_higham$Percentage*100), "%"))


d11_s2_high_plot = ggplot(d11_s2_higham, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 11") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d11_s2_high_plot

#####################################################
#####################################################
#d18
#low
d18_s1_low = subset(d18_S1, Medium == "L")
d18_s1_low_agg = aggregate(Type~Light, d18_s1_low, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_s1_lowa = data.frame(Light=d18_s1_low_agg, d18_s1_low_agg$Type) 

d18_s1_lowam = melt(d18_s1_lowa, id=1, variable.name="Type", value.name="Percentage")
d18_s1_lowam = rename(d18_s1_lowam, c("Light.Light" = "Light"))
d18_s1_lowam = d18_s1_lowam[-c(1,2),]
d18_s1_lowam$label = paste0(sprintf("%.2f", (d18_s1_lowam$Percentage*100), "%"))


d18_S1_low_plot = ggplot(d18_s1_lowam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 18") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d18_S1_low_plot
#medium

d18_s1_med = subset(d18_S1, Medium == "M")
d18_s1_med_agg = aggregate(Type~Light, d18_s1_med, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_s1_meda = data.frame(Light=d18_s1_med_agg, d18_s1_med_agg$Type) 

d18_s1_medam = melt(d18_s1_meda, id=1, variable.name="Type", value.name="Percentage")
d18_s1_medam = rename(d18_s1_medam, c("Light.Light" = "Light"))
d18_s1_medam = d18_s1_medam[-c(1,2),]
d18_s1_medam$label = paste0(sprintf("%.2f", (d18_s1_medam$Percentage*100), "%"))


d18_S1_med_plot = ggplot(d18_s1_medam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 18") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d18_S1_med_plot

#high
d18_s1_high = subset(d18_S1, Medium == "H")
d18_s1_high_agg = aggregate(Type~Light, d18_s1_high, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_s1_higha = data.frame(Light=d18_s1_high_agg, d18_s1_high_agg$Type) 

d18_s1_higham = melt(d18_s1_higha, id=1, variable.name="Type", value.name="Percentage")
d18_s1_higham = rename(d18_s1_higham, c("Light.Light" = "Light"))
d18_s1_higham = d18_s1_higham[-c(1,2),]
d18_s1_higham$label = paste0(sprintf("%.2f", (d18_s1_higham$Percentage*100), "%"))


d18_S1_high_plot = ggplot(d18_s1_higham, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 1, sampling day 18") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d18_S1_high_plot
###############################
d18_s2_low = subset(d18_S2, Medium == "L")
d18_s2_low_agg = aggregate(Type~Light, d18_s2_low, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_s2_lowa = data.frame(Light=d18_s2_low_agg, d18_s2_low_agg$Type) 

d18_s2_lowam = melt(d18_s2_lowa, id=1, variable.name="Type", value.name="Percentage")
d18_s2_lowam = rename(d18_s2_lowam, c("Light.Light" = "Light"))
d18_s2_lowam = d18_s2_lowam[-c(1,2),]
d18_s2_lowam$label = paste0(sprintf("%.2f", (d18_s2_lowam$Percentage*100), "%"))


d18_s2_low_plot = ggplot(d18_s2_lowam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 18") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d18_s2_low_plot
#medium

d18_s2_med = subset(d18_S2, Medium == "M")
d18_s2_med_agg = aggregate(Type~Light, d18_s2_med, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_s2_meda = data.frame(Light=d18_s2_med_agg, d18_s2_med_agg$Type) 

d18_s2_medam = melt(d18_s2_meda, id=1, variable.name="Type", value.name="Percentage")
d18_s2_medam = rename(d18_s2_medam, c("Light.Light" = "Light"))
d18_s2_medam = d18_s2_medam[-c(1,2),]
d18_s2_medam$label = paste0(sprintf("%.2f", (d18_s2_medam$Percentage*100), "%"))


d18_s2_med_plot = ggplot(d18_s2_medam, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 18") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d18_s2_med_plot

#high
d18_s2_high = subset(d18_S2, Medium == "H")
d18_s2_high_agg = aggregate(Type~Light, d18_s2_high, function(x)c(Large_Rods=sum(x=="Rhods")/length(x),Vibrio=sum(x=="Vibrio")/length(x),Large_cocci=sum(x=="Large_cocci")/length(x), Small_rods = sum(x=="Small_rhods")/length(x), Small_cocci=sum(x=="Small_cocci")/length(x)))
d18_s2_higha = data.frame(Light=d18_s2_high_agg, d18_s2_high_agg$Type) 

d18_s2_higham = melt(d18_s2_higha, id=1, variable.name="Type", value.name="Percentage")
d18_s2_higham = rename(d18_s2_higham, c("Light.Light" = "Light"))
d18_s2_higham = d18_s2_higham[-c(1,2),]
d18_s2_higham$label = paste0(sprintf("%.2f", (d18_s2_higham$Percentage*100), "%"))


d18_s2_high_plot = ggplot(d18_s2_higham, aes(x=Light, y=Percentage, fill = Type)) + 
  geom_bar(stat="identity", width=0.5, position=position_stack()) +
  scale_fill_brewer(palette="Set1", name="Morphotype", labels=c("Large rods", "Vibrio-shaped", "Large cocci", "Small rods", "Small cocci")) + 
  scale_y_continuous(labels=percent) + 
  scale_x_discrete(labels=c("Dark", "Light")) + xlab("") + 
  theme_classic() +
  ggtitle("Strain 2, sampling day 18") +
  geom_text(aes(label=label), size = 3, position=position_stack(vjust = 0.5), colour = "white") 
d18_s2_high_plot
########################################################
########################################################
########################################################
#Arrange plots for strain 1
#save legend
legend = get_legend(d3_S1_low_plot)

d3_str1L_nolegend = d3_S1_low_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d11_str1L_nolegend = d11_S1_low_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d18_str1L_nolegend = d18_S1_low_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))

d3_str1M_nolegend = d3_S1_med_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d11_str1M_nolegend = d11_S1_med_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d18_str1M_nolegend = d18_S1_med_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))

d3_str1H_nolegend = d3_S1_high_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d11_str1H_nolegend = d11_S1_high_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d18_str1H_nolegend = d18_S1_high_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))

grid.arrange(d3_str1L_nolegend,d11_str1L_nolegend, d18_str1L_nolegend, legend, d3_str1M_nolegend, d11_str1M_nolegend, d18_str1M_nolegend,legend, d3_str1H_nolegend,d11_str1H_nolegend, d18_str1H_nolegend, legend, ncol=4,widths=c(1.5,1.5,1.5,0.8))

########################################################
#Arrange plots for strain 2
legend = get_legend(d3_s2_low_plot)

d3_str1L_nolegend = d3_s2_low_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d11_str1L_nolegend = d11_s2_low_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d18_str1L_nolegend = d18_s2_low_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))

d3_str1M_nolegend = d3_s2_med_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d11_str1M_nolegend = d11_s2_med_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d18_str1M_nolegend = d18_s2_med_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))

d3_str1H_nolegend = d3_s2_high_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d11_str1H_nolegend = d11_s2_high_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))
d18_str1H_nolegend = d18_s2_high_plot + theme(legend.position = "none") + ylab("") + theme(plot.title = element_text(size = 10))

grid.arrange(d3_str1L_nolegend,d11_str1L_nolegend, d18_str1L_nolegend, legend, d3_str1M_nolegend, d11_str1M_nolegend, d18_str1M_nolegend,legend, d3_str1H_nolegend,d11_str1H_nolegend, d18_str1H_nolegend, legend, ncol=4,widths=c(1.5,1.5,1.5,0.8))

