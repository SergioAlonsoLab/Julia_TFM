#################
# New version ###
################# 

library(data.table)
library(ggplot2)
library(tidyr)

# set the working directory where the data is stored

setwd("~/sandbox/QPCR_Julia/")

# load the expression data

expression_data <- fread("expression.data.csv")

# expected columns
# Experiment:     Label for the experiment
# Sample:         Sample name (cell line, case, etc)
# Treatment:      Applied treatment
# Date:           Date of the experiment
# Target:         Target Gene
# meanCp_Target:  Mean Cp (Ct) of the Target
# errorCp_Target: Standard error of the mean Cp of the Target
# E_Target:       Efficiency of the Target PCR amplification
# Ref:            Reference Gene (in this version of the script, just one)
# meanCp_Ref:     Mean Cp (Ct) of the Reference
# errorCp_Ref:    Standard error of the mean Cp of the Reference
# E_Ref:          Efficiency of the Reference PCR amplification


# Indicate the efficiencies, if necessary

expression_data[Target=="ZNF549",E_Target:=2]
expression_data[Target=="STK33",E_Target:=2]
expression_data[Ref=="TPT1",E_Ref:=2]

# Calculate the relative expression in log2

expression_data[,log2_Rel := meanCp_Ref*log2(E_Ref) - meanCp_Target*log2(E_Target)]
expression_data[,log2_Rel_error := sqrt((errorCp_Ref*log2(E_Ref))^2 + (errorCp_Target*log2(E_Target))^2)]

# Indicate if the Cp are <35, >35 or >40

expression_data[meanCp_Target < 35,Cp:="Cp < 35"]
expression_data[meanCp_Target >= 35,Cp:="Cp > 35"]
expression_data[is.na(meanCp_Target),Cp :="N.D."]



# convert expression lineal scale

expression_data[,Rel := 2^log2_Rel]
expression_data[,Rel_max := 2^(log2_Rel+log2_Rel_error)]
expression_data[,Rel_min := 2^(log2_Rel-log2_Rel_error)]


# Place the samples with no amplification as a fixed log2 value (by default 1e-8)

expression_data[Cp =="N.D.", Rel := 1e-8]

# Create nice labels for the graphs

expression_data[,Label := paste("Experiment",as.numeric(factor(Experiment)))]

# Add methylation information 

methylation <- fread("methylation.csv")
methylation[,Sample:=gsub("( |-)","",Sample)]
expression_data <- merge(expression_data,methylation,by=c("Sample","Target","Treatment"),all.x=T,all.y=T)

# graphic themes

rotate_x_axis <- theme(axis.title.x = element_text(margin=margin(t=18)),
                       axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))
theme_light2 <-  theme_light() +
  theme(strip.background = element_rect(fill="white",color="grey"),
        strip.text.y = element_text(color="black",size=12),
        strip.text.x = element_text(color="black",size=12))


s_shape <- scale_shape_manual(values=c(19,21,4))



# Plot the methylation per Sample

ggplot(expression_data[,list(Sample,Target,Methylation)] %>% unique %>% na.omit) + 
  aes(Sample,Methylation) + 
  geom_hline(yintercept=c(.2,.5,.8),lty=2,lwd=0.2) +
  geom_col(aes(fill=Sample)) +
  facet_wrap(~Target) +
  theme_light2 +
  xlab("") +
  rotate_x_axis
  
# common geoms to plot expression in the y axis

common_geoms <- list(
  geom_hline(yintercept = 1e-8,linewidth = 6,color="grey90"),
  geom_errorbar(aes(ymin=Rel_min,ymax=Rel_max),width=0),
  geom_point(aes(shape=Cp),fill="white",size=2.5),
  theme_light2,
  ylab("Expression relative to TPT1 (log scale)"),
  scale_y_log10(),
  s_shape
)


# Plot the expression before treatment

ggplot(expression_data[Treatment=="DMSO" & !is.na(Cp)]) + 
  aes(Sample,Rel,color=Sample) +
  common_geoms +
  facet_wrap(~Target) +
  xlab("") +
  rotate_x_axis

# Plot methylation vs expression, before treatment

ggplot(expression_data[Treatment=="DMSO" & !is.na(Cp)]) + 
  aes(Methylation,Rel,color=Sample) + 
  geom_smooth(method="lm",color="black",fill="lightblue",lwd=.5) +
  common_geoms +
  facet_wrap(~Target) +
  xlab("Methylation") +
  rotate_x_axis

# Plot changes in expression after treatment

ggplot(expression_data[Sample!="WATER"]) + 
  aes(Label,Rel,color=Treatment) +
  common_geoms +
  facet_grid(cols=vars(Sample),rows=vars(Target)) +
  xlab("") +
  rotate_x_axis
  


