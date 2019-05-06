library(ggplot2)
library(reshape2)
library(dplyr)

data <- read.csv("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#3/Genomic Prediction/all_accuracy_12traits.csv")
colnames(data) <- c("Trait","Day","90","80","70","60","50","40","30","20","10")
new_data<-melt(data=data,id.vars=c("Trait","Day"),variable.name="Ratio",value.name="Accuracy")
head(new_data)
new_data$Ratio <- as.numeric(as.character(new_data$Ratio))
new_data$Day <- as.factor(new_data$Day)
str(new_data)

summary(new_data$Accuracy)
# Faceting

tiff("C:/Users/falk/Google Drive/PhD/PhD Projects/Blue Steel/Paper#3/Genomic Prediction/Prediction_Accuracy_March24.tiff", units = "in" ,height = 7, width = 9, res=200)
ggplot(new_data, mapping = aes(y=Accuracy, x=Ratio, color=Day, group = Day, fill=Day)) + 
  geom_line(size = 1) +
  geom_point(size = 0.5) +
  facet_wrap(~Trait, ncol = 4, labeller = label_parsed) +
  xlab("Training Percent age") +
  ylab("Prediction Accuracy") +
  scale_x_continuous(breaks=c(20,40,60,80))+
  theme_bw() +
  theme(panel.background = element_blank(), strip.background = element_blank())
dev.off()
