library(tidyverse) #loads program for later use

#########################################################################

COHA.prey = read.csv("COHA_preyweights (3).csv") #pulls in COHA Prey Weights
COHA.prey = na.omit(COHA.prey) #omits rows with nothing in them
View(COHA.prey) #just lets us see the data
COHA.prey=COHA.prey[!grepl("na", COHA.prey$Taxa),] #removes the cells with "na" written in them in the "Taxa" column only

#########################################################################

md= COHA.prey[c("freq", "Taxa")] #creates a new data.frame called "md" with only the columns needed for graph
md = aggregate(.~Taxa,data=md,FUN=sum) #adds together all the rows of the same taxa
view(md)

#########################################################################

ggplot(md, aes(x = Taxa, y = freq))+
         geom_col() +
         geom_text(aes(label = freq), vjust = -0.3)+
         labs(y="Count", x="Taxon") #creates the graph