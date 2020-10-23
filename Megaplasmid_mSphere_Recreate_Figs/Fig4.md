**Step 1: Prepare the Environment in R**

These commands will load various libraries into the active R environment. If you receive an error that one of these libraries isn't recognized, simply use the command install.packages()`

```
library(RCurl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(ggpubr)
library(cowplot)
```

**Step 2: Acquire Data to Recreate Figure 1 from Github**

For this one we are going to build four graphs and then combine them together at the end. So first, we need to load data from two different sources. Then we'll pull the data from the second data set
into three different subsets for each figure B through D.

```
Fig4A_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Final_Putida_pMP_comps_Results.csv"))
Fig4BCD_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/mSphere_megaplasmid_2020/pMPPutida_cross-enviroment-relative_fitness.csv"))
Fig4B_data<-subset(Fig4BCD_data,Experiment=="Nal")
Fig4C_data<-subset(Fig4BCD_data,Experiment=="Cipro")
Fig4D_data<-subset(Fig4BCD_data,Experiment=="Temp")

```
**Step 3: Create the Plot Data Order**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Then I'm using 'mutate' to reorder the data on the x-axis. I do this for both data sets for consistency, and call the data for each graph into a different variable, `p` or `q`

```
p<- Fig4A_data %>% mutate(Strain_Pair=recode(Strain_Pair,'C305-1604'="DBL1604/DBL305",'D305-1620'="DBL1620/DBL305")) 
q<- ≈
r<-Fig4C_data %>% mutate(Treatment=recode(Treatment,'Cip 0 ng/uL'="0",'Cip 0.5 ng/uL'="0.5"))
s<-
```

**Step 4: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variables for the graphs called
`p_plots` and `q_plots`, and use `geom_boxplot` and `geom_point` to populate the data in the graph. 
. 
Have set the theme as `theme_bw` with a bunch of modifications to the axis/legend fonts and sizes. These are all done through the `theme()` command using commands like `axis.text.x(element_text=()`. I get rid of grid lines using `theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()`. The lengend and axis titles can all be configured using commdands such as `axis.text.y=element(size=FONTSIZE, family=FONT_FAMILY)` as well as `legend.text=element_text()`.

I've also set the particular colors for each treatment using the `scale_fill_manual` command and by designating colors for each treatment as they appear after the reordering from above. Fill for this is set in the `aes` line above.

I also place the legend inside of the first figure using `legend.position=c(.2,.1)`

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
p_plot<-p %>% ggplot(aes(x=Strain_Pair, y=Comp_Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=13,family = "Helvetica"),legend.position="none",axis.title.x=element_blank())+ylab("Relative Fitness")

q_plot<-q %>% ggplot(aes(x=Treatment, y=Relative.Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=13,family = "Helvetica"),legend.position="none",axis.title.x=element_blank())+ylab("   ")

r_plot<-r %>% ggplot(aes(x=Treatment, y=Relative.Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=12,family = "Helvetica"),axis.title.x=element_text(size=16, family ="Helvetica"),axis.title.y=element_text(size=16,family="Helvetica"),legend.position="none")+ylab("Relative Fitness")+xlab(expression("Ciprofloxaxin Concentration " * mu*"g/mL"))

s_plot<-s %>% ggplot(aes(x=Treatment, y=Relative.Fitness))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=16)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=13,family = "Helvetica"),legend.position="none",axis.title.x=element_blank())+ylab("   ")

```

**Step 5: Combine Graphs Together into One Grid***

```
pq_plot<-plot_grid(p_plots,q_plots,labels=c("A","B"),label_size=20)
```

**Step 6: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
ggsave(file="~/Desktop/Fig3.png",plot=pq_plot,width=8,height=6)
ggsave(file="~/Desktop/Fig3.svg",plot=pq_plot,width=8,height=6)
```


**Step 6: Statistics**
