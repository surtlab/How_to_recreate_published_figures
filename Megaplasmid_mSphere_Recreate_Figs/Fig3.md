**Step 1: Prepare the Environment in R**

These commands will load various libraries into the active R environment. If you receive an error that one of these libraries isn't recognized, simply use the command install.packages()`

```
library(RCurl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
```

**Step 2: Acquire Data to Recreate Figure 1 from Github**

For this one we are going to build two different graphs and then combine them together at the end. So first, we need to load data for the two different graphs:

```
Fig3B_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Final_Putida_Cip_Results.csv"))
Fig3A_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Final_Putida_Nal_Results.csv"))
```
**Step 3: Create the Plot Data Order**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Then I'm using 'mutate' to reorder the data on the x-axis. I do this for both data sets for consistency, and call the data for each graph into a different variable, `p` or `q`

```
p<-Fig3A<-Fig3A_data %>% mutate(Strain=fct_relevel(Strain,"DBL305","DBL759","DBL1604","DBL1620"))
q<-Fig3B<-Fig3B_data %>% mutate(Strain=fct_relevel(Strain,"DBL305","DBL759","DBL1604","DBL1620"))

```

**Step 4: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variables for the graphs called
`p_plots` and `q_plots`, and use `geom_boxplot` and `geom_point` to populate the data in the graph. 
. 
Have set the theme as `theme_bw` with a bunch of modifications to the axis/legend fonts and sizes. These are all done through the `theme()` command using commands like `axis.text.x(element_text=()`. I get rid of grid lines using `theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()`. The lengend and axis titles can all be configured using commdands such as `axis.text.y=element(size=FONTSIZE, family=FONT_FAMILY)` as well as `legend.text=element_text()`.

I've also set the particular colors for each treatment using the `scale_fill_manual` command and by designating colors for each treatment as they appear after the reordering from above. Fill for this is set in the `aes` line above.
Then I use `facet_wrap` to bring together the data into one grid, and label the Y axis title using `ylab()`.
Lastly, I add significance values to each of the plots using `stat_compare_means`.

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
p_plots<-p%>%ggplot(aes(x=Strain, y=Normalized_Area,fill=Treatment))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=14)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=11,family = "Helvetica"),legend.position="none",axis.title.x=element_blank())+ylab("Halo Size Around Naldixic Acid Disc (Normalized)")
```
I do this for both sets of data
```
q_plots<-q %>%ggplot(aes(x=Strain, y=Normalized_Area,fill=Treatment))+geom_boxplot(width=0.4, outlier.shape=NA)+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))+theme_bw(base_size=14)+scale_fill_manual(values=c("grey","white"))+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=11,family = "Helvetica"),legend.position="none",axis.title.x=element_blank())+ylab("Halo Size Around Ciprofloxaxin Disc (Normalized)")
```

**Step 5: Combine Graphs Together into One Grid***

```
pq_plot<-plot_grid(Fig3A,Fig3B,labels=c("A","B"),label_size=20)
```

**Step 6: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
ggsave(file="~/Desktop/Fig3.png",plot=pq_plots,width=8,height=6)
ggsave(file="~/Desktop/Fig3.svg",plot=pq_plots,width=8,height=6)
```


**Step 6: Statistics**
