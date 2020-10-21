# A Basic Guide for Recreating Figure 2 from This Paper

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

```
Fig2_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/Final_Data_Overlay_All_together.csv"))
```
**Step 3: Create the Plot Data Order**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Then I'm using 'mutate' to reorder the data on the x-axis. Then I'm loading the data into ggplot. Them I'm adding the boxplot (adding this before points so that it is behind points. 

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
q<-Fig2_data %>% mutate(Treatment=fct_relevel(Treatment,"wt","pMP"))


```

**Step 4: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variable for the graph called `q_plots`, and use `geom_boxplot` and `geom_point` to populate the data in the graph. 
. 
Have set the theme as `theme_bw` with a bunch of modifications to the axis/legend fonts and sizes. These are all done through the `theme()` command using commands like `axis.text.x(element_text=()`. I get rid of grid lines using `theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank()`. The lengend and axis titles can all be configured using commdands such as `axis.text.y=element(size=FONTSIZE, family=FONT_FAMILY)` as well as `legend.text=element_text()`.

I've also set the particular colors for each treatment using the `scale_fill_manual` command and by designating colors for each treatment as they appear after the reordering from above. Fill for this is set in the `aes` line above.
Then I use `facet_wrap` to bring together the data into one grid, and label the Y axis title using `ylab()`.
Lastly, I add significance values to each of the plots using `stat_compare_means`.

```
q_plots<-q%>%ggplot(aes(x=Treatment, y=Normalized_Area,fill=Treatment))
+geom_boxplot(width=0.4, outlier.shape=NA)
+geom_point(alpha=0.5,fill="grey",pch=21,position=position_jitter(width=0.11))
+theme_bw(base_size=16)
+scale_fill_manual(values=c("grey","white"))
+theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(), strip.text=element_text(size=10,family = "Helvetica"),axis.text.x=element_text(size=15,family = "Helvetica"),legend.position="none",axis.title.x=element_blank())
+facet_wrap(~Genomic_Background)
+ylab("Halo Size Around Naldixic Acid Discs (Normalized)")
+stat_compare_means(method="t.test",label="p.signif",label.y=2,label.x=0.75,size=5.,symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns")))
```
**Step 5: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
ggsave(file="~/Desktop/Fig2.png",plot=q_plots,width=8,height=6)
ggsave(file="~/Desktop/Fig2.svg",plot=q_plots,width=8,height=6)
```


**Step 6: Statistics**

In this case, we are going to be doing a bunch of ttests of the wild type and megaplasmid backgrounds for each strain. We've got to pull the data for each strain into it's own variable:

```
DBL880<-subset(Nal_Other, Strain=="DBL880",select=c(Normalized_Area))
DBL907<-subset(Nal_Other, Strain=="DBL907",select=c(Normalized_Area))
DBL187<-subset(Nal_Other, Strain=="DBL187",select=c(Normalized_Area))
DAB282<-subset(Nal_Other, Strain=="DAB282",select=c(Normalized_Area))
DBL883<-subset(Nal_Other, Strain=="DBL883",select=c(Normalized_Area))
DBL910<-subset(Nal_Other, Strain=="DBL910",select=c(Normalized_Area))
DBL885<-subset(Nal_Other, Strain=="DBL885",select=c(Normalized_Area))
DBL912<-subset(Nal_Other, Strain=="DBL912",select=c(Normalized_Area))
DAB462<-subset(Nal_Other, Strain=="DAB462",select=c(Normalized_Area))
DAB895<-subset(Nal_Other, Strain=="DAB895",select=c(Normalized_Area))

```
There is a statistically significant difference! The megaplasmid sensitizes P. stutzeri strain 28a24 (two versions of it) to Nalidixic acid according to Kirby-Bauer diffusion assays
