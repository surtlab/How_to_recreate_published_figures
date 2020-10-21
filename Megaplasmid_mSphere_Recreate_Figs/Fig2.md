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

We've already included ttest results in the graphs using the code above. However, below you'll find the code to pull these data and perform ttests independently. First we've got to pull the data for each strain into it's own variable:

```
DBL880<-subset(Fig2_data, Strain=="DBL880",select=c(Normalized_Area))
DBL907<-subset(Fig2_data, Strain=="DBL907",select=c(Normalized_Area))
DBL187<-subset(Fig2_data, Strain=="DBL187",select=c(Normalized_Area))
DAB282<-subset(Fig2_data, Strain=="DAB282",select=c(Normalized_Area))
DBL883<-subset(Fig2_data, Strain=="DBL883",select=c(Normalized_Area))
DBL910<-subset(Fig2_data, Strain=="DBL910",select=c(Normalized_Area))
DBL885<-subset(Fig2_data, Strain=="DBL885",select=c(Normalized_Area))
DBL912<-subset(Fig2_data, Strain=="DBL912",select=c(Normalized_Area))
DAB462<-subset(Fig2_data, Strain=="DAB462",select=c(Normalized_Area))
DAB895<-subset(Fig2_data, Strain=="DAB895",select=c(Normalized_Area))

```
Next, we actually perform the ttests like so:

```
> t.test(DBL880,DBL907,var.equal=FALSE)
```

Which yields the following result:

```
	Welch Two Sample t-test

data:  DBL880 and DBL907
t = -11.133, df = 26.414, p-value = 1.779e-11
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.9008174 -0.6202076
sample estimates:
mean of x mean of y 
 1.000011  1.760523 
```
And another significant difference! Just like was plotted for the graph!

> t.test(DAB282,DBL187,var.equal=FALSE)

	Welch Two Sample t-test

data:  DAB282 and DBL187
t = -9.5573, df = 45.812, p-value = 1.755e-12
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.5353400 -0.3490546
sample estimates:
mean of x mean of y 
 1.000032  1.442229 

> t.test(DBL883,DBL910,var.equal=FALSE)

	Welch Two Sample t-test

data:  DBL883 and DBL910
t = -13.514, df = 43.895, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.7128173 -0.5277935
sample estimates:
mean of x mean of y 
 1.000035  1.620341 

> t.test(DBL885,DBL912,var.equal=FALSE)

	Welch Two Sample t-test

data:  DBL885 and DBL912
t = -15.872, df = 21.134, p-value = 3.219e-13
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.4692532 -0.3605660
sample estimates:
mean of x mean of y 
 1.000010  1.414919 

> t.test(DAB462,DAB895,var.equal=FALSE)

	Welch Two Sample t-test

data:  DAB462 and DAB895
t = -4.1488, df = 52.633, p-value = 0.0001224
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1635243 -0.0569289
sample estimates:
mean of x mean of y 
 1.000033  1.110259 
