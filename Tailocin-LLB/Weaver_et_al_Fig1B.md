# A Basic Guide for Recreating Figure 1 from This Paper

**Step 1: Prepare the Environment in R**

These commands will load various libraries into the active R environment. If you receive an error that one of these libraries isn't recognized, simply use the command install.packages()`

```
library(RCurl)
library(forcats)
library(dplyr)
library(showtext)
library(ggplot2)
library(svglite)
```

**Step 2: Set Basic Formats For Figure**

This will load the theme for the graph as well as fonts for the axes. In particular the `font_add_google` function will This function will search the [Google Fonts Repository](https://fonts.google.com/) for a specified family name, download the proper font files, and then add them to sysfonts.

```
font_add_google("Poppins", "Poppins")
font_add_google("Roboto Mono", "Roboto Mono")
showtext_auto()
theme_set(theme_light(base_size=24, base_family="Poppins"))
```

**Step 3: Acquire Data to Recreate Figure 1 from Github**

The raw data underlying Figure 1 is openly available as a .csv file in Github. The RCurl function allows you to import this data from the webpage, however, you need to make sure that the address is correct and you are pulling in in the raw data.
We'll pull this raw data into a variable called `Fig1_data` through the use of the `read.csv` function.

```
Fig1B_data<-read.csv(text=getURL("https://raw.githubusercontent.com/surtlab/data_for_figures/master/LLB_Fig1B_data.csv"))
```


**Step 4: Create the Plot Background**

This step does a couple of different things at once. First, note the `%>%` function which acts as a pipe between commands. Next, note that I am reordering the factors on the x-axis to be shown in the order I want. You can alter the orders here in case you'd like to rearrange how the data is presented. Next I am labelling the axes and getting rid of the legend. All of these commands create a graph which I've called `q`.

*also important to note that if you are copy/pasting the lines below that the '+' cannot be on a new line or the code will not work*

```
q<-Fig1B_data%>%
mutate(Tailocin=fct_relevel(Tailocin,"Ptt-A","PttPEG-A","11-A","Ptt-B","PttPEG-B"))%>%
ggplot(aes(x=Tailocin,y=Area,color=Tailocin))
+theme(legend.position="none",axis.title=element_text(size=24),panel.grid=element_blank())+labs(x="Tailocin", y="Area(mm2)")
```

**Step 5: Add Boxplot and Data Points to the Graph**

I am a huge fan of showing each specific data point as well as a summary of the data using a boxplot. To do this, I create a new variable for the graph called `qbp_jitter_color`, and use `geom_boxplot` and `geom_jitter` to populate the data in the graph. I've set the color for the boxplot as a somewhat dark grey, and I've set the alpha value for geom_jitter (which determines the transparency of the data points). Within this command, I've also set the particular colors for each of the strain's data points through the `scale_color_manual` command and by designating colors for each strain as they appear after the reordering from above.

```
qbp_jitter_color<-q
+geom_boxplot(color = "gray50", outlier.alpha = 0)
+geom_jitter(aes(shape=Replicate),size = 2, alpha = 0.4, width = 0.2)
+scale_color_manual(values=c("blue","springgreen4","mediumpurple","grey38","red"))
```

**Step 6: Export Graph to a Figure File**

The last step here will be to export the graph to a readable figure file using the ggsave command. In this case, I will export as both `.png` files and `.svg` files on my desktop and called `Fig1.png` or `Fig1.svg`.

```
 ggsave(file="~/Desktop/Fig1B.png",plot=qbp_jitter_color,width=10,height=8)
 
 ggsave(file="~/Desktop/Fig1B.svg",plot=qbp_jitter_color,width=10,height=8)
```
**Step 6: Statistics**

OK, one figure at a time for the statistics. First up is Figure 4A, which compares wild type white/blue megaplasmid- strains as well as white megaplasmid- and blue megaplasmid+ strain in the P. putida background.

Set up the ANOVA framework:

```
Tailocin_size<-aov(Area~Tailocin+Tailocin:Replicate, data=Fig1B_data)
```
Check the summary stats of the ANOVA.

```
summary(Tailocin_size)
```
Which yields the following stats:
```
                   Df Sum Sq Mean Sq F value Pr(>F)    
Tailocin            4 1592.2   398.1   298.6 <2e-16 ***
Tailocin:Replicate  5  557.8   111.6    83.7 <2e-16 ***
Residuals          47   62.6     1.3                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
So there are significant effects for both the treatment (Tailocin) as well as a block assay effect. Now I'll set up a Post-hoc Tukey's HSD test to step through the Nal results
```
TukeyHSD(Tailocin_size)
```
which yields (there are values for the assay block effect too, but I'm not going to include those below)
```
 Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Area ~ Tailocin + Tailocin:Replicate, data = Fig1B_data)

$Tailocin
                         diff         lwr         upr     p adj
Ptt-A-11-A          9.6198168   8.2234469  11.0161868 0.0000000
Ptt-B-11-A         -0.9512992  -2.3182686   0.4156702 0.2943747
PttPEG-A-11-A      10.3152374   8.9188675  11.7116074 0.0000000
PttPEG-B-11-A      -1.4865542  -2.8535235  -0.1195848 0.0268581
Ptt-B-Ptt-A       -10.5711160 -11.9380854  -9.2041466 0.0000000
PttPEG-A-Ptt-A      0.6954206  -0.7009493   2.0917906 0.6227797
PttPEG-B-Ptt-A    -11.1063710 -12.4733404  -9.7394016 0.0000000
PttPEG-A-Ptt-B     11.2665366   9.8995673  12.6335060 0.0000000
PttPEG-B-Ptt-B     -0.5352550  -1.8721774   0.8016674 0.7869030
PttPEG-B-PttPEG-A -11.8017916 -13.1687610 -10.4348222 0.0000000

```
Which means that the size of the overlay of the Ptt tailocin (and PEG prepped version of this tailocin) against sensitivity class A strains signficantly differs by post hoc test from the tailocin from strain 011 against sensitivity class A strains!

There does not appear to be a significant difference of PEG treatment on the size of the overlay, and the size of the overlay for strain 011 tailocin against sensitivity class A is no statistically different than the size of the Ptt tailocin overlay against sensitivity class B. 

Lastly, we truly care about the comparison between the Ptt killing activity against sensitivity class A strains compared to all other killing activity (USA011 against A and Ptt against sensitivity class B combined). To set up this ttest we subset the data from these two classes
```
PttA<-subset(Fig1B_data, Ptt=="PttA",select=c(Area))
NonPttA<-subset(Fig1B_data, Ptt==NonPttA",select=c(Area)")
```
and then perform a ttest with these two groups against each other
```
t.test (PttA,NotPttA,var.equal=FALSE)
```
which yields a highly significant result regardless of possible p value corrections for multiple tests
```
	Welch Two Sample t-test

data:  PttA and NotPttA
t = 9.4536, df = 22.506, p-value = 2.681e-09
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  8.436461 13.170264
sample estimates:
mean of x mean of y 
19.237446  8.434083 
	Welch Two Sample t-test
```

