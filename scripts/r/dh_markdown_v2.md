---
title: "**Microbial biofilms facilitate the attachment and survival of harmful algal bloom (HAB)-causing phytoplankton on plastic debris in the marine environment**"
author: "*John Bennedick B. Quijano*"
date: "*12/28/2021*"
output: 
  html_document:
    toc: true
    toc_float: true
    toc_depth: 3
    fig_caption: true
    number_sections: true
    highlight: "haddock"
    keep_md: true
---


```{=html}
<script>
  addClassKlippyTo("pre.r, pre.markdown");
  addKlippy('right', 'top', 'auto', '1', 'Click to copy', 'Done');
</script>
```



<br>

# Preface {-}

This R Markdown document contains the pipeline used for the analysis of the high-throughput microbiome data of **Tolentino *et al*** (in prep). For publication details, see this [**github link**](https://github.com/jbquijano/dh-tol) (to be updated).

For a brief backgroud, this pipeline is limited to making alpha diversity plots, PCoA plots, PERMANOVA test, ANCOM-BC test and relative abundance taxonomy plots. These plots were made to look at the differences in the nano-microeukaryote distribution in plastic (polypropylene, PP; and polyethylene, PE) and seawater samples.

Initial data processing were mostly done in [**QIIME2**](https://qiime2.org/). QIIME2 scripts for data cleaning and initial analyses are included in this [**Jupyter notebook**](https://github.com/jbquijano/dh-tol) (to be updated).

<br>


# Installation
I assume that R is already installed in your system. RStudio is not required, but writing scripts will be easier when working with several functions and packages. The package **pacman** should be installed if you want ease in installing and updating the required packages. You can install pacman with:

```r
# Install from CRAN
install.packages('pacman')

# Load the packages
library('pacman') 
```

<br>

# Alpha diversity
As mentioned, initial analyses were already done in QIIME2 This includes the Shannon diversity, Pielou's evenness, Faith's phylogenetic diversity and observed features.


## Load the required packages
Before making plots, install the required packages with:

```r
# Install and load the package qiime2R
p_load_gh('jbisanz/qiime2R') # For importing .qza files to R

# Install and load the following packages
p_load('tidyverse') # For data cleaning
p_load('ggsci') # For choosing color templates
```


## Import data & initial results
The metadata which contains the important environmental data and sample information can be imported with the following:



```r
# Import metadata
metadata <- read.table(
  file = "data/metadata.tsv", 
  sep = "\t", header = TRUE)

# Make a dataframe with environment type metadata only
env <- metadata %>% 
  select(
    samplename, 
    environment1, 
    environment)
```

On the other hand, the alpha diversity measures from the initial analysis can be imported with:



```r
# Observed features
obs.alpha <- read_qza("data/alpha_div/observed_features_vector.qza")$dat %>%
  rownames_to_column("samplename")

# Shannon diversity index
shannon.alpha <- read_qza("data/alpha_div/shannon_vector.qza")$data %>%
  rownames_to_column("samplename")

# Faith's phylogenetic diversity
faith.alpha <- read_qza("data/alpha_div/faith_pd_vector.qza")$data %>%
  rownames_to_column("samplename")

# Pielou's evenness
even.alpha <- read_qza("data/alpha_div/evenness_vector.qza")$data %>%
  rownames_to_column("samplename")

# Merge the alpha diversity measure dataframes
merged <-
  Reduce(left_join, 
    list(even.alpha, 
      faith.alpha,
      obs.alpha,
      shannon.alpha))

# Transform the table to a long table
merged.long <- merged %>%
  pivot_longer(
    !samplename, 
    names_to = "alphametric",
    values_to = "alphavalue") %>%
  left_join(env)
```

You can view the data frame containing the transformed data by calling the variable name (as shown below) or by double-clicking the variable name stored in the `Environment` pane (if you are using RStudio)

```r
merged.long
```

The final data frame should look like this:

```{=html}
<div id="htmlwidget-71c7fa6f0a9b207c5e43" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-71c7fa6f0a9b207c5e43">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.2.1. <\/strong>\n  Final data table of merged alpha diversity measure and metadata.\n<\/caption>","data":[["DH10","DH10","DH10","DH10","DH11","DH11","DH11","DH11","DH12","DH12","DH12","DH12","DH13","DH13","DH13","DH13","DH14","DH14","DH14","DH14","DH15","DH15","DH15","DH15","DH16","DH16","DH16","DH16","DH2","DH2","DH2","DH2","DH3","DH3","DH3","DH3","DH6","DH6","DH6","DH6","DH8","DH8","DH8","DH8","DH9","DH9","DH9","DH9"],["pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy","pielou_evenness","faith_pd","observed_features","shannon_entropy"],[0.692522487895018,6.9137605656121,46,3.82519086785088,0.671884392553995,5.45763685770561,37,3.50015041010397,0.496107521809405,5.03647672281384,34,2.52392858247014,0.595992816784962,6.06239611493385,59,3.50601300112964,0.714954564725889,10.3511366204359,54,4.11448308915574,0.755422091668014,7.56537475670107,58,4.42524825629074,0.719259627551504,7.13863062154029,67,4.36309305182286,0.609338987347436,5.55391753598521,25,2.82968262806445,0.668189121320275,6.93211415406292,64,4.00913472792165,0.599444583426417,3.71612126028322,31,2.96976614349255,0.800071844049895,9.09314695452091,50,4.51549042930548,0.799984487091207,9.57900905740008,75,4.98295832530047],["Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Free-living","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached"],["PP","PP","PP","PP","PE","PE","PE","PE","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","SW","PP","PP","PP","PP","PE","PE","PE","PE","PP","PP","PP","PP","PE","PE","PE","PE","PE","PE","PE","PE"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>samplename<\/th>\n      <th>alphametric<\/th>\n      <th>alphavalue<\/th>\n      <th>environment1<\/th>\n      <th>environment<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":2}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

## Alpha diversity plots
In this paper, we used violin plots to compare the different alpha diversity measures between the different environment types (i.e., plastic and seawater samples).


### Make alpha diversity violin plots
To make a violin plot with a box plot, run the code below:

```r
ggplot(
  merged.long,
  aes(x = environment1, 
      y = alphavalue, 
      fill = environment1)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white")
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/intial plot-1.png" alt="**Figure 2.3.1.1.** Initial alpha diversity violin plots comparing free-living and plastics-attached nano-microeukaryote communties."  />
<p class="caption">**Figure 2.3.1.1.** Initial alpha diversity violin plots comparing free-living and plastics-attached nano-microeukaryote communties.</p>
</div>

The plot we made, however, is uninformative for several reasons. First, if you are going to review table 2.2.1, four alpha diversity metrics (column `alphametric`) were indicated. As you can observe, the plot does not inform about the alpha diversity metric used. It lumped all the diversity metrics, only discriminating between the free-living and plastics-attached.

Second, we also want to check if plastic type contributes to a difference in alpha diversity, but the plot only grouped the samples to free-living and plastics-attached.

Lastly, the figure contains wrong labels and is not visually pleasing.

### Customize the violin plots
To address these issues, you can run the code below to make basic violin-box plots per alpha diversity metric. This makes a plot for the `environment1` group (free-living vs plastics-attached).

```r
# Assign labels to each alpha diversity measure as it will be grouped using the facet_wrap function
mylab <- c(
  `faith_pd` = "Faith's PD", 
  `observed_features` = "Observed ASVs", 
  `pielou_evenness` = "Evenness", 
  `shannon_entropy` = "Shannon"
  )

# Make 4 violin plots put side-by-side
ggplot(merged.long,
  aes(x = environment1, 
      y = alphavalue, 
      fill = environment1)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white") + 
  facet_wrap(vars(alphametric), 
    scales = "free_y", 
    ncol = 4, nrow = 1,
    labeller = as_labeller(mylab)) # Group alpha diversity indices accordingly
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/fig-custom-1-1.png" alt="**Figure 2.3.2.1.** Violin plots of the four commonly used alpha diversity metrics comparing free-living and plastics-attached nano-microeukaryote communties."  />
<p class="caption">**Figure 2.3.2.1.** Violin plots of the four commonly used alpha diversity metrics comparing free-living and plastics-attached nano-microeukaryote communties.</p>
</div>

To make another violin plot for another grouping [i.e., `environment` (seawater vs polyethylene vs polypropylene)] and put it on the side of the plot above, you can use this:

```r
# Make plots based on environment type (same with the plot above)
p1 <- ggplot(merged.long,
  aes(x = environment1, 
      y = alphavalue, 
      fill = environment1)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white") + 
  facet_wrap(vars(alphametric), 
    scales = "free_y", 
    ncol = 4, nrow = 1,
    labeller = as_labeller(mylab))

# Make new violin plots based on substrate type
p2 <- ggplot(merged.long,
  aes(x = factor(environment, level = c('SW','PE','PP')), 
      y = alphavalue, 
      fill = environment,)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white") + 
  facet_wrap(vars(alphametric), 
    scales = "free_y", 
    ncol = 4, nrow = 1,
    labeller = as_labeller(mylab))

# Arrange and put plots side-by-side
ggpubr::ggarrange(p1, p2,
                  ncol = 2, nrow = 1,
                  legend = "bottom")
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/fig-side-2-1.png" alt="**Figure 2.3.2.2.** Violin plots of the four commonly used alpha diversity metrics comparing free-living and plastics-attached nano-microeukaryote communties based on (**left**) environment type and (**right**) substrate type."  />
<p class="caption">**Figure 2.3.2.2.** Violin plots of the four commonly used alpha diversity metrics comparing free-living and plastics-attached nano-microeukaryote communties based on (**left**) environment type and (**right**) substrate type.</p>
</div>

To make customizations based on your preference, you can consult the documentation/tutorials in [ggplot2]( https://ggplot2.tidyverse.org/) website. In this publication, however, we used:

```r
mycol <- pal_npg()(3)

# Customize left plot
p1 <- p1 +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.title = element_blank()
  ) + # Remove unnecessary labels
  labs (y = "Alpha Diversity Measure") + 
  scale_fill_manual(
    values = mycol,
    breaks = c("Free-living", "Plastics-attached"),
    labels = c("Seawater", "Plastic")
    ) + # Edit labels by changing the "breaks" and "labels" parameters
  guides(fill = guide_legend(title='Environment'))

# Customize right plot
p2 <- p2 +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.title = element_blank()
  ) +
  labs (y = "Alpha Diversity Measure") + 
  scale_fill_manual(
    values = mycol,
    breaks = c("SW", "PE", "PP"),
    labels = c("Seawater", "Polyethylene", "Polypropylene")
    ) +
  guides(fill = guide_legend(title='Environment'))

ggpubr::ggarrange(p1, NULL, p2, 
                  labels = c("A", "", "B"),
                  ncol = 3, nrow = 1,
                  widths = c(0.45, 0.025, 0.5), # Adjust the distance between plots
                  legend = "bottom")
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/fig-side-3-1.png" alt="**Figure 2.3.2.3.** Final violin plots of the four commonly used alpha diversity metrics comparing free-living and plastics-attached nano-microeukaryote communties based on (**A**) environment type and (**B**) substrate type."  />
<p class="caption">**Figure 2.3.2.3.** Final violin plots of the four commonly used alpha diversity metrics comparing free-living and plastics-attached nano-microeukaryote communties based on (**A**) environment type and (**B**) substrate type.</p>
</div>

## Statistical testing
Some statistical tests were also done in QIIME2. You can import the data with:

```r
envi_even <- read_csv("data/alpha_sig/envi/even.csv") %>%
    mutate_if(is.numeric, ~round(., 4)) # Round off to 4 decimal places
```

The table should look like the one presented below. As you can see, the Kruskal-Wallis test results showed high p-value and adjusted p-values (q-value). You can consult [Kruskal and Wallis (1952)]( https://www.jstor.org/stable/2280779) for interpretation of the results.

```{=html}
<div id="htmlwidget-6b444ebbef94781091da" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6b444ebbef94781091da">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.1. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the Pielou's evenness based on substrate type.\n<\/caption>","data":[["PE (n=4)","PE (n=4)","PP (n=3)"],["PP (n=3)","SW (n=5)","SW (n=5)"],[2,0.96,0.2],[0.1573,0.3272,0.6547],[0.4719,0.4908,0.6547]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

The difference of other alpha diversity measures among the substrate types was also tested using the Kruskal-Wallis test. The data can be imported with:

```r
envi_fait <- read_csv("data/alpha_sig/envi/fait.csv") %>%
    mutate_if(is.numeric, ~round(., 4))

envi_obse <- read_csv("data/alpha_sig/envi/obse.csv") %>%
    mutate_if(is.numeric, ~round(., 4))

envi_shan <- read_csv("data/alpha_sig/envi/shan.csv") %>%
    mutate_if(is.numeric, ~round(., 4))
```

The results can be seen in the tables below

```{=html}
<div id="htmlwidget-435fcdd7fb69d18c8c93" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-435fcdd7fb69d18c8c93">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.2. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the Faith's phylogenetic diversity based on substrate type.\n<\/caption>","data":[["PE (n=4)","PE (n=4)","PP (n=3)"],["PP (n=3)","SW (n=5)","SW (n=5)"],[2,0.06,1.8],[0.1573,0.8065,0.1797],[0.2696,0.8065,0.2696]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-1244af8f7c600a0740a1" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-1244af8f7c600a0740a1">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.3. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the observed features based on substrate type.\n<\/caption>","data":[["PE (n=4)","PE (n=4)","PP (n=3)"],["PP (n=3)","SW (n=5)","SW (n=5)"],[3.125,0.06,3.7556],[0.0771,0.8065,0.0526],[0.1156,0.8065,0.1156]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-f6d1a5964d87655391f0" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f6d1a5964d87655391f0">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.4. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the Shannon diversity metric based on substrate type.\n<\/caption>","data":[["PE (n=4)","PE (n=4)","PP (n=3)"],["PP (n=3)","SW (n=5)","SW (n=5)"],[3.125,0.54,1.0889],[0.0771,0.4624,0.2967],[0.2313,0.4624,0.4451]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

In addition, the difference of alpha diversity measures between the environment types was also tested using the Kruskal-Wallis test. The data can be imported with:

```r
envi1_even <- read_csv("data/alpha_sig/envi1/even.csv") %>%
    mutate_if(is.numeric, ~round(., 4))

envi1_fait <- read_csv("data/alpha_sig/envi1/fait.csv") %>%
    mutate_if(is.numeric, ~round(., 4))

envi1_obse <- read_csv("data/alpha_sig/envi1/obse.csv") %>%
    mutate_if(is.numeric, ~round(., 4))

envi1_shan <- read_csv("data/alpha_sig/envi1/shan.csv") %>%
    mutate_if(is.numeric, ~round(., 4))
```

The results can be seen in the tables below:

```{=html}
<div id="htmlwidget-25f37418ae5ce3c9bc24" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-25f37418ae5ce3c9bc24">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.5. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the Pielou's evenness based on environment type.\n<\/caption>","data":[["Free-living (n=5)"],["Particle-attached (n=7)"],[0.1648],[0.6847],[0.6847]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-ef8a0fab9491bf74938a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-ef8a0fab9491bf74938a">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.6. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the Faith's phylogenetic diversity based on environment type.\n<\/caption>","data":[["Free-living (n=5)"],["Particle-attached (n=7)"],[0.3231],[0.5698],[0.5698]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-51df3ed4543a2561da3a" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-51df3ed4543a2561da3a">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.7. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the observed features based on environment type.\n<\/caption>","data":[["Free-living (n=5)"],["Particle-attached (n=7)"],[0.7978],[0.3718],[0.3718]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-bc169def22f75e95df07" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bc169def22f75e95df07">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 2.4.8. <\/strong>\n  Data table of the Kruskal-Wallis one-way ANOVA results for the Shannon diversity metric based on environment type.\n<\/caption>","data":[["Free-living (n=5)"],["Particle-attached (n=7)"],[0.0066],[0.9353],[0.9353]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Group 1<\/th>\n      <th>Group 2<\/th>\n      <th>H<\/th>\n      <th>p-value<\/th>\n      <th>q-value<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

<br>

# Beta diversity
Beta diversity metrics were also computed in QIIME2. The metrics used were Jaccard index, Bray-Curtis dissimilarity, unweighted and weighted Unifrac distances.

## Load the required packages
In this part, statistical tests will be done outside QIIME2 and inside R. To load the required libraries, run the code below:

```r
# Install and load the packages for ecological statistics
pacman::p_load('vegan', 'phyloseq')
```

## Import data & initial results
The metrices were computed in QIIME2 as well as the PCoA vectors. These data should also be imported. I assume that you have followed the tutorial from the start, thus, the library loading and importing of other files are not included in the succeeding codes. Please refer to the previous sections if you are having trouble with proceeding. You can import the required files with:

```r
# Import PCoA data from QIIME2 artifact
bc_pcoa_data <- read_qza(
  "data/beta_pcoa/bray_curtis_pcoa_results.qza")
jac_pcoa_data <- read_qza(
  "data/beta_pcoa/jaccard_pcoa_results.qza")
unwuni_pcoa_data <- read_qza(
  "data/beta_pcoa/unweighted_unifrac_pcoa_results.qza")
wtuni_pcoa_data <- read_qza(
  "data/beta_pcoa/weighted_unifrac_pcoa_results.qza")

# Extract PCoA data from the file
bc_pcoa_vector <- bc_pcoa_data$data$Vectors %>% 
  select(SampleID, PC1, PC2)

# Extract important columns from the metadata
env_pcoa <- metadata %>% 
  select(
    SampleID,
    cluster,
    samplename, 
    environment1, 
    environment)

# Merge metadata and PCoA data
bc_pcoa_full <- bc_pcoa_vector %>%
  left_join(env_pcoa)
```
In the code below, I selected the Bray-Curtis dissimilarity metric as an example. If you followed the same code and no errors occurred, the data frame should look like this:

```{=html}
<div id="htmlwidget-6233de8ac47b01b181f9" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-6233de8ac47b01b181f9">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.2.1. <\/strong>\n  Final data table of Bray-Curtis  measure and metadata.\n<\/caption>","data":[["DH10","DH11","DH12","DH13","DH14","DH15","DH16","DH2","DH3","DH6","DH8","DH9"],[-0.38723961358561,-0.267921577563932,0.378169465103275,0.475125614173056,0.42278986559083,0.451483941596965,0.352961883145723,-0.379481707698388,-0.077140015311253,-0.416464981479033,-0.374833235930041,-0.177449638041593],[-0.331883222003699,-0.553099957423684,0.0305470106910111,0.0306144067449422,-0.02615889739435,0.0218954777074648,-0.00521986412218274,0.264586317894408,0.178390393847774,0.270313923152816,0.262165673346718,-0.142151262441221],["OB","OB","OB","OB","IB","IB","IB","IB","IB","IB","OB","OB"],["DH10","DH11","DH12","DH13","DH14","DH15","DH16","DH2","DH3","DH6","DH8","DH9"],["Plastics-attached","Plastics-attached","Free-living","Free-living","Free-living","Free-living","Free-living","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached","Plastics-attached"],["PP","PE","SW","SW","SW","SW","SW","PP","PE","PP","PE","PE"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>SampleID<\/th>\n      <th>PC1<\/th>\n      <th>PC2<\/th>\n      <th>cluster<\/th>\n      <th>samplename<\/th>\n      <th>environment1<\/th>\n      <th>environment<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[1,2]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

## PCoA plots
Beta diversity metrics are usually viewed using ordination plots. In this example, we used a Principal Coordinates Analysis (PCoA) plots.

### Make PCoA plots
To make a PCoA plot, which is usually in the form of scatter plots, use the code below:

```r
ggplot(bc_pcoa_full, aes(PC1, PC2)) + 
geom_point(aes(fill = environment1, shape = cluster),
           size=5, 
           alpha = 0.8)
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/initial pcoa plot-1.png" alt="**Figure 3.3.1.1.** Initial PCoA plot for Bray-Curtis Dissimilarity."  />
<p class="caption">**Figure 3.3.1.1.** Initial PCoA plot for Bray-Curtis Dissimilarity.</p>
</div>

### Customize the PcoA plot
To make the plot less cluttered and simple, you can customize the plot as follows:

```r
# Get data for grouping based on environment type (used for geom_polygon)
bc_pcoa_env1 <- bc_pcoa_full %>%
  group_by(environment1) %>%
  slice(chull(PC1, PC2))

# Get the proportioned explained
bc_pc1 <- bc_pcoa_data$data$ProportionExplained[[1]] * 100
bc_pc2 <- bc_pcoa_data$data$ProportionExplained[[2]] * 100

# Select color and shapes
pcoa_colors <- pal_npg()(2)
pcoa_shapes <- c(21, 24)

# Make a placeholder for the plot data
p_bc_pcoa_env1 <- ggplot(bc_pcoa_full, 
                         aes(PC1, PC2)) + 
  geom_point(aes(fill = environment1, shape = cluster),
             size=5, 
             alpha = 0.8)

p_bc_pcoa_env1 +
  theme_bw() + # Select ggplot theme
  geom_vline(xintercept = 0, 
             linetype = 'dashed', 
             col = 'gray') + # For y-axis
  geom_hline(yintercept = 0, 
             linetype = 'dashed',
             col = 'gray') + # For x-axis
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) + # For axis title mods
  labs(cluster = "Cluster",
       environment1 = "Environment") +
  xlab(bc_pc1 %>%
         round(2) %>%
         toString() %>%
         paste("PC1 ", "(", ., "%", ")", sep = "")) + # For x-axis label
  ylab(bc_pc2 %>%
         round(2) %>%
         toString() %>%
         paste("PC2 ", "(", ., "%", ")", sep = "")) + # For y-axis label
  geom_polygon(data = bc_pcoa_env1, 
               alpha = 0.35, 
               aes(x = PC1, 
                   y = PC2, 
                   fill = environment1, 
                   group = environment1), 
               show.legend = FALSE) + # For polygon grouping
  scale_fill_manual(values = pcoa_colors,
                    breaks = c("Free-living", "Plastics-attached"),
                    labels = c("Seawater", "Plastic"),
                    "Environment",
                    guide = guide_legend(
                      override.aes = list(shape=15,
                                          color=pcoa_colors))) + # Modify legend
  scale_shape_manual(values = pcoa_shapes,
                     breaks = c("IB", "OB"),
                     labels = c("Inner bay", "Outer bay"),
                     "Cluster",
                     guide = guide_legend(
                       override.aes = list(shape=c(21,24),
                                           fill="white",
                                           color="black"))) # Modify legend
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/pcoa plot custom 1-1.png" alt="**Figure 3.3.2.1.** Bray-Curtis Dissimilarity PCoA plot grouped by environment type and clustering based on environmental parameters."  />
<p class="caption">**Figure 3.3.2.1.** Bray-Curtis Dissimilarity PCoA plot grouped by environment type and clustering based on environmental parameters.</p>
</div>









You can repeat the codes above and change the parameters to pick a different beta diversity metric, I made plots for other metrics (code not shown here) and merged them to get the plot below:
<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/pcoa plot 4-1.png" alt="**Figure 3.3.2.2.** PCoA plots of different beta diversity metrics grouped by environment type and clustering based on environmental parameters."  />
<p class="caption">**Figure 3.3.2.2.** PCoA plots of different beta diversity metrics grouped by environment type and clustering based on environmental parameters.</p>
</div>

## Statistical testing
You can either import the beta diversity metrics from QIIME2 or compute it in R, there is no difference between the results as long as you are using the same OTU/ASV table and other related information. In this portion, we will compute the beta diversity metrics in R.

We will import all related data from QIIME2 to R as a phyloseq object. You can do this by using this code:

```r
# Import metadata for plastic vs seawater communities
env1_perm <- read.table(
  file = "data/metadata.tsv", 
  sep = "\t", header = TRUE) %>%
  arrange(samplename)

# Filter metadata to include PE and PP communities only
env_perm <- env1_perm %>%
  filter(environment != 'SW')

# To compare plastic and seawater communities
physeq_env1 <- qza_to_phyloseq(
  features = "data/rarefied_table.qza",
  tree = "data/rooted-tree.qza",
  taxonomy = "data/asv-taxa-250.qza",
  metadata = "data/metadata.tsv")

# To compare PE vs PP communities
physeq_env <- subset_samples(physeq_env1, environment != "SW")
```

You can compute the distance metrics using the code below:

```r
# Compute distance metrics to compare plastic and seawater samples
jac_dist_env1 <- distance(physeq_env1, method = "jaccard")
bc_dist_env1 <- distance(physeq_env1, method = "bray")
unwuni_dist_env1 <- distance(physeq_env1, method = "unifrac")
wtuni_dist_env1 <- distance(physeq_env1, method = "wunifrac")

# Compute distance metrics to compare PE and PP samples
jac_dist_env <- distance(physeq_env, method = "jaccard")
bc_dist_env <- distance(physeq_env, method = "bray")
unwuni_dist_env <- distance(physeq_env, method = "unifrac")
wtuni_dist_env <- distance(physeq_env, method = "wunifrac")
```

PERMANOVA test comparing plastic and seawater communities can be run using:

```r
# Run PERMANOVA based on environment type (plastic vs seawater)
jac_env1_permanova <- adonis2(jac_dist_env1 ~ environment1, 
        data = env1_perm, 
        permutations = 999)

bc_env1_permanova <- adonis2(bc_dist_env1 ~ environment1, 
        data = env1_perm, 
        permutations = 999)

unwuni_env1_permanova <- adonis2(unwuni_dist_env1 ~ environment1, 
        data = env1_perm, 
        permutations = 999)

wtuni_env1_permanova <- adonis2(wtuni_dist_env1 ~ environment1, 
        data = env1_perm, 
        permutations = 999)
```

The tables for plastic-seawater PERMANOVA test should look like these

```{=html}
<div id="htmlwidget-3c795c750cbc423ba7b0" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-3c795c750cbc423ba7b0">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.1. <\/strong>\n  PERMANOVA test results of Jaccard distance comparing plastics and seawater communities.\n<\/caption>","data":[[1,10,11],[1.20115129679497,3.51896453973708,4.72011583653205],[0.254474961715659,0.745525038284341,1],[3.41336573083147,null,null],[0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-4ec515a7e94e95fbaba3" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-4ec515a7e94e95fbaba3">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.2. <\/strong>\n  PERMANOVA test results of Bray-Curtis dissimilarity comparing plastics and seawater communities.\n<\/caption>","data":[[1,10,11],[1.51831124313204,2.79096755059386,4.3092787937259],[0.352335347934004,0.647664652065996,1],[5.44008920063932,null,null],[0.002,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-e88a01e682e976027ec7" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-e88a01e682e976027ec7">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.3. <\/strong>\n  PERMANOVA test results of unweighted Unifrac comparing plastics and seawater communities.\n<\/caption>","data":[[1,10,11],[0.819412150087299,1.588621212741,2.4080333628283],[0.340282723128419,0.659717276871581,1],[5.15800836294692,null,null],[0.001,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-9fcf0798488e3c1c887f" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-9fcf0798488e3c1c887f">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.4. <\/strong>\n  PERMANOVA test results of weighted Unifrac comparing plastics and seawater communities.\n<\/caption>","data":[[1,10,11],[0.0526092038878806,0.0439645796636355,0.096573783551516],[0.544756578371157,0.455243421628843,1],[11.9662701862234,null,null],[0.002,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

For PERMANOVA test comparing PE and PP communities, the test can done using:

```r
# Run PERMANOVA based on substrate type (PE vs PP)
jac_env_permanova <- adonis2(jac_dist_env ~ environment, 
        data = env_perm, 
        permutations = 999)

bc_env_permanova <- adonis2(bc_dist_env ~ environment, 
        data = env_perm, 
        permutations = 999)

unwuni_env_permanova <- adonis2(unwuni_dist_env ~ environment, 
        data = env_perm, 
        permutations = 999)

wtuni_env_permanova <- adonis2(wtuni_dist_env ~ environment, 
        data = env_perm, 
        permutations = 999)
```

Table of PERMANOVA test results for PE-PP comparison should look like:

```{=html}
<div id="htmlwidget-831feb38b8ec9624f069" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-831feb38b8ec9624f069">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.5. <\/strong>\n  PERMANOVA test results of Jaccard distance comparing PE and PP communities.\n<\/caption>","data":[[1,5,6],[0.445717897905598,1.89909531807718,2.34481321598278],[0.190086739049185,0.809913260950815,1],[1.17350059700238,null,null],[0.323,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-bc75de77a9eaa4716a1b" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-bc75de77a9eaa4716a1b">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.6. <\/strong>\n  PERMANOVA test results of Bray-Curtis dissimilarity comparing PE and PP communities.\n<\/caption>","data":[[1,5,6],[0.417299829793342,1.58116895599186,1.9984687857852],[0.208809781149213,0.791190218850787,1],[1.31959278675432,null,null],[0.287,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-998e5db42799dc68f9ae" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-998e5db42799dc68f9ae">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.7. <\/strong>\n  PERMANOVA test results of unweighted Unifrac comparing PE and PP communities.\n<\/caption>","data":[[1,5,6],[0.19063689229384,0.862071890334827,1.05270878262867],[0.181091765775726,0.818908234224274,1],[1.10569022393131,null,null],[0.332,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


```{=html}
<div id="htmlwidget-f02e177715e8e2421aac" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-f02e177715e8e2421aac">{"x":{"filter":"none","vertical":false,"caption":"<caption style=\"text-align: center;\">\n  <strong>Table 3.4.8. <\/strong>\n  PERMANOVA test results of weighted Unifrac comparing PE and PP communities.\n<\/caption>","data":[[1,5,6],[0.0100760150150939,0.0287590169294313,0.0388350319445252],[0.259456848895791,0.740543151104209,1],[1.75180101597672,null,null],[0.203,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"dom":"t","columnDefs":[{"className":"dt-right","targets":[0,1,2,3,4]}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

# Differential abundance testing
Differential abundance testing will also be done here. There are many differential abundance tests, but we chose Analysis of Compositions of Microbiomes with Bias Correction (ANCOM-BC) of the package *ANCOMBC* (Lin and Peddada, 2020). This test aims to determine the major contributors to the observed dissimilarity in the community composition.

## Load the required packages
To load the essential packages, run the code below:

```r
pacman::p_load(ANCOMBC, microbiome)
```

## Data import and extraction
To run ANCOMBC, a phyloseq data should be imported first. You can do that with this code:

```r
# Prepare feature/ASV/OTU table
asv_short <-
  read.table(file = "data/table-short-dist.csv",
    sep = ",",
    header = TRUE) %>%
  column_to_rownames("asvidshort")

# Prepare taxonomy mapping file
tax_short <- 
  read.table(file = 
    "data/table-short-taxa.csv",
    sep = ",",
    header = TRUE) %>% 
  separate(taxonomy,
    sep = ";", 
    c("Domain", "Kingdom", "Phylum",
    "Class", "Order", "Family",
    "Genus", "Species")) %>%
  mutate(Species = str_c("ASV-", asvidshort)) %>%
  column_to_rownames("asvidshort") %>%
  mutate(Domain = str_sub(Domain, 4, -1)) %>%
  mutate(across(!Domain, ~str_sub(., 5, -1))) %>%
  mutate_all(~replace_na(., "Unclassified")) %>%
  mutate_all(~str_replace(., "_", " ")) %>%
  select(!Kingdom)

# Prepare the metadata
metadata <- read.table(
  file = "data/metadata.tsv", 
  sep = "\t", header = TRUE) %>% 
  select(
    samplename, 
    environment1, 
    environment,
    cluster)  %>%
  column_to_rownames("samplename") %>%
  mutate_all(~str_replace(., "-", "")) %>%
  sample_data()

# Import ASV table & taxonomy mapping file as phyloseq data
asvmat <-data.matrix(asv_short)
taxmat <- as.matrix(tax_short, mode = 'character')

# Modify the phyloseq data
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)

# Merge phyloseq data (feature and taxonomy mapping file)
physeq <- phyloseq(ASV, TAX)

# Merge with metadata
physeq_meta = merge_phyloseq(
  physeq, metadata)

# Select taxonomy level
physeq_phylum = aggregate_taxa(physeq_meta, "Phylum")
```

You can run ANCOMBC and extract important information with:

```r
# Run ANCOMBC
ancombc_data <-
  ancombc(phyloseq = physeq_phylum, 
    formula = "environment1", 
    p_adj_method = "holm", zero_cut = 0.90, 
    lib_cut = 0, group = "environment1", 
    struc_zero = TRUE, neg_lb = TRUE, 
    tol = 1e-5, max_iter = 100, 
    conserve = TRUE, alpha = 0.05, global = TRUE)

# Extract important information
results <- ancombc_data$res

# Rename headers and make placeholders
# Coefficients
tab_coef <- results$beta
colnames(tab_coef) <- "Coefficient"
# Standrd errors
tab_se <- results$se
colnames(tab_se) <- "Standard errors"
# Test statistics
tab_w <- results$W
colnames(tab_w) <- "W"
# p-values
tab_p <- results$p_val
colnames(tab_p) <- "p-value"
# Adjusted p-values
tab_q <- results$q
colnames(tab_q) <- "q-value"
# Differentially abundant taxa
tab_diff <- results$diff_abn
colnames(tab_diff) <- "Differentially abundant?"

# Merge data into a table
results_ancombc <- bind_cols(
  tab_coef, 
  tab_diff, 
  tab_p, 
  tab_q, 
  tab_se, 
  tab_w)

# Export data as a .csv file
write.csv(results_ancombc, "ancom_results.csv")
```

ANCOMBC includes bias-correction. You can do this by:

```r
# Bias-adjusted abundances
samp_frac = ancombc_data$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(abundances(physeq_phylum) + 1) 
# Adjust the log observed abundances
log_obs_abn_adj = t(t(log_obs_abn) - samp_frac) %>%
  as_tibble(rownames = NA)

# Export as a .csv file
write.csv(log_obs_abn_adj, "ancom-results-adjabn.csv")
```

## Plot ANCOMBC heatmap
A heatmap based on the bias-adjusted ANCOM results can be done with the code below:

```r
# Make a long table
abn_adj_long <- log_obs_abn_adj %>%
  rownames_to_column("Phylum") %>%
  pivot_longer(!Phylum, 
    names_to = "samplename", 
    values_to = "log_adj_abundance") %>%
  arrange(desc(log_adj_abundance))

# Adjust position of samples based on clustering 
abn_adj_long$samplename <- factor(abn_adj_long$samplename,
  levels = c("DH12", "DH15", "DH14", "DH13", "DH16",
    "DH3", "DH9", "DH11", "DH10", "DH2",
    "DH8", "DH6"))

# Plot and customize
ggplot(abn_adj_long,
  aes(x = Phylum, y = samplename, fill = log_adj_abundance)) +
  geom_tile(color = "white",
    lwd = 0.5,
    linetype = 1) +
  coord_fixed() +
  scale_fill_distiller(palette = "YlOrRd", 
    limits = c(-0.5, 8),
    direction = 1) +
  labs(fill = element_blank()) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, 
          vjust = 0.5, 
          hjust=1))
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/ancombc heatmap-1.png" alt="**Figure 4.3.1.** Heatmap of log bias-adjusted abundance of identified phyla."  />
<p class="caption">**Figure 4.3.1.** Heatmap of log bias-adjusted abundance of identified phyla.</p>
</div>
You can check the `ancom_results.csv` to check for the adjusted p-values and other important statistical information. The figure can be manually modified to add additional info.

# Taxonomy barplots
Distribution of taxa among sites can be seen with bar plots. In this part, a perecent stacked barplot will be made.

## Load the required packages
Import the essential packages before proceeding. Load the packages using:

```r
pacman::p_load('RColorBrewer')
```

## Data import and cleaning
To make the plot, the data should be imported and cleaned. To do this, use the code below:

```r
# Import and transform to a long table
taxa_bplot_data <- read.table(file = "data/tax_bplot_phylum.csv",
  sep = ",", header = TRUE) %>%
  select(!Total) %>%
  pivot_longer(!Phylum,
    names_to = "samplename",
    values_to = "count")

# Arrange order of taxa
taxa_bplot_data$Phylum <- factor(taxa_bplot_data$Phylum,
  levels = c("Ciliophora", "Dinoflagellata", 
    "Diatomea", "Chlorophyta", "Protalveolata", "Cercozoa", 
    "Holozoa", "Labyrinthulomycetes",
    "Others", "Unclassified"))

# Arrange order of samples based on clustering
taxa_bplot_data$samplename <- factor(taxa_bplot_data$samplename,
  levels = c("DH12", "DH15", "DH14", "DH13", "DH16",
    "DH3", "DH9", "DH11", "DH10", "DH2",
    "DH8", "DH6"))
```

## Plot taxonomy barplot
You can plot a percent stacked taxonomy bar graph using the code below:

```r
# Choose color palette
cust_color <- c(brewer.pal(n = 8, name = 'Set1'),
         "#ABABAB", "#000000")

ggplot(taxa_bplot_data,
  aes(fill = Phylum, y = count, x = samplename)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  scale_fill_manual(values = cust_color) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(0.6, 'cm'))
```

<div class="figure" style="text-align: center">
<img src="dh_markdown_v2_files/figure-html/plot taxonomy barplot-1.png" alt="**Figure 5.3.1.** Heatmap of log bias-adjusted abundance of identified phyla."  />
<p class="caption">**Figure 5.3.1.** Heatmap of log bias-adjusted abundance of identified phyla.</p>
</div>

You can add information on the plot by manually editing it.
