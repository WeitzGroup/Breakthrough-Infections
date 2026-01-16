# on first run, install ggmagnify for easy inset
# install.packages(
#   "ggmagnify",
#   repos = c("https://hughjonesd.r-universe.dev", "https://cloud.r-project.org")
# )

require(deSolve)
require(magrittr)
require(tidyverse)
require(cowplot)
require(usmap)
require(latex2exp)
library(ggrepel)
library(scales)
library(viridis)
library(ggmagnify)

# read in simulation output for different disease parameters
diffdis_df <- read.csv("diffdisease-ode-output.csv") 

# supplemental figure - fV with different disease parameters
pdf("figs/supp-fv.pdf", height=10, width=10)


diffdis_df %>%
  data.frame() %>%
  rename(phi = assortativity) %>%
  # exclude small/NA values (indicative of elimination)
  filter(!is.nan(IV) & !is.nan(IU)) %>%
  filter(IU + IV >= 1) %>%
ggplot(aes(x=coverage, y=fV, color=as.factor(phi))) +
  geom_line()+
  scale_x_continuous(limits=range(c(0,1)), expand=c(0,0))+
  theme_classic(base_size=20)+
  theme(legend.position="bottom",
        legend.title.align = 0.5)+
  ylab(expression(paste("Breakthrough Fraction (", italic(f[V]), ")")))+
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  guides(color = guide_legend(expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"))+
  scale_color_viridis_d()+
  # dynamic labelling of facets based on variable values
  facet_grid(rows = vars(vaccine_failure),
             cols = vars(R0),
             labeller = label_bquote(rows = epsilon ==.(vaccine_failure),
                                     cols = R[0] == .(R0)))+
  theme(panel.spacing.x = unit(3, "lines"))
dev.off()


# supplemental figure - fV with different disease parameters
pdf("figs/supp-Iv.pdf", height=10, width=10)
diffdis_df %>%
  data.frame() %>%
  rename(phi = assortativity) %>%
  ggplot(aes(x=coverage, y=incidence_V, color=as.factor(phi))) +
  geom_line()+
  # set limits to avoid space below 0 infections
  scale_x_continuous(limits=range(c(0,1)), expand=c(0,0))+
  coord_cartesian(ylim = c(0, NA), expand=FALSE)+ 
  theme_classic(base_size=20)+
  theme(legend.position="bottom",
        legend.title.align = 0.5)+
  ylab("Breakthrough Infections (annual per 100K)")+
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  guides(color = guide_legend(expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"))+
  # dynamic labelling of facets based on variable values
  facet_grid(rows = vars(vaccine_failure),
             cols = vars(R0),
             labeller = label_bquote(rows = epsilon ==.(vaccine_failure),
                                     cols = R[0] == .(R0)))+
  scale_color_viridis_d()+
  theme(panel.spacing.x = unit(3, "lines"))
  
dev.off()

# supplemental figure - herd immunity threshold with mild assortativity

# read in and combine two dataframes here (inset gives smaller steps in coverage at high values)
lowphi_df <- read.csv("lowphi-ode-output.csv") %>% data.frame()
lowphi_df <- read.csv("lowphi-inset-ode-output.csv") %>% 
              data.frame() %>%
              rbind(lowphi_df)


# plot larger figure with coverage from 0 to 1
lowphi_main <- lowphi_df %>%
  data.frame() %>%
  ggplot(aes(x=coverage, y=incidence_V, color=as.factor(phi))) +
  geom_line()+
  theme_classic(base_size=20)+
  theme(legend.position="bottom",
        legend.title.align = 0.5)+
  ylab("Breakthrough Infections (annual per 100K)")+
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  guides(color = guide_legend(expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"))+
  coord_cartesian(ylim = c(0, max(lowphi_df$incidence_V)+1), xlim=c(0, 1.1), expand=FALSE)+ 
  scale_color_viridis_d()


pdf("figs/supp-lowphiHIT.pdf", height=8, width=8)

# use geo_magnify to add an inset
lowphi_main +
  geom_magnify(from=c(xmin=.95, xmax=1, ymin=0, ymax=10),
               to=c(xmin=.3, xmax=.85, ymin=1, ymax=8),
               axes="x")
dev.off()

# solve for phat across different assortativity values
hatp_df <- read.csv("ode_out.csv") %>%
  data.frame() %>%
  rename(p = coverage,
         annual_CU = incidence_U,
         annual_CV = incidence_V,
         phi=assortativity)  %>% 
  group_by(phi) %>%
  slice_max(annual_CV) %>%
  rename(p_hat = p)


# read in csv of county-level MMR vaccine coverage
pdf("figs/county_map.pdf", height=10, width=8)

map1 <- read.csv("mmr_data_us_counties.csv") %>%
  rename(fips = FIPS) %>%
  # manually bin values based on estimated p_c and \hat{p}
  mutate(category = case_when(
    SY2022_23 <= .77 ~ "1",
    SY2022_23 > .77 & SY2022_23 <= .996 ~ "2",
    SY2022_23 > .996 ~ "3"
  )) %>%
  plot_usmap(data=., values="category", color="grey", size=.1)+
  theme_void(base_size=20)+
  theme(legend.position="bottom")+
  scale_fill_viridis_d(name="",
                      option="A",
                      # allow latex special characters in map
                      labels = c("3" = TeX("$p > p_{c}"), 
                                 "2" = TeX("$\\hat{p} < p \\leq \\p_{c}$"), 
                                 "1" = TeX("p \\leq \\hat{p}$")),
                      guide=guide_legend(reverse=FALSE),
                      na.value="white")

county_data <- read.csv("mmr_data_us_counties.csv")

# bottom panel, histogram
map2 <- ggplot()+
  geom_histogram(data=county_data, aes(x=SY2022_23), fill="grey70", closed="left")+
  ylab("Count (Counties)")+
  xlab("Vaccine Coverage (p)")+
  # indicate values for 
  geom_segment(aes(x=.77, y=0, yend=470))+
  geom_segment(aes(x=.997, y=0, yend=470))+
  geom_segment(data=hatp_df %>% filter(phi %in% c(.3, .6)), 
               aes(x=p_hat, y=0, yend=470),
               linetype="dotted")+
  scale_x_continuous(expand = expansion(mult = c(0, .1)),
                     breaks=seq(from=0, to=1, by=.1))+
  scale_y_continuous(limits=c(0, 650),expand = expansion(mult = c(0, 0)))+
  theme_classic(base_size=20)+
  annotate("text", x=.996, y=550, size=8, label="p[c]", parse=TRUE)+
  annotate("text", x=.77, y=550, size=8, label="hat(p)", parse=TRUE)+
  theme(axis.text.x = element_text(hjust = 0),
        axis.ticks.x  = element_blank())

plot_grid(map1, map2, rel_heights=c(6, 2), nrow =2)
dev.off()


# get values for number of counties above/below threshold value for coverage
read.csv("mmr_data_us_counties.csv") %>% 
  filter(SY2022_23 >= .997) %>%
  nrow()

read.csv("mmr_data_us_counties.csv") %>% 
  filter(SY2022_23 < .997 & SY2022_23 > .77) %>%
  nrow()
