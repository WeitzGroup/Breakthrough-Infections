require(deSolve)
require(magrittr)
require(tidyverse)
require(cowplot)
require(usmap)
require(latex2exp)
library(ggrepel)
library(viridis)
library(scales)

# helper function to add legend to multiplot figures
get_legend_35 <- function(plot) {
  # return all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  nonzero <- vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE)
  idx <- which(nonzero)
  # return first non-zero legend if exists, and otherwise first element (which will be a zeroGrob) 
  if (length(idx) > 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}



# read in and transform the outbreak data
outbreak_df <- data.frame(State = state.name, Abb = state.abb) %>%
  right_join(read.csv("recent-outbreaks.csv")) %>%
  # if unknown and unvaccinated are not reported separately, assume half of combined category are unvaccinated
  mutate(Unknown = ifelse(is.na(Unknown), .5*Unknown.unvaccinated, Unknown),
         Unvaccinated = ifelse(is.na(Unvaccinated), .5*Unknown.unvaccinated, Unvaccinated)) %>%
  mutate(fv_mid = At.least.one.dose/(At.least.one.dose+Unvaccinated)) %>%
  # uncertainty: add upper and lower bound based on assumption of vaccination rate in known vs unknown populations
  mutate(fv_lower = (At.least.one.dose+fv_mid/3*Unknown)/(Unknown.unvaccinated+At.least.one.dose),
         fv_upper = (At.least.one.dose+fv_mid*3*Unknown)/(Unknown.unvaccinated+At.least.one.dose)) %>%
  rename(Coverage = State.Coverage) %>%
  mutate(Name = paste0(County, ", ", Abb))

# read in simulation output
df <- read.csv("ode_out.csv") %>%
  data.frame() %>%
  rename(p = coverage,
         annual_CU = incidence_U,
         annual_CV = incidence_V,
         phi = assortativity) 

# solve for HIT values (will use to draw caps on p vs f_V plot)
df_HIT <- df %>%
  filter(IU+IV >= 1) %>%
  group_by(phi) %>%
  slice_max(p) %>%
  rename(pHIT = p) %>%
  rename(fvHIT = fV) 

# interpolate the value of p for observed breakthrough fractions 
outbreak_df <- outbreak_df %>%
                mutate(infer_p = approx(df_phi0$fV, df_phi0$p, xout = outbreak_df$fv_mid, method = "linear")$y,
                       infer_p_lower = approx(df_phi0$fV, df_phi0$p, xout = outbreak_df$fv_lower, method = "linear")$y,
                       infer_p_upper = approx(df_phi0$fV, df_phi0$p, xout = outbreak_df$fv_upper, method = "linear")$y)


#plot figure 1
pdf("figs/fig1.pdf", height=6, width=12)
# panel A
p1a <- ggplot()+
  geom_line(data=df %>% filter(phi==0), aes(x=p, y=fV), color="#440154")+
  geom_point(data=outbreak_df, aes(x=infer_p, y=fv_mid, shape=Abb), size=3)+
  # uncertainty
  geom_linerange(data=outbreak_df, aes(x=infer_p, ymin=-Inf, ymax=fv_mid), alpha=.8, lty="dashed")+
  # state labels
  geom_text(data=outbreak_df %>% filter(Abb %in% c("ND", "SC")), aes(x=infer_p, y=fv_mid, label=Abb), nudge_y=.015, nudge_x=.02)+
  geom_text(data=outbreak_df %>% filter(Abb=="AZ"), aes(x=infer_p, y=fv_mid, label=Abb), nudge_y=.015, nudge_x=-.02)+
  geom_text(data=outbreak_df %>% filter(infer_p >= .55), aes(x=infer_p, y=fv_mid, label=Abb), nudge_y=.015, nudge_x=-.03)+
  # set limits, labels, and theme
  scale_x_continuous(limits=range(c(-.02,1.1)), expand=c(0,0))+
  scale_y_continuous(limits=range(c(-.01,.43)), expand=c(0,0))+
  theme_classic(base_size=20,)+
  theme(legend.position="none",
        plot.margin = margin(t = 20, r = 5, b = 10, l = 5, unit = "pt"))+
  ylab(expression(paste("Breakthrough fraction (", italic(f[V]), ")")))+
  xlab(expression(paste("Model-predicted vaccine coverage (", italic(p),")"))) +
  # manually set shapes for different states
  scale_shape_manual(values=c("MI"=15, "ND"=19,"UT"=17, "TX"=18, "NM"=8, 
                              "AZ"=7, "SC"=18))
# panel B
p1b <- ggplot()+
  geom_point(data=outbreak_df, aes(x=infer_p, y=Coverage, shape=Abb), size=3)+
  geom_segment(data=outbreak_df, aes(x=infer_p_lower, xend=infer_p_upper, y=Coverage))+
  # state labels
  geom_text(data=outbreak_df %>% filter(!Abb%in%c("ND", "AZ")), aes(x=infer_p, y=Coverage, label=Abb), nudge_y=.03)+
  geom_text(data=outbreak_df %>% filter(Abb=="ND"), aes(x=infer_p, y=Coverage, label=Abb), nudge_y=.03, nudge_x=.01)+
  geom_text(data=outbreak_df %>% filter(Abb=="AZ"), aes(x=infer_p, y=Coverage, label=Abb), nudge_y=-.03)+
  # x = y line
  geom_abline(linetype="dashed", color="grey30")+
  # set limits, labels, and theme
  theme_classic(base_size=20)+
  theme(legend.position="none",
        plot.margin = margin(t = 20, r = 5, b = 10, l = 25, unit = "pt"))+
  scale_x_continuous(limits=range(c(-.02,1.1)), expand=c(0,0))+
  scale_y_continuous(limits=range(c(-.02,1.1)), expand=c(0,0))+
  xlab("Model-predicted vaccine coverage")+
  ylab("Surveyed kindergarten vaccine coverage")+
  # manually set shapes for different states
  scale_shape_manual(values=c("MI"=15, "ND"=19,"UT"=17, "TX"=18, "NM"=8, 
                              "AZ"=7, "SC"=18))
  
# combine two panels, output
plot_grid(p1a, p1b, labels="AUTO", label_size=20, hjust=-1, label_x=-.03)
dev.off()

# Figure 2
pdf("figs/fig2.pdf", height=6, width=6)
ggplot(outbreak_df)+
  #predict relationship between p and fV by phi
  geom_line(data = df %>% filter(phi %in% c(0,.3, .6, .9, .98)), aes(x=p, y=fV, color=as.factor(phi)))+
  # caps at herd immunity
  geom_point(data = df_HIT %>% filter(phi %in% c(0,.3, .6, .9, .98)), aes(x=pHIT, y=fvHIT, color=as.factor(phi)))+
  # outbreak data as points
  geom_linerange(aes(x=Coverage, ymin=fv_lower, ymax=fv_upper))+
  geom_point(aes(x=Coverage, y=fv_mid, shape=Abb), size=3)+
  # state labels
  geom_text(data=outbreak_df %>% filter(Abb %in% c("TX", "NM")), aes(x=Coverage, y=fv_mid, label=Abb), hjust=0, nudge_x=.02)+
  geom_text(data=outbreak_df %>% filter(!Abb %in% c("TX", "NM", "SC")), aes(x=Coverage, y=fv_mid, label=Abb), hjust=0, nudge_x = -.06, nudge_y=.006)+
  geom_text(data=outbreak_df %>% filter(Abb %in% c("SC")), aes(x=Coverage, y=fv_mid, label=Abb), hjust=0, nudge_x=.01, nudge_y=.006)+
  # limits, labels, and theme
  scale_x_continuous(limits=range(c(0,1.1)), expand=c(0,0))+
  scale_y_continuous(limits=range(c(-.02,.43)), expand=c(0,0))+
  scale_color_viridis_d()+
  theme_classic(base_size=20)+
  guides(color = guide_legend(override.aes=list(shape=NA),
                              expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"),
         shape="none")+
  theme(legend.position = c(0.05, 0.9), 
        legend.justification = c("left", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, hjust=0))+
  ylab(expression(paste("Breakthrough Fraction (", italic(f[V]), ")")))+
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) 
  # manually set shapes for different states
  scale_shape_manual(values=c("MI"=15, "ND"=19,"UT"=17, "TX"=18, "NM"=8, 
                              "AZ"=7, "SC"=18))
dev.off()

# identify \hat{p} to indicate as points on fig 4
df_peak <- df %>% 
            group_by(phi) %>%
            slice_max(annual_CV) 

# figure 4
pdf("figs/fig4.pdf", height=7, width=12)
# panel A: log-scale cases in vaccinated & unvaccinated
p4a <- df %>%
  # allows us to plot cases with linetype indicating vaccination status
  pivot_longer(c(annual_CU, annual_CV), names_to="names", values_to="Cases") %>%
  mutate(`Vaccination Status` = ifelse(names=="annual_CU", "Unvaccinated", "Vaccinated (Breakthrough)")) %>%
  ggplot()+
  geom_line(aes(x=p, y=Cases, color=as.factor(phi), linetype=`Vaccination Status`))+
  # limits, labels, theme, and pseudo-log transformation
  scale_y_continuous(trans = pseudo_log_trans(base = 10),
                     breaks=c(0, 10, 100,1000, 10000)) +
  theme_classic(base_size=20)+
  theme(legend.position="bottom",legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, hjust=.5),
        plot.margin = margin(t = 20, r = 5, b = 10, l = 25, unit = "pt"))+
  scale_color_viridis_d()+
  ylab(expression(paste("Infections (log-scaled, annual per 100K)"))) + 
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  scale_x_continuous(limits=range(c(0,1.1)), expand=c(0,0))+
  scale_linetype_manual(values=c("Vaccinated (Breakthrough)"="solid", "Unvaccinated"="dashed")) +
  # customize legend and add special characters
  guides(color = guide_legend(expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"),
         linetype = guide_legend(title.position="top"))

# panel B: linear-scaled, just breakthroughs
p4b <- ggplot()+
  geom_line(data = df , aes(x=p, y=annual_CV, color=as.factor(phi)))+
  # point to indicate \hat{p}
  geom_point(data = df_peak, aes(x=p, y=annual_CV, color=as.factor(phi)), size=3)+
  # legend, limits, theme
  theme_classic(base_size=20)+
  theme(legend.position="bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, hjust=.5),
        plot.margin = margin(t = 20, r = 5, b = 10, l = 25, unit = "pt"))+
  scale_color_viridis_d()+
  ylab(expression(paste("Infections (annual per 100K)"))) + 
  xlab(expression(paste("Vaccine coverage (", italic(p),")"))) +
  scale_x_continuous(limits=range(c(0,1.1)), expand=c(0,0))+
  # customize legend and add special characters
  guides(color = guide_legend(override.aes=list(shape=NA),
                              expression(paste("Assortativity (", italic(phi), ")")),
                              title.position = "top"))

# extract legend from panel A to then move to bottom of the combined plot
leg<-get_legend_35(p4a)
plot_grid(p4a+theme(legend.position="none"), 
          p4b+theme(legend.position="none"), 
          nrow=1,
          labels="AUTO") %>%
  plot_grid(., leg, nrow=2, rel_heights = c(9,1), label_size=20, hjust=-1, label_x=-.03)
dev.off()


# calculate p_HIT and p_hat for example values
eps <- .03
g <- .1
m<-1/(365*80)
b <- 15*(g+m)
r0 <- b/(g+m)

fv_eq_calc <- function(p, epsilon=eps){
  fv <- (p*epsilon)/(1-p*(1-epsilon))
  return(fv)
}

p_HIT <- 1/(1-eps)*(1-1/r0)
p_hat <- 1/(1-eps)*(1-1/sqrt(r0))
fv_eq_calc(.9)


