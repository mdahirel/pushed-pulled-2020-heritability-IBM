### KEY QUESTION: REINCLUDE THE 0 EXPANSION DATAPOINTS FOR SPEED?
### IF YES: would even more strongly justify putting as suppl



mod_speed_fix <- brm(bf(edge_post|trials(120)~ (start_d0 * start_Kslope) * Mort * Fec  + 
                           (1|scenario)),
                      family=binomial,
                      data=tab_fix,backend="cmdstanr",
                      seed=42,
                      prior = prior_binomial
)


mod_speed_noevo <- brm(bf(edge_post|trials(120)~ (start_d0 * start_Kslope) * Mort * Fec  + 
                           (1|scenario)),
                      family=binomial,
                      data=tab_noevo,backend="cmdstanr",
                      seed=42,
                      prior = prior_binomial
)

mod_speed_evo <- brm(bf(edge_post|trials(120)~ (start_d0 * start_Kslope) * Mort * Fec  + 
                         (1|scenario)),
                    family=binomial,
                    data=tab_evo,backend="cmdstanr",
                    seed=42,
                    prior = prior_binomial
)

fits_fix_speed=select(tab_fix,Mort,Fec) |> distinct() |> 
  expand_grid(start_Kslope=c(-5:5)/10, start_d0=c(0:5)/10) |> 
  add_epred_draws(mod_speed_fix,re_formula=NA) |> 
  mean_hdi()

ggplot(fits_fix_speed)+
  geom_contour_filled(aes(x=start_d0,y=start_Kslope, z = .epred/120),
                      breaks=c(0:10)/10)+
  geom_contour(aes(x=start_d0,y=start_Kslope, z = .epred/120),
               breaks=c(0:10)/10,col="black")+
  geom_polygon(data=tibble(y=c(-0.5,-0.5,0),x=c(0.5,0,0)),aes(x=x,y=y),fill="white")+
  geom_polygon(data=tibble(y=c(0.5,0.5,0),x=c(0,0.5,0.5)),aes(x=x,y=y),fill="white")+
  #geom_point(data=tab_fix,aes(x=start_d0,y=start_Kslope),col="grey",alpha=0.2)+
  facet_grid(rows=vars(paste0("*r<sub>0</sub>* = log(",Fec,")")),cols=vars(paste("*m* =",Mort)))+
  scale_y_continuous("initial dispersal density-dependence (*&Delta;<sub>K-0</sub>*)")+ scale_color_continuous(type="viridis")+
  scale_x_continuous("initial *d<sub>0</sub>*")+
  ggtitle(label="",
          subtitle="A) replicates with no individual variation in dispersal traits")+
  labs(fill="predicted speed\n(patches/generation)")+
  theme_bw()+
  theme(
    strip.text.x = element_markdown(),
    strip.text.y = element_markdown(angle = 0),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    plot.subtitle = element_markdown(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )




fits_noevo_speed=select(tab_noevo,Mort,Fec) |> distinct() |> 
  expand_grid(start_Kslope=c(-5:5)/10, start_d0=c(0:5)/10) |> 
  add_epred_draws(mod_speed_noevo,re_formula=NA) |> 
  mean_hdi()

ggplot(fits_noevo_speed)+
  geom_contour_filled(aes(x=start_d0,y=start_Kslope, z = .epred/120),
                      breaks=c(0:10)/10)+
  geom_contour(aes(x=start_d0,y=start_Kslope, z = .epred/120),
               breaks=c(0:10)/10,col="black")+
  geom_polygon(data=tibble(y=c(-0.5,-0.5,0),x=c(0.5,0,0)),aes(x=x,y=y),fill="white")+
  geom_polygon(data=tibble(y=c(0.5,0.5,0),x=c(0,0.5,0.5)),aes(x=x,y=y),fill="white")+
  #geom_point(data=tab_noevo,aes(x=start_d0,y=start_Kslope),col="grey",alpha=0.2)+
  facet_grid(rows=vars(paste0("*r<sub>0</sub>* = log(",Fec,")")),cols=vars(paste("*m* =",Mort)))+
  scale_y_continuous("initial dispersal density-dependence (*&Delta;<sub>K-0</sub>*)")+ scale_color_continuous(type="viridis")+
  scale_x_continuous("initial *d<sub>0</sub>*")+
  ggtitle(label="",
          subtitle="B) individual variation in dispersal traits but no evolution (*h<sup>2</sup>* = 0)")+
  labs(fill="predicted speed\n(patches/generation)")+
  theme_bw()+
  theme(
    strip.text.x = element_markdown(),
    strip.text.y = element_markdown(angle = 0),
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    plot.subtitle = element_markdown(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )


fits_evo_speed=select(tab_evo,Mort,Fec) |> distinct() |> 
  expand_grid(start_Kslope=c(-5:5)/10, start_d0=c(0:5)/10) |> 
  add_epred_draws(mod_speed_evo,re_formula=NA) |> 
  mean_hdi()

ggplot(fits_evo_speed)+
  geom_contour_filled(aes(x=start_d0,y=start_Kslope, z = .epred/120),
                      breaks=c(0:10)/10)+
  geom_contour(aes(x=start_d0,y=start_Kslope, z = .epred/120),
               breaks=c(0:10)/10,col="black")+
  geom_polygon(data=tibble(y=c(-0.5,-0.5,0),x=c(0.5,0,0)),aes(x=x,y=y),fill="white")+
  geom_polygon(data=tibble(y=c(0.5,0.5,0),x=c(0,0.5,0.5)),aes(x=x,y=y),fill="white")+
  #geom_point(data=tab_evo,aes(x=start_d0,y=start_Kslope),col="grey",alpha=0.2)+
  facet_grid(rows=vars(paste0("*r<sub>0</sub>* = log(",Fec,")")),cols=vars(paste("*m* =",Mort)))+
  scale_y_continuous("initial dispersal density-dependence (*&Delta;<sub>K-0</sub>*)")+ scale_color_continuous(type="viridis")+
  scale_x_continuous("initial *d<sub>0</sub>*")+
  ggtitle(label="",
          subtitle="C) individual variation in dispersal traits *and* evolution possible (*h<sup>2</sup>* = 1)")+
  labs(fill="predicted speed\n(patches/generation)")+
  theme_bw()+
  theme(
        strip.text.x = element_markdown(),
        strip.text.y = element_markdown(angle = 0),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        plot.subtitle = element_markdown(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  )






