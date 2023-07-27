tab_evo |>  
  ggplot()+
  geom_point(aes(y=start_d0,x=start_Kslope, col=edge_pre/120))+
  geom_abline(intercept = 0.25, slope = -1) + 
  geom_hline(yintercept= 0.25) + 
  facet_grid(rows=vars(paste0("*r<sub>0</sub>* = log(",Fec,")")),cols=vars(paste("*m* =",Mort)))+
  scale_x_continuous("initial dispersal density-dependence (*&Delta;<sub>K-0</sub>*)")+ scale_color_continuous(type="viridis")+
  scale_y_continuous("initial *d<sub>0</sub>*")+
  ggtitle(label="",
          subtitle="A) replicates with no individual variation in dispersal traits")+
  theme_bw()+
  theme(legend.position = "none",
        strip.text.x = element_markdown(),
        strip.text.y = element_markdown(angle = 0),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown(),
        plot.subtitle = element_markdown()
  )

mod_speed_evoB <- brm(bf(edge_post|trials(120)~ (start_d0 + start_Kslope) * Mort * Fec  + 
                         (1|scenario)),
                    family=binomial,
                    data=tab_evo,backend="cmdstanr",
                    seed=42,
                    prior=c(
                      set_prior("normal(0,1.5)",class="Intercept"),
                      set_prior("normal(0,1)",class="b"),
                      set_prior("normal(0,1)",class="sd")
                    )
)


fits_evo_speed=select(tab_evo,Mort,Fec) |> distinct() |> 
  expand_grid(start_Kslope=c(-5:5)/10, start_d0=c(0:5)/10) |> 
  add_epred_draws(mod_speed_evoB,re_formula=NA) |> 
  mean_hdi()

ggplot(fits_evo_speed)+
  geom_contour_filled(aes(x=start_d0,y=start_Kslope, z = .epred/120),
                      bins=10)+
  geom_contour(aes(x=start_d0,y=start_Kslope, z = .epred/120),
               bins=10,col="black")+
  geom_polygon(data=tibble(y=c(-0.5,-0.5,0),x=c(0.5,0,0)),aes(x=x,y=y),fill="white")+
  geom_polygon(data=tibble(y=c(0.5,0.5,0),x=c(0,0.5,0.5)),aes(x=x,y=y),fill="white")+
  geom_point(data=tab_evo,aes(x=start_d0,y=start_Kslope),col="grey",alpha=0.2)+
  facet_grid(rows=vars(paste0("*r<sub>0</sub>* = log(",Fec,")")),cols=vars(paste("*m* =",Mort)))+
  scale_y_continuous("initial dispersal density-dependence (*&Delta;<sub>K-0</sub>*)")+ scale_color_continuous(type="viridis")+
  scale_x_continuous("initial *d<sub>0</sub>*")+
  ggtitle(label="",
          subtitle="xxxx")+
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






