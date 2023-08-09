ggplot(tab_evo |> 
         filter(scenario %in% 
                  c(
                    "F=5_M=0.5_h2=1_slope=6_midpoint=1_VP.DDD=1_VP_dmax=0",
                    "F=5_M=0.5_h2=1_slope=4_midpoint=1_VP.DDD=1_VP_dmax=0",
                    "F=5_M=0.5_h2=1_slope=2_midpoint=1_VP.DDD=1_VP_dmax=0"
                    )
                )|> 
         expand_grid(N=c(0:100)/100))+
  geom_line(aes(N,0.5/(1+exp(-start_slope*(N-start_midpoint))),group=scenario))+
  geom_line(aes(N,mean_dmax/(1+exp(-mean_slope*(N-mean_midpoint))),group=replicateID),col="red")


data |>  
  filter(ticks ==121 & pxcor < 5) |> 
  select(seedID,replicateID, scenario, 
         dispersal_mortality, heritability, fecundity,
         K,
         is.VP_dmax, is.VP_DDD,
         start_slope, start_midpoint,
         start_d0, start_dK,        
         start_slopeA_0_K, N_predispersal) |>
  ggplot()+
  geom_boxplot(aes(factor(fecundity),N_predispersal))
