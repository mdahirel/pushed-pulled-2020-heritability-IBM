alpha <- c(-6,-4,-2,0,2,4,6)
beta <- c(0,0.25,0.5,0.75,1)
dmax <- 0.5

parsamples <- expand_grid(alpha,beta,dmax) |> 
  expand_grid(
    tibble(zscores_a=rnorm(1000,0,1),
           zscores_b=rnorm(1000,0,1),
           zscores_dmax=rnorm(1000,0,1))) |> 
  mutate(alpha_sample = alpha + (zscores_a*6),
         beta_sample = beta + (zscores_b*0.5),
         dmax_sample=plogis(zscores_dmax*1.5))

## version with phenotypic var in DDD

parsamples |> 
  mutate(
    d0 = dmax/(1+exp(-alpha * (0 - beta))),
    dK = dmax/(1+exp(-alpha * (1 - beta))),
    d0_sample = dmax/(1+exp(-alpha_sample * (0 - beta_sample))),
    dK_sample = dmax/(1+exp(-alpha_sample * (1 - beta_sample)))
  ) |> 
  mutate(delta_sample = dK_sample - d0_sample,
         delta = dK - d0) |> 
  group_by(alpha,beta,dmax) |> 
  summarise(delta=mean(delta), meandelta_sample = mean(delta_sample),
            d0=mean(d0),meand0_sample = mean(d0_sample)) |> 
  ggplot() + 
  geom_contour_filled(aes(x=alpha,y=beta,z=meandelta_sample))


parsamples |> 
  mutate(
    d0 = dmax/(1+exp(-alpha * (0 - beta))),
    dK = dmax/(1+exp(-alpha * (1 - beta))),
    d0_sample = dmax/(1+exp(-alpha_sample * (0 - beta_sample))),
    dK_sample = dmax/(1+exp(-alpha_sample * (1 - beta_sample)))
  ) |> 
  mutate(delta_sample = dK_sample - d0_sample,
         delta = dK - d0) |> 
  mutate(i=1:35000) |> 
  group_by(i) |> 
  nest() |> 
  mutate(vstar = map(
    .x=data,
    .f = \(x) vstar(
      fec=1.5,
      disp=x$d0,
      mort=0)
  )
  )|> 
  mutate(vstar_sample = map(
    .x=data,
    .f = \(x) vstar(
      fec=1.5,
      disp=x$d0_sample,
      mort=0)
  )
  )|> 
  unnest(data) |> 
  unnest(vstar) |> 
  unnest(vstar_sample) |> 
  group_by(alpha,beta,dmax) |> 
  summarise(vstar=mean(vstar), meanvstar_sample = median(vstar_sample)) |>  
  ggplot()+geom_point(aes(vstar,meanvstar_sample))+
  geom_abline(intercept=0,slope=1)
  
