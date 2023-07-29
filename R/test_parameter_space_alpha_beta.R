alpha <- c(-6,-4,-2,0,2,4,6)
beta <- c(0,0.25,0.5,0.75,1)
dmax <- 0.5

expand_grid(alpha,beta,dmax) |> 
  mutate(d0 = dmax/(1+exp(-alpha * (0 - beta))),
         dK = dmax/(1+exp(-alpha * (1 - beta)))
  ) |> 
  ggplot() + 
  geom_contour_filled(aes(x=alpha,y=beta,z=(dK-d0)))






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
  mutate(d0 = dmax/(1+exp(-alpha_sample * (0 - beta_sample))),
         dK = dmax/(1+exp(-alpha_sample * (1 - beta_sample)))
  ) |> 
  mutate(delta = dK-d0) |> 
  group_by(alpha,beta,dmax) |> 
  summarise(meandelta=mean(delta)) |> 
  ggplot() + 
  geom_contour_filled(aes(x=alpha,y=beta,z=meandelta))

## version with all phenotypic var

parsamples |> 
  mutate(d0 = dmax_sample/(1+exp(-alpha_sample * (0 - beta_sample))),
         dK = dmax_sample/(1+exp(-alpha_sample * (1 - beta_sample)))
  ) |> 
  mutate(delta = dK-d0) |> 
  group_by(alpha,beta,dmax) |> 
  summarise(meandelta=mean(delta)) |> 
  ggplot() + 
  geom_contour_filled(aes(x=alpha,y=beta,z=meandelta))

## version with only phenotypic var in dmax

parsamples |> 
  mutate(d0 = dmax_sample/(1+exp(-alpha * (0 - beta))),
         dK = dmax_sample/(1+exp(-alpha * (1 - beta)))
  ) |> 
  mutate(delta = dK-d0) |> 
  group_by(alpha,beta,dmax) |> 
  summarise(meandelta=mean(delta)) |> 
  ggplot() + 
  geom_contour_filled(aes(x=alpha,y=beta,z=meandelta))




## version with all phenotypic var

parsamples |> 
  mutate(d0 = dmax_sample/(1+exp(-alpha_sample * (0 - beta_sample))),
         dK = dmax_sample/(1+exp(-alpha_sample * (1 - beta_sample)))
  ) |> 
  mutate(delta = dK-d0) |> 
  group_by(alpha,beta,dmax) |> 
  summarise(meandelta=mean(delta)) |> 
  mutate(fix_delta = 0.5/(1+exp(-alpha * (1 - beta))) - 
           0.5/(1+exp(-alpha * (0 - beta)))) |> 
  ggplot() + 
  geom_point(aes(x=fix_delta,y=meandelta))+
  lims(x=c(-0.5,0.5),y=c(-0.5,0.5))

#including var in alpha and beta reduces the range of var in DDD, through weird averaging stuff


