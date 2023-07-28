alpha <- c(-6:6)
beta <- c(-1,0,0.25,0.5,0.75,1,2)

expand_grid(alpha,beta) |> 
  mutate(d0 = 0.5/(1+exp(-alpha * (0 - beta))),
         dK = 0.5/(1+exp(-alpha * (1 - beta)))
  ) |> 
  ggplot() + 
  geom_polygon(data=tibble(y=c(-0.5,-0.5,0),x=c(0.5,0,0)),aes(x=x,y=y),fill="white")+
  geom_polygon(data=tibble(y=c(0.5,0.5,0),x=c(0,0.5,0.5)),aes(x=x,y=y),fill="white")+
  geom_point(aes(d0,dK-d0,col=beta),size=3)+
  geom_line(aes(d0,dK-d0,group=alpha))
