
vf = function(fec,d){ 
  # this will plateau a few % above 1 rather than 1 (optim errors?)
  # and will also give - values of vf if net fecundity < 1 (when vf should be 0)
  # so some post-proc needed
  test = optimise(
  function(K,rho=fec,m=d){
    log(rho*(1+ m *(cosh(K)-1)))/K
    },
    lower=0,upper=707 #cosh(K>709) = overflow, estimated as Inf; upper a few orders lower by safety
)

return(test$objective)
}


demo=expand_grid(
  fec = c(1:10)/2,
  d = c(0:50)/100, 
  vf=NA
) 

for(i in 1:dim(demo)[1]){
  demo$vf[i] = vf(demo$fec[i],demo$d[i])
}

demo=demo |> 
  mutate(vf=round(vf,4)) |> 
  mutate(vf_corrected = case_when(vf>1~1,
                                  vf<0~0,
                                  T~vf))
          

demo |>
  ggplot()+geom_contour_filled(aes(x=(d),y=fec,z=(vf_corrected)))




mod=brm(edge_pre|rate(vf_corrected*120)~poly(start_Kslope,2)*Mort*Fec + (1|scenario),
        data=filter(tab_fix,vf_corrected>0),
        family=poisson,
        backend="cmdstanr")

tab_fix |> select(Mort,Fec,start_Kslope) |> distinct() |> mutate(vf_corrected=1/120) |> 
  add_epred_draws(mod,re_formula=NA) |> 
  ggplot()+stat_lineribbon(aes(start_Kslope,.epred))+
  geom_point(data=tab_fix,aes(start_Kslope,edge_pre/(vf_corrected*120)))+facet_wrap(~Mort+Fec)



tab_fix  |> 
  ggplot()+
  geom_point(aes(vf_corrected,edge_pre/120,fill=start_Kslope),pch=21,size=2.5)+
  geom_abline(intercept = 0,slope=1)+
  scale_fill_gradient2("&Delta;<sub>K-0</sub>",high="#f1a340",mid="#f7f7f7",low="#998ec3",midpoint = 0)+
  scale_x_continuous("Theoretical velocity if density-independent dispersal (i.e. if &Delta;<sub>K-0</sub> = 0)")+
  scale_y_continuous("Observed velocity")+
  theme_bw()+
  theme(legend.title=element_markdown(),axis.title.x = element_markdown())

tab_noevo  |> 
  ggplot()+
  geom_point(aes(vf,edge_pre/120,fill=start_Kslope),pch=21,size=2.5)+
  geom_abline(intercept = 0,slope=1)+
  scale_fill_gradient2("&Delta;<sub>K-0</sub>",high="#e66101",mid="#f7f7f7",low="#5e3c99",midpoint = 0)+
  scale_x_continuous("Theoretical velocity if density-independent dispersal (i.e. if &Delta;<sub>K-0</sub> = 0)")+
  scale_y_continuous("Observed velocity")+
  theme_bw()+
  theme(legend.title=element_markdown(),axis.title.x = element_markdown())+
  coord_cartesian(xlim=c(0,1))

tab_evo  |> filter(is.VP_dmax==0) |> 
  ggplot()+
  geom_point(aes(vf,edge_pre/120,fill=start_Kslope),pch=21,size=2)+
  geom_abline(intercept = 0,slope=1)+
  scale_fill_gradient2("&Delta;<sub>K-0</sub>",high="#f1a340",mid="#f7f7f7",low="#998ec3",midpoint = 0)+
  scale_x_continuous("Theoretical velocity if density-independent dispersal (i.e. if &Delta;<sub>K-0</sub> = 0)")+
  scale_y_continuous("Observed velocity")+
  theme_bw()+
  theme(legend.title=element_markdown(),axis.title.x = element_markdown())+
  coord_cartesian(xlim=c(0,1))
