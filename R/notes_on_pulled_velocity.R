
vf = function(fec,d){ 
  # this will plateau a few % above 1 rather than 1 (optim errors?)
  # and will also give - values of vf if net fecundity < 1 (when vf should be 0)
  # so some post-proc needed
  test = optimise(
  function(K,rho=fec,m=d){
    log(rho*(1+ m *(cosh(K)-1)))/K
    },
    lower=0,upper=500 # upper cannot be too high or else get Inf/NA warnings
)

return(test$objective)
}


demo=expand_grid(
  fec = c(1.5,5),
  d = c(1:50)/100,
  m = c(0.1,0.5,0.9), 
  vf=NA
) |> 
  mutate(net_fec = fec*(1-d*m), # remove the individuals lost to disp mortality from fec
         net_disp = d*(1-m))    # remove the individuals lost to disp mortality


for(i in 1:dim(demo)[1]){
  demo$vf[i] = vf(demo$net_fec[i],demo$net_disp[i])
}

demo |> 
  mutate(vf_corrected = case_when(vf>1~1,
                                  vf<0~0,
                                  T~vf)) |> 
  ggplot()+geom_point(aes(x=d,y=vf_corrected,col=net_fec>1))
