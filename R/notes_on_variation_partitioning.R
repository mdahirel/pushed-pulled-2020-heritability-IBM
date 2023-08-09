#partition the variation at the between-scenario level* into what's explained by fixed effects and what's not
#* so ignoring the "residual" within-scenario variation
#* 
varpart_mod=function(mod,sourcedata){
  
  tab <- subset(sourcedata, status == "expanded") |> 
    add_linpred_draws(mod,re_formula=NA) |> 
    group_by(.draw) |> 
    summarise(VF=var(.linpred))
  
  tab$VR <- VarCorr(mod,summary=FALSE)$scenario$sd[,"Intercept"]
  
  tab$ratio = tab$VF/(tab$VF+tab$VR)
  
  return(mean_hdi(tab$ratio))
  
}

varpart_mod(mod_gen_fixBB,tab_fix |> filter(status=="expanded"))
varpart_mod(mod_gen_noevoBB,tab_noevo |> filter(status=="expanded"))
varpart_mod(mod_gen_evoBB,tab_evo |> filter(status=="expanded"))
varpart_mod(mod_slope,tab_evo_DDD)
varpart_mod(mod_slope,tab_evo_dmax)