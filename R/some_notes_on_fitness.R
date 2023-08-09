K = 250
Fec = c(2, 5)
Mort = 0.5

# fitness of an individual staying in a K/2 patch

Fres = exp(log(Fec)*(1-0.9/1))

# fitness of a individual deciding to disperse to an empty patch
Fdisp = exp(log(Fec*(1-Mort))*(1-0/1))

Fdisp/Fres
