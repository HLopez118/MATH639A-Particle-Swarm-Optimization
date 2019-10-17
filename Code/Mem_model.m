function dV = Mem_model(t,v,par)
vl = -60;
vca = 120;
i = par;
gl = 2;
gca = 4;
c = 20;
v1 = -1.2;
v2 = 18;

minf = .5*(1+tanh((v-v1)/v2));

dV = (i+gl*(vl-v)+gca*minf*(vca-v))/c;


end