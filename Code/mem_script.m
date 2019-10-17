clear all

x0 = -10;

par = -100;
t = 1;
rsize = 10;
sigma = 200;
damping = 0.0000001;
vmax_c = .5;

tic

for ii = -215:-215
    par = ii;
    C(ii+216,:) = SPSO_Gen(@Mem_model,x0,par,rsize,sigma,damping,vmax_c);
    xax = par*ones(1,length(C(ii+216,:)));
    plot(xax,C(ii+216,:),'b.')
    hold on

end
toc
