clear all

x0 = 0;

par = -1;
t = 1;
rsize = 10;
sigma = 3;
damping = 0.0000001;
vmax_c = .5;

tic
par0 = -2;
for ii = 1:40
    par = par0 + ii/10;
    C(ii,:) = SPSO_Gen(@Pitch_Bif,x0,par,rsize,sigma,damping,vmax_c);
    xax = par*ones(1,length(C(ii,:)));
    plot(xax,C(ii,:),'b.')
    title('SPSO: Subcritical Pitchform Normal Form')
    xlabel('\mu')
    ylabel('x')
    hold on

end
toc

