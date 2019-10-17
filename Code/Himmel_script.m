clear all

x0 = [0;0];

par = -1;
t = 1;
rsize = 20;
sigma = 7;
damping = 0.0000001;
vmax_c = .05;
% tic
C = SPSO_Gen(@Himmelblau,x0,par,rsize,sigma,damping,vmax_c);
% toc