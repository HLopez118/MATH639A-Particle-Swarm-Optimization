function C = PSO_Gen(func,x0,par,rsize,sigma,damping,vmax_c)

xlength = length(x0);
x = zeros(xlength,rsize^2);
vmax = vmax_c*(sigma/2);


for ii = 1:xlength
    x(ii,:) = sigma*(rand(1,rsize^2)-rand(1,rsize^2));
end

for tt = 1:rsize^2
    x(:,tt) =  x(:,tt) + x0;
end

c1 = 2;
c2 = 2;

crit = x;
pbest = crit;
velocity = zeros(xlength,rsize^2,2);

x3 = zeros(1,rsize^2);
for jj = 1:rsize^2
    x3(jj) = merit_fun(func, crit(:,jj),par);
end

itermax = 20000;

count = 1;
peval = zeros(1,rsize^2);
while mean(x3) > 1e-5 && count < itermax
    for kk = 1:rsize^2
        peval(kk) = merit_fun(func,pbest(:,kk),par);
    end
    
    [~,I] = min(peval);
    gbest = crit(:,I); %sets gbest
    for nn = 1:rsize^2
        
        omega = damping*10*rand; %damping parameter
        
        velocity(:,nn,1) = omega*velocity(:,nn,1) +...
            c1*(eye(xlength,xlength)*rand)*(pbest(:,nn)-crit(:,nn)) +...
            c2*(eye(xlength,xlength)*rand)*(gbest-crit(:,nn));
        for qq = 1:xlength %sets max velocity
            if velocity(qq,nn,1) > vmax
                velocity(qq,nn,1) = vmax;
            elseif velocity(qq,nn,1) < -vmax
                velocity(qq,nn,1) = -vmax;
            else
                velocity(qq,nn,1) = velocity(qq,nn,1);
            end       
        end
        velom(count,nn) = sum(abs((velocity(:,nn,1))));
        crit(:,nn) = crit(:,nn) + velocity(:,nn,1); %updates position
    end
    

    for ll = 1:rsize^2
        x3(ll) = merit_fun(func, crit(:,ll),par);
    end

    
    comp = [peval;x3]; %sets new p-best
    for iii = 1:rsize^2
        [~,II] = min(comp(:,iii));
        if II == 2
            pbest(:,iii) = crit(:,iii);
        end
    end

    count = count + 1;
end

% if count >= 1000  %use for bifucation diagrams
%   crit = NaN*crit;  
% end

% for yy = 1:xlength
%     C(yy) = mean(crit(yy,:));
% end

C = crit;



function Z = merit_fun(func,x,p) 
t = 1;
Z = .5*norm(feval(func,t,x,p))^2;

end 


end