function C = SPSO_Gen(func,x0,par,rsize,sigma,damping,vmax_c)

xlength = length(x0);
x = zeros(xlength,rsize^2);
vmax = vmax_c*(sigma/2); %maximum velocity. helps prevent blow up. 
radius0 = sigma;
radius = radius0;

for hh = 1:xlength
    x(hh,:) = sigma*(rand(1,rsize^2)-rand(1,rsize^2));
end

for tt = 1:rsize^2
    x(:,tt) =  x(:,tt) + x0;
end

c1 = 2;
c2 = 2;
psi = c1+c2;

crit = x;
pbest = crit;
velocity = zeros(xlength,rsize^2);

x3 = zeros(1,rsize^2);
for jj = 1:rsize^2
    x3(jj) = merit_fun(func, crit(:,jj),par);
end
smax = 10;
itermax = 10000;
iter = 1;
peval = zeros(1,rsize^2);

while mean(x3) > 1e-8 && iter < itermax
    vx3(iter,:) = x3;
    vcrit(:,:,iter) = crit;
    
    if mod(iter,100) == 0 %contracts species radius
        radius = radius*.9;
    end
    [x3,SI] = sort(x3,'ascend');
    crit = crit(:,SI);
    velocity = velocity(:,SI);
    pbest = pbest(:,SI);
    
    ii = 1;
    S = [];
    for kk = 1:rsize^2
        peval(kk) = merit_fun(func,pbest(:,kk),par);
    end
    
    while ii < length(x3) %finds species seeds
        Seed = crit(:,ii);
        if ii == 1
           S = [S Seed]; 
           ii = ii + 1;
        else
            count = 0;
            for j = 1:length(S(1,:))
                if norm(Seed-S(:,j)) < radius
                    count = count + 1;
                end
            end
            if count > 0
                ii = ii + 1;
                Seed = [];
            else
                ii = ii + 1;
                S = [S Seed];
            end
        end
    end
    
    if length(S(1,:)) > smax %sets maximum species seeds
        S = S(:,1:smax);
    end
    
    for w = 1:length(S(1,:))
        for ww = 1:length(x3)
            if norm(S(:,w)-crit(:,ww)) < radius %updates species member position
                gbest = S(:,w); %species g-best
                omega = damping;
                
                velocity(:,ww) = omega*velocity(:,ww) +...
                    c1*(eye(xlength,xlength)*rand)*(pbest(:,ww)-crit(:,ww)) +...
                    c2*(eye(xlength,xlength)*rand)*(gbest-crit(:,ww));
                for qq = 1:xlength %sets a max velocity
                    if velocity(qq,ww) > vmax
                        velocity(qq,ww) = vmax;
                    elseif velocity(qq,ww) < -vmax
                        velocity(qq,ww) = -vmax;
                    else
                        velocity(qq,ww) = velocity(qq,ww);
                    end       
                end
                crit(:,ww) = crit(:,ww) + velocity(:,ww); %updates position
            else
                count2 = 0 ;
                for www = 1:length(S(1,:)) %checks if point is inside species
                    if norm(crit(:,ww)-S(:,www)) < radius
                        count2 = count2 + 1;
                    end
                end
                if count2 == 0 %if not inside species, initiates random pos.
                    crit(:,ww) = sigma*(rand(1,xlength)-rand(1,xlength));
                end
            end
        end
    end
    

    for ll = 1:rsize^2
        x3(ll) = merit_fun(func, crit(:,ll),par);
    end
    
    comp = [peval;x3]; %updates p-best
    for iii = 1:rsize^2
        [~,II] = min(comp(:,iii));
        if II == 2
            pbest(:,iii) = crit(:,iii);
        end
    end
    iter = iter + 1;
   
end
%     x = linspace(-7,7,50); %for nice plots
%     y = x;
%     [X,Y] = meshgrid(x,y);
%     for i1 = 1:length(x)
%         for i2 = 1:length(x)
%             Z(i1,i2) = Himmelblau(t,[X(i1,i2);Y(i1,i2)],par);
%         end
%     end
        
%     v = VideoWriter('SPSO_Himmelblau_C.avi');  %for nice video
%     open(v);
%     fignum = [1 10 50 100 200 500 700 900 iter-1];
% % for ppp = 1:length(fignum)
% for ppp = 1:iter-1
%     jjj = ppp;
% %     jjj = fignum(ppp);
%     X1 = reshape(vcrit(1,:,jjj),rsize,rsize);
%     X2 = reshape(vcrit(2,:,jjj),rsize,rsize);
%     X3 = reshape(vx3(jjj,:),rsize,rsize);
% %     meshc(X,Y,Z,C)
%     contour(X,Y,Z,100)
%     hold on
% %     Z0 = Ackley(0,0);
% %     plot3(0,0,Z0,'ro')
% %     plot(0,0)
% 
%     plot3(X1,X2,X3,'k.')
% %     plot(X1,X2,'k.')
%     xlim([-7 7])
%     ylim([-7 7])
% %     view(0,90)
%     xlabel('x_{1}')
%     ylabel('x_{2}')
%     title('Himmelblau Function')
%     colormap(jet)
%     P(jjj) = getframe;
%     writeVideo(v,P(jjj));
% 
%     hold off
% %     saveas(gcf,['Himmel' num2str(jjj) '.png'])
% end
% 
% close(v);



if iter == itermax %deletes non-zero convergence.
    for gg = 1:rsize^(2)
        fx3 = merit_fun(func, crit(:,gg),par);
        if fx3 > 1e-8 
            crit(:,gg) = NaN*zeros(size(crit(:,gg)));
        end
    end
end


C = crit;



function Z = merit_fun(func,x,p)
t = 1;
Z = .5*norm(feval(func,t,x,p))^2;

end 


end