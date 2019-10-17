clear all

N = 95; %grid size

x = linspace(-50,50, N);

y = x;

[X,Y] = meshgrid(x);

Z = Ackley(X,Y);

C = gradient(Z);

rsize = 10;
sigma = 50;
damping = 0.00001;
vmax = 0.1*(sigma)/2;

x1 = sigma*(rand(1,rsize^2)-rand(1,rsize^2));
x2 = sigma*(rand(1,rsize^2)-rand(1,rsize^2));

c1 = 0;
c2 = 2;

crit = [x1;x2];
xlength = (length(crit(:,1)));
pbest = crit;
velocity = zeros(2,rsize^2,2);

x3 = Ackley(crit(1,:),crit(2,:));
count = 1;
tic
while norm(x3) > 1e-5
    vx3(count,:) = x3;
    vcrit(:,:,count) = crit;

    peval = Ackley(pbest(1,:),pbest(2,:));
    [~,I] = min(peval);
    gbest = crit(:,I); %sets gbest
    
    for nn = 1:rsize^2
        rr1 = rand;
        omega = .005*(rr1/2); %random damping
        velocity(:,nn,1) = omega*velocity(:,nn,1) +...
            c1*(eye(2,2)*rand)*(pbest(:,nn)-crit(:,nn)) +...
            c2*(eye(2,2)*rand)*(gbest-crit(:,nn));
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
        crit(:,nn) = crit(:,nn) + velocity(:,nn,1);
    end

    x3 = Ackley(crit(1,:),crit(2,:));

    
    comp = [peval;x3]; %sets new pbest
    for iii = 1:rsize^2
        [~,II] = min(comp(:,iii));
        if II == 2
            pbest(:,iii) = crit(:,iii);
        end
    end

    count = count + 1;
end

toc

% v = VideoWriter('PSO_Ackley_C.avi');
% open(v);

% for jjj = 1:count-1
% fignum = [1 3 5 10 15 20 30];
for ppp = 1:length(fignum)
    jjj = fignum(ppp);
    X1 = reshape(vcrit(1,:,jjj),rsize,rsize);
    X2 = reshape(vcrit(2,:,jjj),rsize,rsize);
    X3 = reshape(vx3(jjj,:),rsize,rsize);
%     meshc(X,Y,Z,C)
    contour(X,Y,Z,20)
    hold on
    Z0 = Ackley(0,0);
    plot3(0,0,Z0,'ro')
%     plot(0,0)

    plot3(X1,X2,X3,'k.')
%     plot(X1,X2,'k.')
    xlim([-50 50])
    ylim([-50 50])
%     view(0,90)
    xlabel('x_{1}')
    ylabel('x_{2}')
    title('Ackleys function')
    colormap(jet)
%     P(jjj) = getframe;
%     writeVideo(v,P(jjj));

    hold off
%     saveas(gcf,['stuff' num2str(jjj) '.png'])
end

% close(v);



function Z = Ackley(x,y)

Z = -20*exp((-0.2)*sqrt(.5*(x.^2 + y.^2)))...
    - exp(.5*(cos(2*pi*x) + cos(2*pi*y))) + 20 + exp(1);

end