% plotDistribution(cEps, countEps)
% clear all;
% close all;
% clc;
% 
% 
% col='*b';
% 
% useSVD=0;
% 
% eps=0.000001;
% 
% 
% 
% cEpsVect=0.0:0.1:2.5;
%cd('Data');
%indices=1:23;%length(cEpsVect);
%sg=zeros(1,length(indices));
%epsF=zeros(1,length(indices));

count=0;

ZpOrig =sqrt(2*pi);

cm=colormap('hot');
[val, idx]=max(sg);
%%
X=-3:0.01:3;

for i = indices(end:-1:1)
    
    
    count=count+1;
    cEps=cEpsVect(i);
    
    %PBC [0,2pi], V(q)=cos(q)
    s=@(x) -6*(x-eps).^5./(cEps-eps).^5+15*(x-eps).^4./(cEps-eps).^4-10*(x-eps).^3./(cEps-eps).^3+1;
    sD=@(x)-30*(x-eps).^4./(cEps-eps).^5+60*(x-eps).^3./(cEps-eps).^4-30*(x-eps).^2/(cEps-eps).^3;
    
    
    %kinetic energy
    K=@(x)0.5*x.^2;
    KD=@(x)x;
    
    
    %Z(p) in the spline part - p zeta(p)+p.^2/2 zeta'(p)
    Zfunction_splinePart=@(p) p.*s(K(p))+p.^2/2.*sD(K(p)).*p;
    %p.*(-6*(p-eps).^5./(cEps*eps-eps).^5+15.*(p-eps).^4/(cEps*eps-eps).^4-10*(p-eps).^3./(cEps*eps-eps).^3+1)+(1/2).*p.^2.*(-30.*((1/2).*p.^2-eps).^4.*p/(cEps*eps-eps).^5+60*((1/2).*p.^2-eps).^3.*p./(cEps*eps-eps).^4-30.*((1/2).*p.^2-eps).^2.*p./(cEps*eps-eps).^3);
    
    
    
    %partition function
    GZp1=@(x)exp(-0.5*x.^2);
    GZp2=@(x)exp(-0.5*(x.^2).*(1-s(K(x))));
    
    %ZpI1=sqrt(2*pi) - integral(GZp1, -sqrt(2*cEps*eps),sqrt(2*cEps*eps));
    %ZpI2a=integral(GZp2, sqrt(2*eps),sqrt(2*cEps*eps));
    %ZpI2b=integral(GZp2, -sqrt(2*cEps*eps),-sqrt(2*eps));
    
    ZpI1=sqrt(2*pi) - integral(GZp1, -sqrt(2*cEps),sqrt(2*cEps));
    ZpI2a=integral(GZp2, sqrt(2*eps),sqrt(2*cEps));
    ZpI2b=integral(GZp2, -sqrt(2*cEps),-sqrt(2*eps));
    
    
    Zp=2*sqrt(2*eps) + ZpI1 +ZpI2a+ZpI2b;
    
    
    
    
    
    GZq=@(q)exp(-cos(q));
    Zq = integral(GZq, 0, 2*pi);
    
    Z = Zp*Zq;
    
    f=zeros(1,length(X));
    for j =1:length(f)
        if(abs(X(j))> sqrt(2*cEps))
            f(j)=GZp1(X(j));
        else
            f(j)=GZp2(X(j));
        end
    end
    
    %round(i*length(cm)./length(indices))
    
    figure(2);
    hold on
    if(i == idx || i == 1)
        plot(X, f./Zp, 'Color',cm( round((i)*(length(cm)-1)./length(indices)),:), 'LineWidth',5)
    else
        plot(X, f./Zp, 'Color',cm( round((i)*(length(cm)-1)./length(indices)),:), 'LineWidth',1.2)
    
    end
    
    
end



fh=gca(figure(2));
hc=colorbar();
title(hc, '\epsilon', 'FontSize', 14)
%locate = get(hc,'title');
%pos = get(locate,'west'); %it gives a position of 0.0500 2.900 1.0001

% 
%  figure(2);
%  hold on
%  plot(X, exp(-0.5*X.^2)/ZpOrig,'-b','LineWidth',5.5 )
    

xlabel(fh, 'p')
ylabel(fh,'exp(-\beta K(p)/Z_p)')
set(fh, 'FontSize',14)

saveas(figure(2), 'distribution', 'fig')
print(figure(2), 'distribution', '-depsc')
print(figure(2), 'distribution', '-dpng')