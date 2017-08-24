%clear all;
%close all;
%clc;

cEpsVect=0.0:0.1:2.5;
cd('Data');
%cd('Results1')
indices=1:23;%
%indices=2:2:18;%length(cEpsVect)-1;
sg=zeros(1,length(indices));
epsF=zeros(1,length(indices));

count=0;
vrc=zeros(1,length(indices));

for i = indices
    
    count=count+1;
    load(['variance', num2str(i)]);
    sg(count)=real(Lcomplex(2));%Lsorted(2);%spectralGap;
    epsF(count)=cEps;
    fprintf('epsF = %f\nvariance is %f\n', cEps, variance)
    vrc(count)=real(variance);
    
end
%%
f1=figure(1);
hold on
plot(epsF, sg, '-*b', 'LineWidth',2)
xlabel('\epsilon')
ylabel('Spectral gap')
set(gca, 'FontSize',14)

saveas(f1, 'spectralGap', 'fig')
print(f1, 'spectralGap', '-depsc')

%%

f10=figure(10);
plot(epsF, vrc, '-*b', 'LineWidth',2)
xlabel('\epsilon')
ylabel('Variance of cos(q)')
set(gca, 'FontSize',14)
%%
cd ..
plotDistribution;
