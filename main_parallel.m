clear all;
clc;
close all;

%mkdir('Data')

cEpsVect =[0.0 : 0.5 : 1 ];
minEW = zeros(1, length(cEpsVect));
countEps=0;

parfor  countEps= 1:length(cEpsVect)
    
    cEps = cEpsVect(countEps);
    
    computeSpectralGap(cEps, countEps);
    %minEW(countEps) = Lsorted(2);
    
end

% figure(999)
% plot(cEpsVect, minEW, '-*b')
% xlabel('\epsilon_f')
% ylabel('Spectral gap')
