clear all;
clc;
close all;

%mkdir('Data')
cd('src')

cEpsVect =[0.0 : 0.2 : 5];
minEW = zeros(1, length(cEpsVect));
countEps=0;

%par
parfor  countEps= 1:length(cEpsVect)
    
    cEps = cEpsVect(countEps);
    
    %computeSpectralGap(cEps, countEps);
    computeSpectralGapExp(cEps, countEps);
    %minEW(countEps) = Lsorted(2);
    
end

cd ..

% figure(999)
% plot(cEpsVect, minEW, '-*b')
% xlabel('\epsilon_f')
% ylabel('Spectral gap')
