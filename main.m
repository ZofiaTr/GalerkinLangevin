clear all;
clc;
close all;

%mkdir('Data')

cEpsVect =1.1:0.1:3;%[0.0 : 1 : 5 ];
minEW = zeros(1, length(cEpsVect));

vrc = zeros(1, length(cEpsVect));
countEps=0;


cd('src');
%par
for  countEps= 1:length(cEpsVect)
    
    cEps = cEpsVect(countEps);
    
    %computeSpectralGap(cEps, countEps);
    %[spectrum, variance]=computeSpectralGapExp(cEps, countEps);
    [spectrum, variance]=computeSpectralGapGeneral(cEps, countEps);
    [Ls, idx]=sort(real(spectrum));
    Lsorted=spectrum(idx);
    minEW(countEps) = Lsorted(2);
    vrc(countEps) = variance;
    
end

cd ..

%%
figure(999)
plot(cEpsVect, minEW, '-*b')
xlabel('\epsilon_f')
ylabel('Spectral gap')


figure(998)
plot(cEpsVect, vrc, '-*b')
xlabel('\epsilon_f')
ylabel('Variance cos(q)')