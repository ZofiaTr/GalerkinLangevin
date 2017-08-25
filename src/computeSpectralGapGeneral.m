function [L, variance]= computeSpectralGapGeneral(cEps, countEps)
%%

N_G = 6;%12;
N_H = 2*N_G-1;


col='*b';

useSVD=0;

eps=0.00001;


% ------------ A(q)=cos(q)

Aobservable=@(q)cos(q);
%------------------------
%----------------------------
%------------- A(q)=G_1(q)=e^(iq)
%averA=0;
%Aobservable=@(q)exp(1i*q);




%cEps=1;

plotSpectrum=0;
computeQ=1;
computeA=1;
computeS=1;
computeQoriginal=0;

saveFigureSpectrum=0;
%useSVD=1;



if(computeQoriginal == 1)
    
    main_galerkin_cosPotential;
    
    varianceOriginal=variance;
    Qoriginal=Q;
    Soriginal=S;
    Aoriginal=A;
    
end


%PBC [0,2pi], V(q)=cos(q)




%kinetic energy
power=cEps; % must be parny
K=@(x)abs(x).^(power)./power;%0.5*x.^2 + cEps*cos(x);% cEps.*exp(-0.5*x.^2);
KD=@(x)sign(x).*abs(x).^(power-1);% x -cEps* sin(x);% cEps.*x.*exp(-0.5*x.^2);


%partition function
GZp1=@(x)exp(-K(x));

Zp=integral(GZp1, -inf,inf);




ZpOrig =sqrt(2*pi);

GZq=@(q)exp(-cos(q));
Zq = integral(GZq, 0, 2*pi);

Z = Zp*Zq;


GA=@(q)  Aobservable(q).*exp(-cos(q));
averA = integral(GA, 0, 2*pi)/Zq;




Mj = 2*N_G-1;
Mn = N_H;

totalBasisSize=Mj*Mn-1;

Q = zeros(totalBasisSize,totalBasisSize);
gamma =1;



if(computeA==1)
    J=zeros(1,totalBasisSize+1);
    
    %doesExist=exist('Glk');
    %if(doesExist==0)
    
    Glk=zeros(Mj,Mj);
    k=0;
    j=0;
    
    
    
    for kIndex=-N_G+1:1:N_G-1
        
        j=0;
        
        for jIndex=-N_G+1:1:N_G-1
            
            
            
            G1=@(q) exp(1i*q*(jIndex-kIndex)).*exp(-cos(q));
            
            Glk(j+1,k+1)= integral(G1, 0, 2*pi)/Zq;
            
            j=j+1;
            
            
            
        end
        k=k+1;
    end
    
    %end
    
    
    for n =0:N_H-1
        k = 0;
        
        for kIndex=-N_G+1:1:N_G-1
            
            if(kIndex==0 && n==0)
                indexA00=k*Mn+n+1;
            end
            
            
            Hn =  H(n);
            
            
            JSpline = @(p) Hn(p).*exp(-K(p));
            
            Jtmp=integral(JSpline, -inf, inf);
            J(k*Mn+n+1) = Jtmp;
            
            
            
            k = k+1;
            
            
            
        end
        
    end
    
    %  J=J([1:indexA00-1,indexA00+1:end]);
end

if(computeQ == 1)
    
    for n =0:N_H-1
        k = 0;
        
        for kIndex=-N_G+1:1:N_G-1
            
            
            if(kIndex==0 && n==0)
                indexQ00=k*Mn+n+1;
            end
            
            
            for m=0: N_H-1
                
                j = 0;
                
                
                for jIndex =-N_G+1:1:N_G-1
                    
                    
                    if(jIndex==0 && m==0)
                        
                        indexQ00m=j*Mn+m+1;
                    end
                    
                    Hn=H(n);
                    Hm=H(m);
                    HnMinus1=H(n-1);
                    % HmPlus1=H(m+1);
                    % HmMins1=H(m-1);
                    
                    
                    pHnHmFull = @(p)p.*Hn(p).*Hm(p).*exp(-K(p));
                    
                    HnHmFull = @(p)Hn(p).*Hm(p).*exp(-K(p));
                    HnHmMinus1Full=@(p)Hm(p).*HnMinus1(p).*exp(-K(p));
                    
                    I1= integral(pHnHmFull,-inf,inf);% + integrateInMiddleDynamics(pMinusZpHnHmSpline,eps,cEps);
                    
                    I2=  integral(HnHmFull,-inf,inf) ;
                    
                    I3= 0.0;% integral(pHnHmMinus1Fu,-inf, inf) ;
                    
                    I4=integral(HnHmMinus1Full,-inf,inf);
                    
                    
                    Q1= -1i*kIndex*Glk(j+1,k+1)*I1/Zp;
                    
                    Q2=-n*gamma*Glk(j+1,k+1)*I2/Zp;
                    
                    Q3= gamma*sqrt(n)*Glk(j+1,k+1)*I3/Zp;
                    
                    Q4=-1i*(jIndex-kIndex)*sqrt(n)*Glk(j+1, k+1)*I4/Zp;
                    
                    %Qtmp = - ((sqrt(n+1)*KroneckerDelta(n+1,m)+sqrt(n)*KroneckerDelta(n-1,m)) *(-1i*kIndex * KroneckerDelta(jIndex,kIndex)) - gamma*n*KroneckerDelta(m,n)*KroneckerDelta(jIndex,kIndex));
                    
                    Q(k*Mn+n+1, j*Mn+m+1) = -( Q1 + Q2 + Q3 + Q4);
                    
                    
                    
                    
                    
                    j=j+1;
                    
                end
                
                
            end
            
            k = k+1;
        end
    end
    
    
end

fprintf('Q done\n')



if(computeA == 1)
    
    A=zeros(Mj*Mn-1,1);
    
    
    for n =0:N_H-1
        k = 0;
        
        for kIndex=-N_G+1:1:N_G-1
            
            if(kIndex==0 && n==0)
                indexA00=k*Mn+n+1;
            end
            
            G=@(q)(Aobservable(q)-averA).*exp(-1i.*q*kIndex-cos(q) );
            intAGk = integral(G, 0, 2*pi);
            
            IHn=J(k*Mn+n+1);
            
            % Atmp=IHn*pi*(KroneckerDelta(kIndex,1)+KroneckerDelta(kIndex,-1))/Z;
            Atmp=IHn*intAGk/Z;
            
            
            A(k*Mn+n+1) = Atmp;
            
            
            k = k+1;
            
            
        end
        
    end
    
end

fprintf('A done\n')

if(computeS == 1)
    
    S = zeros(Mj*Mn-1, Mj*Mn-1);
    
    for n =0:N_H-1
        k = 0;
        
        for kIndex=-N_G+1:1:N_G-1
            
            if(kIndex==0 && n==0)
                indexS00=k*Mn+n+1;
            end
            
            
            for m=0: N_H-1
                
                j = 0;
                
                for jIndex =-N_G+1:1:N_G-1
                    
                    if(jIndex==0 && m==0)
                        
                        indexS00m=j*Mn+m+1;
                    end
                    
                    Hn=H(n);
                    Hm=H(m);
                    
                    HnHmFull = @(p)Hn(p).*Hm(p).*exp(-K(p));
                    
                    
                    %I2 from Q
                    IHnHm =  integral(HnHmFull,-inf,inf) ;
                    
                    
                    S(k*Mn+n+1, j*Mn+m+1) = IHnHm*Glk(j+1,k+1)/Zp;
                    
                    j=j+1;
                    
                end
            end
            
            k = k+1;
        end
    end
    
    
end
fprintf('S done\n')

L=eig(S\Q);

if(plotSpectrum==1)
    
    
    epsStr= num2str(eps);
    
    if(computeQoriginal==1)
        
        Lorig=eig(Soriginal\Qoriginal);
        
        figure(162)
        plot(real(Lorig), imag(Lorig),'*r',real(L), imag(L),'*b')
        legend({['\epsilon=',epsStr],'\epsilon=0'})
        
    else
        figure(12)
        hold on
        plot(real(L), imag(L),'*b')
        
    end
    
    if(saveFigureSpectrum == 1)
        
        NStrH=num2str(N_H);
        NStrG=num2str(N_G);
        strEps=num2str(eps);
        NStr=['N_H_',NStrH,'N_G',NStrG,'_',strEps];
        saveas(gcf, NStr,'png')
        
    end
end



if(useSVD==0)
    
    
    SnonProj=S;
    QnonProj=Q;
    AnonProj=A;
    
    S=S([1:indexS00-1,indexS00+1:end],[1:indexS00m-1,indexS00m+1:end]);
    Q=Q([1:indexQ00-1,indexQ00+1:end],[1:indexQ00m-1,indexQ00m+1:end]);
    A=A([1:indexA00-1,indexA00+1:end]);
    
    
    if(computeA==1 && computeS==1)
        % A_ml= (S')\A;
        
        
        if(det(Q)==0)
            fprintf('det(Q)=0  !!!! ')
        else
            F=Q'\A;
        end
        
        %
        % F_ml= (S')\F;
        
        
        variance = 2* F'*A;
        %            varianceInNhNG(NHcount,NGcount)=variance;
        
        fprintf(': N_H=%d, N_G=%d \n ARPS epsilon=%f using proj 0-EV: Variance is %f +i%f \n',N_H, N_G, eps, real(variance), imag(variance));
        
        Ssvd=SnonProj;
        Qsvd=QnonProj;
        Asvd=AnonProj;
        
        [Usvd,Dsvd,Vsvd] = svd(Qsvd);
        
        tmp= (Usvd'*(Asvd))./diag(Dsvd);
        dsvd=diag(Dsvd);
        
        indexCount=0;
        for index_i=length(dsvd):-1:1
            if(abs(dsvd(index_i))<10^(-12))
                tmp(index_i)=0;
                indexCount=indexCount+1;
                
            end
        end
        
        %     if(indexCount>1)
        %
        %           Qsvd=Qsvd([1:indexQ00-1,indexQ00+1:end],[1:indexQ00m-1,indexQ00m+1:end]);
        %           Asvd=Asvd([1:indexA00-1,indexA00+1:end]);
        %
        %     end
        %
        %     [Usvd,Dsvd,Vsvd] = svd(Q);
        %
        %      tmp= (U'*(A))./diag(D);
        %     d=diag(D);
        %
        %      for index_i=length(d):-1:1
        %         if(abs(d(index_i))<10^(-12))
        %             tmp(index_i)=0;
        %
        %         end
        %     end
        
        Fsvd= Vsvd*tmp;
        variance_svd=2*Fsvd'*Asvd;
        
        
        fprintf('N_H=%d, N_G=%d \n ARPS epsilon=%f using SVD:  Variance is %f +(%f i) \n', N_H, N_G,eps, real(variance_svd), imag(variance_svd));
        
    end
    
else
    
    QnotProjected=Q;
    AnotProjected=A;
    %
    %     S=SnonProj;
    %  Q=QnonProj;
    %  A=AnonProj;
    
    [U,D,V] = svd(Q);
    
    tmp= (U'*(A))./diag(D);
    d=diag(D);
    
    indexCount=0;
    for index_i=length(d):-1:1
        if(abs(d(index_i))<10^(-12))
            tmp(index_i)=0;
            indexCount=indexCount+1;
            
        end
    end
    
    if(indexCount>1)
        
        Q=Q([1:indexQ00-1,indexQ00+1:end],[1:indexQ00m-1,indexQ00m+1:end]);
        A=A([1:indexA00-1,indexA00+1:end]);
        
    end
    
    [U,D,V] = svd(Q);
    
    tmp= (U'*(A))./diag(D);
    d=diag(D);
    
    for index_i=length(d):-1:1
        if(abs(d(index_i))<10^(-12))
            tmp(index_i)=0;
            
        end
    end
    
    F= V*tmp;
    variance=2*F'*A;
    
    
    fprintf('N_H=%d, N_G=%d \n ARPS epsilon=%f using SVD:  Variance is %f +(%f i) \n', N_H, N_G, cEps,real(variance), imag(variance));
end


Lsorted = sort(real(L));
spectralGap = Lsorted(2);

save(['DataTmp/variance', num2str(countEps)], 'variance', 'Lsorted', 'spectralGap', 'eps', 'cEps', 'N_H', 'N_G')

%%
end