function computeSpectralGap(cEps, countEps)

N_G = 12;
N_H = 2*N_G-1;


col='*b';

useSVD=0;

eps=0.000001;


% ------------ A(q)=cos(q)

Aobservable=@(q)cos(q);
%------------------------
%----------------------------
%------------- A(q)=G_1(q)=e^(iq)
%averA=0;
%Aobservable=@(q)exp(1i*q);

PBC2=4*pi;


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


%PBC [0,PBC2], V(q)=cos(q)



%zeta spline
%s=@(x) -6*(x-eps).^5./(cEps*eps-eps).^5+15*(x-eps).^4./(cEps*eps-eps).^4-10*(x-eps).^3./(cEps*eps-eps).^3+1;
%sD=@(x)-30*(x-eps).^4./(cEps*eps-eps).^5+60*(x-eps).^3./(cEps*eps-eps).^4-30*(x-eps).^2/(cEps*eps-eps).^3;
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



ZpOrig =sqrt(2*pi);

GZq=@(q)exp(-cos(q));
Zq = integral(GZq, 0, PBC2);

Z = Zp*Zq;


GA=@(q)  Aobservable(q).*exp(-cos(q));
averA = integral(GA, 0, 4*pi)/Zq;




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
            
            Glk(j+1,k+1)= integral(G1, 0, PBC2)/Zq;
            
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
            
            JFull=@(p)Hn(p).*exp(-p.^2/2);
            JSpline = @(p) Hn(p).*exp(-(1-s(K(p))).*p.^2/2);
            JRestr = @(p) Hn(p);
            
            Jtmp=sqrt(2*pi)*KroneckerDelta(n,0)  - integrateInFullDynamics(JFull,eps,cEps) + integrateInRestrainedDynamics(JRestr,eps)+integrateInMiddleDynamics(JSpline,eps,cEps);
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
                    
                    
                    pHnHmFull = @(p)p.*Hn(p).*Hm(p).*exp(-p.^2/2);
                    pMinusZpHnHmSpline = @(p)(p-Zfunction_splinePart(p)).*Hn(p).*Hm(p).*exp(-(1-s(K(p))).*p.^2/2);
                    
                    HnHmFull = @(p)Hn(p).*Hm(p).*exp(-p.^2/2);
                    HnHmRestrained = @(p)Hn(p).*Hm(p);
                    HnHmSpline= @(p)Hn(p).*Hm(p).*exp(-(1-s(K(p))).*p.^2/2);
                    
                    pHnHmMinus1Restrained= @(p)p.*Hm(p).*HnMinus1(p);
                    ZHnHmSpline= @(p)Zfunction_splinePart(p).*Hn(p).*Hm(p).*exp(-(1-s(K(p))).*p.^2/2);
                    ZHnHmMinus1Spline= @(p)Zfunction_splinePart(p).*HnMinus1(p).*Hm(p).*exp(-(1-s(K(p))).*p.^2/2);
                    
                    HnHmMinus1Restrained= @(p)Hm(p).*HnMinus1(p);
                    HnHmMinus1Midlde=@(p)Hm(p).*HnMinus1(p).*exp(-(1-s(K(p))).*p.^2/2);
                    HnHmMinus1Full=@(p)Hm(p).*HnMinus1(p).*exp(-p.^2/2);
                    
                    I1=sqrt(2*pi)*(sqrt(n+1)*KroneckerDelta(n+1,m)+sqrt(n)*KroneckerDelta(n-1,m)) - integrateInFullDynamics(pHnHmFull,eps,cEps) + integrateInMiddleDynamics(pMinusZpHnHmSpline,eps,cEps);
                    
                    I2= integrateInRestrainedDynamics(HnHmRestrained,eps) + sqrt(2*pi)*KroneckerDelta(n,m) - integrateInFullDynamics(HnHmFull,eps,cEps) + integrateInMiddleDynamics(HnHmSpline,eps,cEps);
                    
                    I3= integrateInRestrainedDynamics(pHnHmMinus1Restrained,eps) + integrateInMiddleDynamics(ZHnHmMinus1Spline,eps,cEps);
                    
                    I4=integrateInRestrainedDynamics(HnHmMinus1Restrained, eps) + integrateInMiddleDynamics(HnHmMinus1Midlde, eps, cEps) + sqrt(2*pi)*KroneckerDelta(n-1,m)-integrateInFullDynamics(HnHmMinus1Full,eps,cEps);
                    
                    
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
            intAGk = integral(G, 0, PBC2);
            
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
                    
                    HnHmFull = @(p)Hn(p).*Hm(p).*exp(-p.^2/2);
                    HnHmRestrained = @(p)Hn(p).*Hm(p);
                    HnHmSpline= @(p)Hn(p).*Hm(p).*exp(-(1-s(K(p))).*p.^2/2);
                    
                    
                    %I2 from Q
                    IHnHm =  integrateInRestrainedDynamics(HnHmRestrained,eps) + sqrt(2*pi)*KroneckerDelta(n,m) - integrateInFullDynamics(HnHmFull,eps,cEps) + integrateInMiddleDynamics(HnHmSpline,eps,cEps);
                    
                    
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


[Lsorted,idxreal] = sort(real(L));
spectralGap = Lsorted(2);
Lcomplex = L(idxreal);

save(['Data/variance', num2str(countEps)], 'variance', 'Lcomplex','Lsorted', 'spectralGap', 'eps', 'cEps', 'N_H', 'N_G')


end
