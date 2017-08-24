function val=integrateInFullDynamics(f,eps,cEps)

%val=integral(f, -sqrt(2*eps*cEps),sqrt(2*eps*cEps));
val=integral(f, -cEps,cEps);
end