function val=integrateInMiddleDynamics(f,eps,cEps)

val=integral(f, -sqrt(2*cEps),-sqrt(2*eps))+ integral(f, sqrt(2*eps),sqrt(2*cEps));

end