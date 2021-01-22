function psi2 = Psi2(l,x)

psi2=-Psi(l,x).*(1-l*(l+1)./(x.^2));

end