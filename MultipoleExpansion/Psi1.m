function psi1 = Psi1(l,x)

psi1=Psi(l-1,x)-l.*sphericalBessel(l,x);

end

function j = sphericalBessel(l,x)

j=(pi./(2.*x)).^(1/2).*besselj(l+0.5,x);

end