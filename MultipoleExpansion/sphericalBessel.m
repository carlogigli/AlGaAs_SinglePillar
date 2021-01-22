function j = sphericalBessel(l,x)

j=(pi./(2.*x)).^(1/2).*besselj(l+0.5,x);

end
