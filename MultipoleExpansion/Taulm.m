function taulm=Taulm(theta,l,m)
%     taulm=(l-m+1)*(1./sin(theta)).*Plm(cos(theta),l+1,m)'-(l+1)*(cos(theta)./sin(theta)).*Plm(cos(theta),l,m);
    taulm=(l*cos(theta).*Plm(cos(theta),l,m)-(l+m).*Plm(cos(theta),l-1,m))./sin(theta);
end