function plm=Plm(x,l,m)

    if abs(m)>l
        plm=zeros(1,numel(x));
    elseif m>=0
        Pl=legendre(l,x);
        plm=Pl(m+1,:);
            
    else        
        Pl=legendre(l,x);
        plm=(-1)^(-m)*gamma(l+m+1)/gamma(l-m+1).*Pl(-m+1,:);
    end
    
end