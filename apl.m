function[w9]=apl(w9,mu,X9,E9)

P = X9 * E9 ;
d = norm(X9' * P)^2 ;

if d == 0 
    d = 1;
end

w9 = w9 + mu * ( (norm(P)^2 ) / ( d ) ) * ( P ) ;

end
