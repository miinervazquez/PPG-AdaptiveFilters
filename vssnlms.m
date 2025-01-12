function[w3]=vssnlms(e3,xtdl,w3,mumax,C,alpha,P1)

x2 = norm(xtdl)^2;

if x2==0
    x2=1;
end

P1 = alpha .* P1 +  (1-alpha) .* xtdl .* e3 ./ x2 ;

P = norm(P1)^2;

mu = mumax * ( P / ( P + C) );

w3=w3+mu*( xtdl / x2 ) * e3 ;
end