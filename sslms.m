function[w3]=sslms(e3,xtdl,w3,mu3)

x=mean(xtdl);

if x > 0
   sgn = 1;
elseif x < 0
   sgn = -1;
elseif x == 0
   sgn = 0;
end

if e3 > 0
   sgne = 1;
elseif e3 < 0
   sgne = -1;
elseif e3 == 0
   sgne = 0;
end

w3=w3+2*mu3*sgne*sgn;


end