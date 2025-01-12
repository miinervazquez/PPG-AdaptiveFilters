function[w3]=slms(e3,xtdl,w3,mu3)

if e3 > 0
   sgn = 1;
elseif e3 < 0
   sgn = -1;
elseif e3 == 0
   sgn = 0;
end

w3=w3+2*mu3*sgn*xtdl;

end