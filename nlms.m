function[w3]=nlms(e3,xtdl,w3,mu3,gamma)

mu = (mu3)/(gamma+(xtdl'*xtdl));
w3=w3+mu*e3*xtdl;

end