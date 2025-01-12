function[w9]=apa(w9,mu9,X9,E9,M,gamma)

w9=w9+mu9*X9*inv(X9'*X9+gamma*eye(M))*E9;

end