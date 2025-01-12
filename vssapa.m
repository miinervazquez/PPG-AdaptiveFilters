function[w2]=vssapa(w2,mu2max,X2,E2,C2,alpha2,P2,M,gamma)

A = (X2 * inv( X2' * X2 +gamma*eye(M))) * E2 ;
P2 = alpha2 * P2 + (1 - alpha2) * A ;
p = norm( P2 ) ^2 ;
mu = mu2max * ( p / ( p + C2) ) ;
w2 = w2 + mu * A ;
end

% 0 < alpha < 1
% C is a positive constant.
% mu2max is chosen less than 2