function [aic,P,rmsd] = scheme_All_DONG6(x,y,knots)
% knots = .05
knot1 = knots(1);
L = length(x);

y1 = y-0; x1 = x-knot1;    
a = regress( y1,x1,1); 

% y = a*(x - knot1); y0 = -a*knot1
p = [a -a*knot1];


yhat = polyval(p,x);

P = p; n = length(yhat);

rmsd = sum( (yhat - y).^2 );

sep = 1;
if sep ==3;
    a = 1./max(x,.05); weis = a./sum(a)*length(a);
    rmsd = sum( (yhat - y).^2.*weis );
end
    
% aic =4 -2*( -n/2*log(2*pi) -n/2*log(error_var) - 1/2/error_var*rmsd);

k = 2;
aic = n*log(rmsd/n) +  (2*k + 2*k*(k+1)/( n - k-1))*0;

if knot1<.0; rmsd = 10^6; aic = 10^6; end
