
function [aic,P,rmsd] = scheme_All_DONG2(x,y,knots)
% knots = [.05 .2]; a = mean(P2); knots = [ -a(2)./a(1) a(end)];
knot1 = knots(2);
L = length(x);

theta_w = knots(1);

x1 = x(x<knot1); y1 = y(x<knot1);
x2 = x(x>= knot1 ); y2 = y(x>= knot1 );

Y = [y1;  y2];

if length( x1)>2 & length(x2)>2
    
    a = regress( y1,x1-theta_w);
    p = [a -a*theta_w];
    
    y1f = polyval(p,x1);
    y2f = ( p(1)*knot1 + p(2) ) * ones(size(x2));
    
    xhat = [x1; x2]; yhat = [y1f; y2f]; 
    n = length(yhat);
    
    P = [p knot1 ]; n = length(yhat);
    
    rmsd = sum( (yhat - Y).^2 );
    
    
    
    sep = 1;
    if sep ==3;
        a = 1./max(xhat,.05); weis = a./sum(a)*length(a);
        rmsd = sum( (yhat - Y).^2.*weis );
    end
    
    if sep ==2;
            n1 = length(y1);  k = 2;  rmsd1 = sum( (y1 - y1f).^2 )/n1;
            n1 = length(y2);  k = 1;  rmsd2 = sum( (y2 - y2f).^2 )/n1;
            rmsd = (rmsd1 + rmsd2)*length(yhat)/2;
    end

    

    k = 3;  
    aic = n*log(rmsd/n) + (2*k + 2*k*(k+1)/( n - k-1))*0;

%     aic = n/2*log(2*pi) + 0.5*rmsd + 2*k;
    
    if theta_w<.0; rmsd = 10^6; aic = 10^6; end
    if knot1>.4; rmsd = 10^6; aic = 10^6; end
    
else
    aic = nan; P = nan(1,3); rmsd = nan; metrics = nan(1,2);
end

