
function [aic,P,rmsd] = scheme_All_DONG5(x,y,knots)
% knots = knots_all(1,:);
knot1 = knots(1);

x1 = x(x<knot1); y1 = y(x<knot1);
x2 = x(x>= knot1 ); y2 = y(x>= knot1 );


Y = [y1;  y2]; 

if length( x1)>2 & length(x2)>2
    
    y1f = mean(y1)*ones(size(x1));
    p = regress( (y2 - y1f(end)), (x2 - knot1));
    y2f = (x2 - knot1)*p + y1f(end);
    
    
    xhat = [x1; x2]; yhat = [y1f; y2f];
    n = length(yhat);
    
    P = [p  mean(y1f) knot1 ]; n = length(yhat);
    
    rmsd = sum( (yhat - Y).^2 );
    rmsd_all = rmsd;
    
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
    aic = n*log(rmsd/n) + ( 2*k + 2*k*(k+1)/( n - k-1))*0;
    
%     rmsd = rmsd_all;

    if knot1>.4; rmsd = 10^6; aic = 10^6; end
    


%     aic =8 -2*( -n/2*log(2*pi) -n/2*log(error_var) - 1/2/error_var*rmsd);
    
else
    aic = nan; P = nan(1,3); rmsd = nan; metrics= nan(1,2);
end

