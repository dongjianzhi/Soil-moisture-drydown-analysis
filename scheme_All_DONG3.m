
function [aic,P,rmsd] = scheme_All_DONG3(x,y,knots)
% knot1 = .20; knot2 = 0.3; knots =knots_all(1,:);
% a = mean(P3); knots = [-a(2)./a(1) a(5) a(6)]
% knots = knots_all(end,:)

knot1 = knots(2); knot2 = knots(3);
theta_w = knots(1);

x1 = x(x<knot1); y1 = y(x<knot1);
x2 = x(x>= knot1 & x< knot2); y2 = y(x>= knot1 & x< knot2);
x3 = x(x>=  knot2); y3 = y(x>= knot2);

Y = [y1;  y2; y3];


if knot1<min(x) | knot2>max(x) | knot2<knot1
    aic = nan; P = zeros(1,6); rmsd = nan;
else

    if length( x1)>2 & length(x2)>2 & length(x3)>2
        
        % fit a line for cases < theta_s
        
        a = regress( y1,x1-theta_w);
        p = [a -a*theta_w];
        
        y1f = polyval(p,x1);
        y2f = ( p(1)*knot1 + p(2) ) * ones(size(x2));
        
        p2 = regress( y3 - y2f(1), x3 - knot2); y3f = p2*(x3 - knot2) + y2f(1);
        
        xhat = [x1; x2;x3]; yhat = [y1f; y2f; y3f];
        n = length(yhat);

        P = [p p2 mean(y2f) knot1 knot2];
        
        rmsd = sum( (yhat - Y).^2 );
        
        rmsd_all1 = rmsd;
        
        sep = 1;
        if sep ==3;
            a = 1./max(xhat,.05); weis = a./sum(a)*length(a);
            rmsd = sum( (yhat - Y).^2.*weis );
        end
        
        
        if sep ==2;
            n1 = length(y1);  k = 2;  rmsd1 = sum( (y1 - y1f).^2 )/n1;
            n1 = length(y2);  k = 1;  rmsd2 = sum( (y2 - y2f).^2 )/n1;
            n1 = length(y3);  k = 2;  rmsd3 = sum( (y3 - y3f).^2 )/n1;
            rmsd = (rmsd1 + rmsd2 + rmsd3)*length(yhat)/3;
        end
        
%         r = corr(yhat,Y).^2;
%         radj = 1 - ( (1-r)*(n-1) )./(n-k-1);
        
        k = 5;
        aic = n*log(rmsd/n) + (2*k + 2*k*(k+1)/( n - k-1))*0;
%         aic = n*log(rmsd/n) + 2*k;
        
%         aic = n/2*log(2*pi) + 0.5*rmsd + 2*k;
        
%         rmsd = rmsd_all1;
        
        
        if theta_w<.0; rmsd = 10^6; aic = 10^6; end
        if knot1>.4 | knot2>.4; rmsd = 10^6; aic = 10^6; end
        
%     metrics = [corr(yhat,Y),  mean( (yhat - Y).^2).^.5];

         
%         aic =6*log(n) -2*( -n/2*log(2*pi) -n/2*log(error_var) - 1/2/error_var*rmsd); % use bic instead
%         aic =12 -2*( -n/2*log(2*pi) -n/2*log(error_var) - 1/2/error_var*rmsd);

    else
        aic = nan; P = nan(1,6); rmsd = nan; metrics = nan(1,2);
    end
end


