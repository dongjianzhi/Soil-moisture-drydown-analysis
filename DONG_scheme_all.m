function [xhat,yhat] = DONG_scheme_all(case_use,y,P)

xhat = [0:.001:.5]';
% xhat = sort(rand(100,1)*.5);

if case_use ==1;  yhat = xhat*P(1) + P(2); end
if case_use ==6;  yhat = xhat*P(1) + P(2); end


if case_use ==2;
    knot1 = P(3);
    x1 = xhat(xhat<knot1); 
    x2 = xhat(xhat>= knot1 ); 
    yhat = xhat*P(1) + P(2);
    
    y1f = polyval(P(1:2),x1);
    y2f = ( P(1)*knot1 + P(2) ) * ones(size(x2));
    
    xhat = [x1; x2]; yhat = [y1f; y2f];
end

if case_use ==3;
    knot1 = P(end-1); knot2 = P(end); plat = P(4);
    
    x1 = xhat(xhat<knot1); 
    x2 = xhat(xhat>= knot1 & xhat< knot2); 
    x3 = xhat(xhat>=  knot2); 
    
    p = P(1:2); p2 = P(3);
    y1f = polyval(p,x1);
    y2f = ( p(1)*knot1 + p(2) ) * ones(size(x2));
    y3f = p2*(x3 - knot2) + y2f(1);
    
    xhat = [x1; x2; x3]'; yhat = [y1f; y2f; y3f]';
end

if case_use ==4;   yhat = mean(y)*ones(size(xhat)); end


if case_use == 5
    knot1 = P(3); yflat = P(2); p = P(1);
    x1 = xhat(xhat<knot1);
    x2 = xhat(xhat>= knot1);
    
    y1f = yflat*ones(size(x1));
    y2f = (x2 - knot1)*p + y1f(end);
    xhat = [x1; x2;]'; yhat = [y1f; y2f;]';

end

xhat = reshape( xhat,length(xhat),1);
yhat = reshape( yhat,length(yhat),1);





