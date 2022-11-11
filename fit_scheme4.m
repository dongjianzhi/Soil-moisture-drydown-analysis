function [r] = fit_scheme4(x,y,nc_size)
% class B
for i = 1:nc_size
    L = length(x);
    Lall = 1:length(x); r_train = sort( randsample( Lall, ceil(L*.9)) );
    r_test = setdiff(Lall,r_train);
    
    x_train = x(r_train); y_train = y(r_train);
    x_test = x(r_test); y_test = y(r_test);
    
    yhat = mean(y_train);
    
    r(i) = mean( (y_test - yhat).^2 ).^.5;
 
end