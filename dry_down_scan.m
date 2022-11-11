function [sm11,sm12,ts1,ts2,rmsd,lambda] = dry_down_scan(sm1,ts_org,sig)

thres_range = range(sm1)/100;

sm11 = sm1(1:end-1)'; sm12 = sm1(2:end)'; 
ts1 = ts_org(1:end-1)'; ts2 = ts_org(2:end)';

det_size = 1;

while abs(det_size) > 0
    
    td = ts2 - ts1;
    lambda1 = 0; lambda2 = 1;
    
    L1 = length(sm12);
    
    for i = 1:10
        lambda = (lambda1 + lambda2)/2;
        J = sum( ( sm12 -  sm11.*lambda.^td ).*( td.*sm11.^(td-1)) );
        if J >0
            lambda1 =lambda;
        else
            lambda2 = lambda;
        end
    end
    
    y = sm11.*lambda.^td;
    err = std(sm12 - y);
    err_dis = sm12 - y;
    thres = norminv(sig,0,err);
    
    k = find(err_dis<thres);
    sm11 = sm11(k); sm12 = sm12(k); ts1 = ts1(k); ts2 = ts2(k);
    
    % remove sm change smaller than the thresholds
    d = sm11 - sm12;  t = (d - thres_range)/err*sqrt(length(sm11));
    k = find(t<-1.75); sm11(k) = []; sm12(k) = []; ts1(k) = []; ts2(k) = [];
    
    L2 = length(sm12);
    det_size = L1 - L2;
end

y = sm11.*lambda.^td;
rmsd = mean( ( y - sm12).^2 ).^.5;



