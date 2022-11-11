function AIC_calculation(mm,dd,case_use,P)

case_use = 3; P_est= mean(P3);

case_use = 1; P_est= mean(P1);

[xhat,yhat] = DONG_scheme_all(case_use,0,P_est);
y = interp1(xhat,yhat,mm);
rmsd = nanmean( (y - dd).^2 );
n = length(y); K = size(P_est,2);

aic = n*log(rmsd) + 2*K