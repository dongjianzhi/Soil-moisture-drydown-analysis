function  [P, rall]= fit_All_scheme6(x,y,LOOP,nc_size,error_var)

x = x(~isnan(x)); y = y(~isnan(y));

tlow = 0.15; tup = .40;% tup = min( 0.25, prctile(x,30));

knots_all = rand(nc_size*1.3,1) * (tup - tlow) + tlow;
knots_prior = knots_all;

for J = 1:LOOP
    rall = nan(length(knots_all),1);
    P = nan(length(knots_all),2);
    rmsd = rall;
    
%     tic
    for i = 1:length(knots_all)
        try
            [rall(i),P(i,:),rmsd(i)] = scheme_All_DONG6(x,y,knots_all(i,:));
        end
    end
%     toc
    
    k = find(isnan(rall) | isnan(P(:,1))) ;
    rall(k) = []; P(k,:) = []; knots_all(k,:) = []; rmsd(k) = [];
    
    % calculate likelihood
%     k = find(rmsd<10^5); rm = min(rmsd(k));
%     lik = exp( -rmsd);     lik = max(lik,10^-99);
%     
%     if exp(-rm)<10^-50; sz = 100; lik = exp(-rmsd/sz); lik = max(lik,10^-99); end
%     wei = lik./sum(lik);
    k = find(isnan(rmsd)); rmsd(k) = 10^6;
    wei = likli_adjust(rmsd);
    
    % resample
    cwei = cumsum(wei); L = length(cwei);
    u = linspace( 0.01,.999,L);
    clear ki
    for i = 1:L
        k = find( cwei> u(i),1);
        knots_all_resample(i,:) = knots_all(k,:);
        kk(i) = k;
    end
    k = randsample(L,L,'true');
    knot_rec = knots_all;
    knots_all = knots_all_resample*.9 + knots_all_resample(k,:)*.1 ...
        + (max(x) - min(x))*randn(size(knots_all_resample))*.01;
    
    % test convergence
    t = ttest_mine(knots_all,knots_prior);
    knots_prior = knots_all;
    if abs(t)<1.732 & J>5
        knots_all = knot_rec;
        return
    end

    clear knots_all_resample
end  