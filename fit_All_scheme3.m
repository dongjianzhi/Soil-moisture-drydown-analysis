function [P,rall,metrics] = fit_All_scheme3(x,y,LOOP,nc_size)

x = x(~isnan(x)); y = y(~isnan(y));
th_s_fc_diff = 0.075; % threshold for theta_fc and theta_w differences

% sample theta_w, theta_s, theta_fc
knots_all = nan(1,3);
tlow = 0.0; tup = .25; % tup = min( 0.25, prctile(x,10));
theta_w = rand(nc_size*3,1) * (tup - tlow) + tlow;

a = sort(x);  as = a(4); ae = a(find( a> a(end-4)-.075,1)); ae2 = a(end-4);
s1 =(rand( round(nc_size*3),1) )*( ae - as)+as; 
s2 = s1 + ( rand(size(s1)) .*(ae2 - th_s_fc_diff -s1) +th_s_fc_diff);
knots_all = [theta_w s1 s2];

% knots_all(end,:) = P_true(end,1:3);

knots_prior = knots_all;
min_loop = 3;
clear knots_all_resample

for J = 1:LOOP
    rall = nan(length(knots_all),1);
    P = nan(length(knots_all),6);
    rmsd = rall; 
    
%     tic
%     knots_all(i,:) = [0.0968    0.2314    0.3377];
    for i = 1:length(knots_all)
        try
            [rall(i),P(i,:),rmsd(i)] = scheme_All_DONG3(x,y,knots_all(i,:));
        end
    end
%     toc
    
    k = find(isnan(rall) | isnan(P(:,1))) ;
    
    if size(k,1) == size(rall,1)
        break
    end
    
    
    rall(k) = []; P(k,:) = []; knots_all(k,:) = []; rmsd(k,:) = [];
    k = find(P(:,1)<0); rmsd(k,:) = 10^6;
    k = find(P(:,5)<prctile(x,10) | P(:,6)> prctile(x,95) ); rmsd(k,:) = 10^6;
    k = find(P(:,6) - P(:,5)<.05 |  P(:,3)<P(:,1)/10  ); rmsd(k,:) = 10^6;
%     k = find(P(:,3)<P(:,1)); rmsd(k,:) = 10^6;
    
    % calculate likelihood
%     lik = exp( -rmsd);     lik = max(lik,10^-99);
%     wei = lik./sum(lik);
%     
    % modified DONG Aug 31
%     k = find(rmsd<10^5); rm = min(rmsd(k));
%     if exp(-rm)<10^-50; sz = 100; lik = exp(-rmsd/sz); lik = max(lik,10^-99); end
%     k = find(P(:,1)<0 | P(:,3)<0); lik(k) = 0;
%     
%     wei = lik./sum(lik);
    k = find(isnan(rmsd)); rmsd(k) = 10^6;
    wei = likli_adjust(rmsd);

    
%     if std(wei) <.001; lik = exp(-rmsd/3); lik = max(lik,10^-99); wei = lik./sum(lik);end
    % resample
    cwei = cumsum(wei); L = length(cwei);
    u = linspace( 0.01,.999,L);
    clear ki
    for i = 1:L
        k = find( cwei> u(i),1);
        knots_all_resample(i,:) = knots_all(k,:);
        ki(i) = k;
    end
    k = randsample(L,L,'true');
    a = knots_all;
    knots_all = knots_all_resample*.9 + knots_all_resample(k,:)*.1 ...
        + (max(x) - min(x))*randn(size(knots_all_resample))*.01;
    
    % test convergence
    t1 = ttest_mine(knots_all(:,1),knots_prior(:,1));
    t2 = ttest_mine(knots_all(:,1),knots_prior(:,1));
    
    knots_prior = knots_all;
    if abs(t1)<1.732 & abs(t2)<1.732 & J>min_loop
        knots_all = a;
        return
    end
    
    clear knots_all_resample
end







