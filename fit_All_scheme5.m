function  [P, rall]= fit_All_scheme5(x,y,LOOP,nc_size)

x = x(~isnan(x)); y = y(~isnan(y));

tup = min( prctile(x,90), 0.37);
knots_all =  rand(nc_size*2,1) * (tup - prctile(x,10)) +prctile(x,10);

knots_prior = knots_all;

for J = 1:LOOP
    rall = nan(length(knots_all),1);
    P = nan(length(knots_all),3);
    rmsd = rall;
    
    %     tic
    for i = 1:length(knots_all)
        try
            [rall(i),P(i,:),rmsd(i)] = scheme_All_DONG5(x,y,knots_all(i,:));
        end
    end
    %     toc
    
    k = find(isnan(rall) | isnan(P(:,1))) ;
    rall(k) = []; P(k,:) = []; knots_all(k,:) = []; rmsd(k) = [];
    k = find(P(:,end)<prctile(x,5) | P(:,end)> prctile(x,95) | P(:,1)<1 ); rmsd(k) = 10^6;
    
    % calculate likelihood
    % modified DONG Aug 31
%     k = find(rmsd<10^5); rm = min(rmsd(k));
%     lik = exp( -rmsd);     lik = max(lik,10^-99);
%     
%     
%     if exp(-rm)<10^-50; sz =100; lik = exp(-rmsd/sz); lik = max(lik,10^-99); end
%     k = find(knots_all>.4); lik(k) = 0;
%     wei = lik./sum(lik);
    k = find(isnan(rmsd)); rmsd(k) = 10^6;
    wei = likli_adjust(rmsd);

    
%     lik = exp( -rmsd);     lik = max(lik,10^-99);
%     
%     k = find(knots_all>.4); lik(k) = 0;
%     wei = lik./sum(lik);
%     if std(wei) <.001; 
%         lik = exp(-rmsd/3); lik = max(lik,10^-99); 
%         k = find(knots_all>.4); lik(k) = 0;
%         wei = lik./sum(lik);
%     end
    
    
    % resample
    cwei = cumsum(wei); L = length(cwei);
    u = linspace( 0.001,.999,L);
    clear ki knots_all_resample
    for i = 1:L
        k = find( cwei> u(i),1);
        knots_all_resample(i,:) = knots_all(k,:);
    end
    k = randsample(L,L,'true');
    a = knots_all;
    knots_all = knots_all_resample*.9 + knots_all_resample(k,:)*.1 ...
        + (max(x) - min(x))*randn(size(knots_all_resample))*.01;
    
    % test convergence
    t = ttest_mine(knots_all,knots_prior);
    knots_prior = knots_all;
    if abs(t)<1.732 & J>2
        knots_all = a;
        return
    end
    
    clear knots_all_resample
end
