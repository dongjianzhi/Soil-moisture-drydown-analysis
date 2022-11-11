function  [P, rall]= fit_All_scheme2(x,y,LOOP,nc_size)

x = x(~isnan(x)); y = y(~isnan(y));
knots_all = rand(nc_size*2,1) * (max(x) - prctile(x,5)) +prctile(x,5);

tlow = 0.0; tup = .25; % tup = min( 0.25, prctile(x,10));
theta_w = rand(nc_size*2,1) * (tup - tlow) + tlow;

tup = min( prctile(x,90), 0.37);
theta_fc =  rand(nc_size*2,1) * (tup - prctile(x,10)) +prctile(x,10);

knots_all = [theta_w, theta_fc];
knots_prior = knots_all;

clear knots_all_resample

for J = 1:LOOP
    rall = nan(length(knots_all),1);
    P = nan(length(knots_all),3);
    rmsd = rall;
    
    %     tic
    for i = 1:length(knots_all)
        try
            [rall(i),P(i,:),rmsd(i)] = scheme_All_DONG2(x,y,knots_all(i,:));
        end
    end
    %     toc
    
    k = find(isnan(rall) | isnan(P(:,1))) ;
    rall(k) = []; P(k,:) = []; knots_all(k,:) = []; rmsd(k) = [];
    
    k =find(P(:,1)<0); rmsd(k) = 10^6;
    k = find(P(:,3)> prctile(x,95) | P(:,3)< prctile(x,5)); rmsd(k) = nan;
    
    % modified DONG Aug 31
%     k = find(rmsd<10^5); rm = min(rmsd(k));
%     lik = exp( -rmsd);   
%     
%     if exp(-rm)<10^-75; sz = rm/75;  lik = exp(-rmsd/sz); end
%     wei = lik./sum(lik);
    k = find(isnan(rmsd)); rmsd(k) = 10^6;
    wei = likli_adjust(rmsd);

    % calculate likelihood
%     
%     
%     k = find(knots_all(:,2)>.4); lik(k) = 0;
%     wei = lik./sum(lik);
%     if std(wei) <.001; 
%         lik = exp(-rmsd/3); lik = max(lik,10^-99); 
%         k = find(knots_all(:,2)>.4); lik(k) = 0;
%         wei = lik./sum(lik);
%     end
    
   
    
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
    if abs(t)<1.732 & J>3
        knots_all = knot_rec;
        return
    end

    clear knots_all_resample
end  
