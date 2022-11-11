function wei = likli_adjust(rmsd)


% modified DONG Aug 31
k = find(rmsd<10^5); rm = min(rmsd(k));
lik = exp( -rmsd);

% scale_num = 5;
scale_num = 5;

% if exp(-rm)<10^-scale_num; sz = rm/scale_num;  lik = exp(-rmsd/sz); end
% if exp(-rm)<10^-99; sz = rm/scale_num;  lik = exp(-rmsd/sz); end

if exp(-rm)<10^-99; sz = 100;  lik = exp(-rmsd/sz); end

if min(rmsd>10^5); lik = ones(size(lik)); end

wei = lik./sum(lik);