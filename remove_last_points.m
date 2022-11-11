function [mm,dd] = remove_last_points(mm,dd)

% if the SM within the 75 - 100th percentile range is lower than -3*std of
% the 50-75 percentile mean, remove it

thres = -1.8;

rates = prctile(mm,[25,50,75]);

d1 = dd( mm <= rates(1));
d2 = dd( mm > rates(1) & mm<= rates(2) );
d3 = dd( mm > rates(2) & mm<= rates(3) );
d4 = dd( mm >  rates(3) );

t1 = ones(size( d1));

dd2 = d2 - mean(d1);
t2 = dd2/std(d1)*sqrt( length(d1));

dd3 = d3 - mean(d2);
t3 = dd3/std(d2)*sqrt( length(d1));

dd4 = d4 - mean(d3);
t4 = dd4/std(d3)*sqrt( length(d1));

tall = [t1; t2; t3; t4];

%% remove upper outliers

dd2 = d1 - mean(d2);
t1 = dd2/std(d2)*sqrt( length(d1));

dd3 = d2 - mean(d3);
t2 = dd3/std(d3)*sqrt( length(d1));

dd4 = d3 - mean(d4);
t3 = dd4/std(d4)*sqrt( length(d1));

t4 = ones(size(d4));

tall2 = [t1; t2; t3; t4];




%%
k = find(tall <thres | tall2>-thres | dd < range(mm)*.01);

mm(k) = []; dd(k) = [];



