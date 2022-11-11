function [mm,dd] = rm_outlier(mm,dd)

m1 = mean(dd); s1 = std(dd);

thres = m1 + s1*3;
k = find( dd>thres & mm< prctile(mm,75)); 

mm(k) = []; dd(k) = [];
end