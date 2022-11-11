% the previous version forgots the pure driange case

function [case_use,P_est,Rall] = DONG_class_fit_All(x,y,LOOP,nc_size)
% x = mm; y = dd; LOOP = 5; nc_size = 300;

% fit 3: ABC
% if max SM is below 0.2, I will think it is very dry and is not an ABC
% form

% added 0923
% k = find( y> std(y)*3);
% jk = x(k)< prctile(x,75);
% y(k(jk)) = []; x(k(jk)) = [];
% 
if max(x)<.2
    P3 = randn(100,6);
    rall3 = 10^6; 
    case_1_flag = 1; % too dry, try case 1
else
    [P3,rall3] = fit_All_scheme3(x,y,LOOP,nc_size);
    
    % Added DONG 0908. Scheme 3 may not converge for some cases
%     if range(P3(:,1))./mean(P3(:,1))>2 |  prctile(P3(:,1),25)<0
%         [P3,rall3] = fit_All_scheme3(x,y,LOOP,nc_size);
%     end
 
    
    % updated July 0715
    k = find( rall3 ==1000000); rall3(k) = nan; 
    a = mean(P3);
    
    if a(5)<prctile(x,5) | a(6)> prctile(x,95) | a(3)<a(1)/10 | a(1)<0 | a(6) - a(5)<.05 | a(3) < 5; % added by DONG 0208
        rall3 = 1000000*ones(size(rall3));
    end
%     if mean(P3(:,6)) - mean(P3(:,5)) <.075; rall3 = 10^6; end
    case_1_flag = 0;
end
% fit case 6
[P6,rall6] = fit_All_scheme6(x,y,LOOP,nc_size);
a = mean(P6);
if a(1)<0 | -a(2)./a(1)<.15; rall6 = 10^6*ones(size(rall6)); end


T_thres = .0;
k = find( isnan(rall3) ); rall3(k) = []; P3(k,:) = [];

% Test if the first slope is significant
Lt = sqrt(length(P3));
a = mean(P3)./std(P3)*Lt;
aslope = a(1);
 
% Test if the second slope is significant
bslope = a(3);


% fit 2: AB
[P2, rall2]= fit_All_scheme2(x,y,LOOP,nc_size);
k = find( rall2>10^5); rall2(k) = []; P2(k,:) = [];
% updated 0715
if mean(P2(:,3))> prctile(x,95)-.02 & max(x)<.25; rall2 = ones(size(rall2))*10^6;end
if prctile( x,5)>.2; rall2 = rall2 + 5;end

if aslope<1.75 & bslope<1.75; case_use = 4; end
n = length(y); k = 1;
P4 = nan; rall4 =n*log(var(y)) + 2*k;% + 2*k*(k+1)/( n - k-1);


% added 0930
if prctile(P3(:,1),5) > 0 & prctile(P2(:,1),5)>0 & prctile(rall2 - rall4,95)<0 & prctile(rall3 - rall4,95)<0
    rall4 = 10^6; 
end




% fit 1: A
[P1, rall1] = fit_All_scheme1(x,y,LOOP,nc_size);
k = find( rall1>10^5); rall1(k) = []; P1(k,:) = [];
% newly added 0917
if min(x)>.25; rall1 = ones(size(rall1))*10^6; end


% fit 5: BC
if max(x)<.25
    rall5 = 10^6;
    P5 = nan;
else
    [P5, rall5]= fit_All_scheme5(x,y,LOOP,nc_size);
    if mean(P5(:,end))< max(prctile(x,5),  range(x)/10); rall5 = 10^6* size(rall5);end
    if mean(P5(:,1))< .5; rall5 = 10^6* size(rall5);end
    
    % added Aug 31, set a 1mm threshold for the case 5 PET
    if mean(P5(:,2))<.5; rall5 = 10^6* size(rall5);end
    if min(x) <.2; rall5 = rall5 + 5; end
end

% AIC score, larger is worse
% merge case 1 and case 6 here
if mean(rall6)< mean( rall1);
    P1 = P6; rall1 = rall6;
end


T21 = ttest_mine(rall2,rall1);  T51 = ttest_mine(rall5,rall1);  T31 = ttest_mine(rall3,rall1);
T32 =  ttest_mine(rall3,rall2); T35 = ttest_mine(rall3,rall5); T52= ttest_mine(rall5,rall2);


a = mean(P1);theta_w = -a(2)./a(1);
if min([T21,T51,T31])>-1.732 & or(theta_w>.15, theta_w<=0.15 & max(x)<.35) % if all complex cases < linear case
    if theta_w>.15; case_use = 6; rall6 = rall1; P6 = P1; end
    if theta_w<=0.15 & max(x)<.35; case_use = 1;end
else
    % second slope is significantly >0, and case 3 is not significantly worse than 32 or 5
    if bslope>1.732 &T32 < T_thres & T35 < T_thres & or( std(P3(:,1))<std(P2(:,1))*2, std(P3(:,1))./mean(P3(:,1)) <.5) % added 0908
        case_use = 3;
    else % slope <0 or T3 is not better than 2 or 5, case 3 is out
        case_use = 4; % case 4 is used as the default
        % if the slope of case 2>0 and Radjust > 0
        %             if prctile(P2(:,1),5)>0 & prctile(rall2,5) >0; case_use = 2; end
        a = mean(P2)./std(P2)*Lt;
        if a(1)>0; case_use = 2; end % the slope of case 2 is significant
        
        % if the slope of case 5>0 and Radjust > 0 & better than case 2
        %             if T25 < -T_thres & prctile(rall5,5) >0 & prctile(P5(:,1),5)>0; case_use = 5; end
        a = mean(P5)./std(P5)*Lt;
        if T52 < T_thres & a(1)>0; case_use = 5; end % if the slope of case 5 is significant and better than case 2
        
    end
end



% updated 0715
if rall4 < min( [mean(rall2),mean(rall3),mean(rall5),mean(rall1)]) & max(x)>.1 & median(y)>.5
    case_use = 4; 
end


eval(['P_est = nanmedian(P',num2str(case_use),',1);'])
eval(['Rall = nanmedian(rall',num2str(case_use),');'])






