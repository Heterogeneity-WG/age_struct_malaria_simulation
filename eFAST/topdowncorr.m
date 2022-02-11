%% Calculate of top-down correlation coefficient
% R[POIs, N]: matrix of ranks; (data replace by their ranks - from large to small, tiedrank in MATLAB)
% n = # of POIs; k = # of different sample sizes N
% S: Savage scores, rT: TDCC, sign: p-value
% based on Iman, R.L. and Conover W. J., A Measure of Top-Down
% Correlation. Technometrics, August 1987, vol. 29 (3)
function [S, rT, sign]=topdowncorr(R)
[n, k]=size(R);
S=zeros(n,k);
for s=1:k % for each column
    [R_temp,x]=sortrows(R,s);  % sort rows
    for i=1:n
        %        i
        for j=i:n
            S(i,s);
            S(i,s)=S(i,s)+1/R_temp(j,s); % Savage Scores
        end
    end
    S(x,s)=S(:,s);
end
if k<3 % only two columns/different sample size N
    [rT]=corr(S);
    rT=rT(2) % the coefficient ( 2-by-2 matrix, (2,1) entry -> corr between column 1 and 2)
    p=normpdf(sqrt(n-1)*rT);
    sign=p
else % for k>2
    a1=sum(S,2);
    S1=mean(max(S));
    a2=sum(a1.*a1);
    b=k;
    CT=((a2-n*b^2)/(b^2*(n-S1)))
    p=chi2pdf(b*(n-1)*CT,n-1);
    sign=p
end

