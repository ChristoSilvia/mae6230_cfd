
n = 128; 
d = .176;
%construct chebfuns representing eqns
y1 = chebfun(@(x) x.^n-1, [1 1.2]); 
y2 = chebfun(@(x) 1/(d)*(x.^(.7*n)-1), [1 1.2]); 

% Chebfun finds two roots, we ignore r=1. 
r = roots(y1-y2);r = r(2); 

%find delta y1
dy1 = (r-1)./(r^n-1); 


%check results
j = ceil(.7*n); alpha = 1.05; d = 0.176;
err1 = 1- dy1*(r^128-1)/(r-1)
err2 = d- dy1*(r^(.7*128)-1)/(r-1)
dy1 = (r-1)./(r^n-1); 

%%%%%%%%%%
%check

j = ceil(.7*n); 

err1 = 1-dy1*(r^128-1)/(r-1)

err2 = d- dy1*(r^(.7*128)-1)/(r-1)


