% functions to solve

clear

x0 = [0.3 0.3];
options=optimset('Display','iter','LargeScale','off','TolFun',.0001,'MaxIter',100000,'MaxFunEvals',10000);

[x,fval] = fsolve(@root2d,x0,options);
x(1)*1000 
x(2)*1000 

