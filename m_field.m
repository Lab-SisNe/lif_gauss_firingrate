%----------------------------------------------------------------------%
% This code solves the equation for average firing-rate of a LIF neuron
% driven by white Gaussian noise as given by Brunel (2000) (and others).
%
% Assumptions: Network has sparse random connectivity (C << N )
% Thus, correlations of the fluctuating part of the synaptic inputs 
% of different neurons are neglected in the limit C/N -> 0.
%
% How it solves: fsolve tries to find the roots from Eq.(21) from Brunel
% (2000) iteratively. The conditions such as number of iterations or step 
% are set by optimset.
%----------------------------------------------------------------------%
% by: rodrigo pena
% pena@njit.edu / rfdop20@gmail.com                             
%----------------------------------------------------------------------%

clear

%---------------------------------------------%
%       Initial guess for firing-rate E/I
%---------------------------------------------%
x0 = [0.3, 0.3];
%---------------------------------------------%

%---------------------------------------------%
%   Defining iteration rules and running
%---------------------------------------------%
options=optimset('Display','iter','LargeScale','off','TolFun',.0001,'MaxIter',100000,'MaxFunEvals',10000);
[x,fval] = fsolve(@equations,x0,options);
%---------------------------------------------%

%---------------------------------------------%
%                  Results
%---------------------------------------------%
['firing-rate E = ' num2str(x(1)*1000) 'Hz'] 
['firing-rate I = ' num2str(x(2)*1000) 'Hz']


function F = equations(in)

    ve =in(1);  %Firing-rate E from last iteration
    vi =in(2);  %Firing-rate I from last iteration
    g  = 4.0;   %relative inhibition I->E
    gi = 4.0;   %relative inhibition I->I
    tau_e=20;   %membrane time constante E
    tau_i = 20; %membrane time constante I
    Ce=1000;    %number of excitatory connections
    gama=0.25;  %relative number of inhibitory connections (Ci=gama*Ce)
    H=10;       %voltage reset
    teta=20;    %voltage threshold
    tau_0=2;    %refractory period in [ms]
    J=0.1;      %synaptic efficacy

    %---------------------------------------------%
    %           Eq.(4) - Brunel (2000)
    %---------------------------------------------%
    uext = 30;
    ue = uext + tau_e*Ce*J*(ve - vi*g*gama);
    ui = uext + tau_i*Ce*J*(ve - vi*gi*gama);
    %---------------------------------------------%
    

    %---------------------------------------------%
    %           Eq.(5) - Brunel (2000)
    %---------------------------------------------%
    sigma_e = sqrt( Ce*(J^2) * tau_e * (ve + gama*g^2*vi) );
    sigma_i = sqrt( Ce*(J^2) * tau_i * (ve + gama*gi^2*vi) );
    %---------------------------------------------%
    
    
    %---------------------------------------------%
    %      function that will be integrated
    % note that integral is solved using quadgk
    % (adaptive Gauss-Kronrod quadrature method)
    %---------------------------------------------%
    fi = @(x) exp(x.^2).*(1+erf(x));
    %---------------------------------------------%
    

    %---------------------------------------------%
    %           Eq.(21) - Brunel (2000)
    %---------------------------------------------%
    func = @(u,sigma,tau,H,teta,tau_0) (tau_0 + sqrt(pi)*tau*quadgk(fi,(H-u)/sigma,(teta-u)/sigma))^(-1);
    %---------------------------------------------%
    
    
    %---------------------------------------------%
    %         Output for both populations
    %---------------------------------------------%
    F(1) = ve - func(ue,sigma_e,tau_e,H,teta,tau_0);
    F(2) = vi - func(ui,sigma_i,tau_i,H,teta,tau_0);

end
