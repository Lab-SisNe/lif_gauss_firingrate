function F = root2d(in)

ve =in(1);
vi =in(2);
g  = 4.0; %in(3);
gi = 4.0;
tau_i = 20;

tau_e=20;
Ce=1000; gama=0.25;
H=10; teta=20; tau_0=2;
J=0.1;

uext = 30;
ue = uext + tau_e*Ce*J*(ve - vi*g*gama);
ui = uext + tau_i*Ce*J*(ve - vi*gi*gama);

sigma_e = sqrt( Ce*(J^2) * tau_e * (ve + gama*g^2*vi) );
sigma_i = sqrt( Ce*(J^2) * tau_i * (ve + gama*gi^2*vi) );

fi = @(x) exp(x.^2).*(1+erf(x));

%---------------%
func = @(u,sigma,tau,H,teta,tau_0) (tau_0 + sqrt(pi)*tau*quadgk(fi,(H-u)/sigma,(teta-u)/sigma))^(-1);


F(1) = ve - func(ue,sigma_e,tau_e,H,teta,tau_0);
F(2) = vi - func(ui,sigma_i,tau_i,H,teta,tau_0);

end