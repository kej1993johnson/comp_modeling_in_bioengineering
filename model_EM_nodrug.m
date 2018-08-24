function [ Cc] = model_EM_nodrug(Ce_init, p, tsamp)
ge = p(1); % growth rate of epithelial cells
gm = p(2); % growth rate of mesenchymal cells
% in Maclean paper, this was varied from 0 to 0.01 { 0, 0.0001, 0.001,
% 0.0025, 0.01}
kem = p(3); % baseline transition rate E to M 
kme = p(4); % transition rate from M to E



f = @(t,Cc) [ge*Cc(1)- kem.*Cc(1) + kme*Cc(2);  % dE/dt
             gm*Cc(2)+ kem*Cc(1) - kme*Cc(2)]; % dM/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:2);
[t,Cc]=ode45(f, tsamp,Ce_init, options);

E= Cc(:,1);
M = Cc(:,2);

end