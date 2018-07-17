function [ Nmodel,E, M ] = EMT_model( p, N0, tsamp)
% p = [ ge; gm; kem; kme; de; dm; fracE];

ge = p(1);
gm = p(2);
kem = p(3);
kme = p(4);
de = p(5);
dm = p(6);
fracE= p(7);

C_init(1) = N0*fracE;
C_init(2) = N0*(1-fracE);

f = @(t,C) [ ge.*C(1) - kem.*C(1) + kme.*C(2)- de.*C(1); % dE/dt
             gm.*C(2) - kme.*C(2) + kem.*C(1)- dm.*C(2)]; % dM/dt

options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:2);
[t,C]=ode45(f, tsamp,C_init, options);
E=C(:,1);
M = C(:,2);
Nmodel = E+M;

end

