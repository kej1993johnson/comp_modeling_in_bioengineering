function[err] = fit_EMT(pfit, pgiven, tsamp, Nmeas, N0)

% p = [ ge; gm; kem; kme; de; dm; fracE];
kem = pfit(1);
de = pfit(2);
dm = pfit(3);
fracE= pfit(4);

ge = pgiven(1);
gm = pgiven(2);
kme = pgiven(3);

p = vertcat(pgiven(1:2), pfit(1), pgiven(3), pfit(2:4));
[Nmodel, E, M]=EMT_model(p, N0, tsamp);

err =Nmeas-Nmodel;
end
