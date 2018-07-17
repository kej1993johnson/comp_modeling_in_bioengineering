function [ nmodel] = simmodelg(g, t, N0)

nmodel = N0*exp(g*t);


end
