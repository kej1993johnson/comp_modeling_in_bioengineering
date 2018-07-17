function [ nmodel] = simmodeld(d, t, N0)

nmodel = N0*exp(-d*t);


end

