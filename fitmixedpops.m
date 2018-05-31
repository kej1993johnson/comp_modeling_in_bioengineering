function diff = fitmixedpops( params0, dose, viability, nsamp )
n = length(viability);
ns = n/nsamp; % number of data points per sample

switch nsamp
    case 4
        


    %LD50res, sloperes, LD50sens, slopesens, fres
    LD50res = params0(1);
    sloperes = params0(2);
    LD50sens = params0(3);
    slopesens = params0(4);
    fres1 = params0(5);
    fres2 = params0(6);
    fres3 = params0(7);
    fres4 = params0(8);


    fvec1 = zeros([n 1]);
    fvec1(1:ns) = fres1;
    fvec1(ns+1:2*ns) = fres2;
    fvec1(2*ns+1:3*ns) = fres3;
    fvec1(3*ns+1:4*ns) = fres4;
    diff = (1.*(fvec1)./( 1 + exp(sloperes.*(dose - LD50res))) + ((1-fvec1)./(1 + exp(slopesens.*(dose - LD50sens))))) - viability;

    
    case 2
    %LD50res, sloperes, LD50sens, slopesens, fres
    LD50res = params0(1);
    sloperes = params0(2);
    LD50sens = params0(3);
    slopesens = params0(4);
    fres1 = params0(5);
    fres2 = params0(6);
    


    fvec1 = zeros([n 1]);
    fvec1(1:ns) = fres1;
    fvec1(ns+1:2*ns) = fres2;
    diff = (1.*((fvec1)./( 1 + exp(sloperes.*(dose - LD50res))) + ((1-fvec1)./(1 + exp(slopesens.*(dose - LD50sens)))))) - viability;

end 
  
 
end