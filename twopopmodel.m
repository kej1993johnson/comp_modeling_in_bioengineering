function model = twopopmodel( params0, d, nsamp )

ns = length(d); % number of data points per sample
n = length(d)*nsamp;
d2 = repmat(d,nsamp,1);
dose = reshape(d2,n,1);

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
    model = (1.*(fvec1)./( 1 + exp(sloperes.*(dose - LD50res))) + ((1-fvec1)./(1 + exp(slopesens.*(dose - LD50sens)))));

    
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
    mod = (1.*((fvec1)./( 1 + exp(sloperes.*(dose - LD50res))) + ((1-fvec1)./(1 + exp(slopesens.*(dose - LD50sens))))));

end 
    model = reshape(mod, length(d), nsamp);
  
 
end