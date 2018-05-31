function model = twopopmodel2( params0, fres1, dose, nsamp )
LD50res = params0(1);
    sloperes = params0(2);
    LD50sens = params0(3);
    slopesens = params0(4);

    
    model = (1.*((fres1)./( 1 + exp(sloperes.*(dose - LD50res))) + ((1-fres1)./(1 + exp(slopesens.*(dose - LD50sens))))));
end