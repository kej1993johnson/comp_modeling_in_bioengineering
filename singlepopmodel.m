function model = singlepopmodel( params0, dose)
n = length(dose);

% Computes elongated matrix of max viability to be multiplied by each
% component of the sigmoidal curve 


        


    %LD50res, sloperes, LD50sens, slopesens, fres
    LD50 = params0(1);
    slope = params0(2);

    model = (1./( 1 + exp(slope.*(dose - LD50))));
