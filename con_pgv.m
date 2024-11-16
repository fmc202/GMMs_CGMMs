function [PGV,PGV_sigma] = con_pgv(M,R,Vs30,Sa,Sa_sigma,mechanism)
% PGV in cm/s
% relation of Tpgv and M
% crustal: T = exp(-4.09 + 0.66*M)
% this code is consisitent with the functional form in peer report
% which is the most updated

% this function if fully vectorized
switch mechanism   
    case 'crustal'
        a1 = 5.39;
        a2 = 0.799;
        a3 = 0.654;
        a4 = 0.479;
        a5 = -0.062;
        a6 = -0.359;
        a7 = -0.134;
        a8 = 0.023; 
        tau = 0.16;
        phi = 0.29;
        stdd = 0.33;
end




% magnitude term
m_term1 = a2;
m_term2 = a2 + (a3 - a2)*(M - 5)/2.5;
m_term3 = a3;
fm = (M < 5).*m_term1 + (M >= 5).*(M <= 7.5).*m_term2 + (M > 7.5).*m_term3; 

% median model
lnPGV = a1 + fm.*log(Sa) + a4*(M - 6) + a5*(8.5 - M).^2 + a6*log( R + 5*exp(0.4*(M - 6)) )+...
        (a7 + a8*(M - 5)) .* log(Vs30/425);
PGV = exp(lnPGV);

% prop of error
PGV_sigma = sqrt( stdd^2 + (fm.^2).*(Sa_sigma.^2) );

end

