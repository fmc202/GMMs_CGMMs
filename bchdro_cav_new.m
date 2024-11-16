function [cav,tau,phi,sigs] = bchdro_cav_new(M,Rrup,Vs30,ZTOR,mechanism,region)

% input : M, R,vs30,Ztor,mechanism,region----- all  column vector
% mechanism must be a column of 0 or 1
% mechansim = 1 for intraslab
% return CAV in m/s


switch mechanism
case 'interface'
    Fs = 0;
case 'intraslab'
    Fs = 1;
end   

 
Vlin = 400;
n  = 1.18;
b = -1.186;
c  = 1.88;
vssmall = Vs30<=Vlin;
vsbig = Vs30>Vlin;
Vs30_new = min(Vs30,1000);


[MB_slab,MB_face,regid] = recode_reg(region);   



C1 = Fs.*MB_slab + (1-Fs).*MB_face;
C1slab_C1face = MB_slab - MB_face;
Zsmall = ZTOR <= 100;
Zbig = 1 - Zsmall;


%% regionalized coefficients
a1 = zeros([8,1]);
a6 = zeros([8,1]);
a12 = zeros([8,1]);



a1(1) = 5.8296580672347025;
a1(2) = 5.417233486580535;
a1(3) = 5.490012590600625;
a1(4) = 6.533105272959904;
a1(5) = 6.010065946790736;
a1(6) = 5.987430888220086;
a1(7) = 5.8209731208585405;
a1(8) = 5.862000290023683;


a6(1) = -0.0018302197902628334;
a6(2) = -0.0017838052139808507;
a6(3) = -0.000653888537631098;
a6(4) = -0.003426219378279003;
a6(5) = -0.00026332949359995727;
a6(6) = -0.0013159887027458614;
a6(7) = -0.002142987310449599;
a6(8) = mean(a6(1:7));



a12(1) = 1.5261529578118391;
a12(2) = 1.1604407953128222;
a12(3) = 1.8085335610174127;
a12(4) = 0.9802893533021287;
a12(5) = 1.5106176902369122;
a12(6) = 2.1425799075174594;
a12(7) = 1.3009089912966283;
a12(8) = 1.492901608848976;
 
%% global coeff
a4  =  0.7681282143627861;
a14 = -0.10108804263877584;
a2  = -1.178415320010986;
a3  = 0.03284090916243754;
a9  = 0.08558363591082635;
C4  = 10;
a10 = 1.0728532795044654;
a13 =  -0.07493843250944249;
a5  = 0.49738057877883707;
a11 = 0.005288538731072816;
phi1 =  0.5968277386018089;
tau1 = 0.38099561590236064;

for i = 1:length(M)
    for j =1:length(Rrup)

%% prediction

  Msmall = M(i) <= C1;
  Mbig = M(i) > C1;
  intercept = a1(regid) + a4*C1slab_C1face .*Fs;

  f_geom =  (a2 + a14*Fs + a3*(M(i)-7.8)) .* log(Rrup(j) + C4*exp( (M(i)-6)*a9 ) );
  
  f_attn = a6(regid) .*Rrup(j) + a10*Fs;

  f_mag = Msmall.*(a4*(M(i)-C1) + a13*(10-M(i)).^2 ) + Mbig.*(a5*(M(i)-C1) + a13*(10-M(i)).^2 );
    
  f_depth = Zsmall*a11 .*(ZTOR-60) .*Fs + Zbig*a11*(100-60) .*Fs ;
 
  siterock = (a12(regid) + b*n)*log(1100/Vlin);
  IMrock  = exp(intercept + f_mag + f_geom + f_depth  + f_attn  + siterock);
    

  f_site = vssmall .* a12(regid) .*log(Vs30_new/Vlin) - b*log(IMrock +c) +...
   b*log(IMrock + c*(Vs30_new/Vlin).^n)   + vsbig .*(a12(regid) +b*n) .*log(Vs30_new/Vlin);   


  lncav = intercept + f_mag + f_geom + f_depth + f_site + f_attn;
  
    cav(j,i) = exp(lncav)*9.81;
    tau(j,i) = tau1;
    phi(j,i) = phi1;
    sigs(j,i) = sqrt(tau(j,i)^2 + phi(j,i)^2);

    end
end










%% functions
function [MB_slab,MB_face,reg] = recode_reg(reg_str)

switch reg_str 
    case '1_Alaska'
        MB_slab = 7.5;
        MB_face  = 7.8;
        reg = 1;
    case '2_Cascadia'
        MB_slab = 7.2;
        MB_face  = 8.2;
        reg = 2;
    case '3_CentralAmerica&Mexico'
        MB_slab = 7.5;
        MB_face  = 7.8;
        reg = 3;
    case '4_Japan'
        MB_slab = 7.5;
        MB_face  = 7.8;
        reg = 4;
    case '5_NewZealand'
        MB_slab = 7.5;
        MB_face  = 7.8;
        reg = 5;
    case '6_SouthAmerica'
        MB_slab = 7.5;
        MB_face  = 7.8;   
        reg = 6;
    case '7_Taiwan'
        MB_slab = 7.5;
        MB_face  = 7.8;
        reg = 7;
        
    case '0_global'
        MB_slab = 7.5;
        MB_face  = 7.8;
        reg = 8;
    otherwise
        disp('warning! invalid region string')
end
end

end