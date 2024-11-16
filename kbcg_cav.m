function [cav,tau,phi,sigma] = kbcg_cav(M,Rrup,vs30,ZTOR,mechanism,region)
% CAV models using KBCG functional form

% region
% 'Alaska': 1, 
% 'Cascadia': 2, 
% 'CentralAmerica&Mexico': 3, 
% 'Japan': 4,
% 'NewZealand': 5, 
% 'SouthAmerica': 6, 
% 'Taiwan': 7


% input : M, R,vs30,Ztor,mechanism,region
% return CAV in m/s

p = readtable('./map.csv');
p{:,:}(isnan(p{:,:})) = 0;

k1 = 400;
k2 = -1.311;
c = 1.88;
n = 1.18;

theta10 = 0;
ZB_if = 30;
ZB_slab =80;
vsflag1 = vs30<=k1;
vsflag2 = vs30>k1;

switch region
    
    case '1_Alaska'
        MB_slab = 7.2;
        MB_if  = 8.4;
        reg = 1;
    case '2_Cascadia'
        MB_slab = 6.4;
        MB_if  = 7.7;
        reg = 2;
    case '3_CentralAmerica&Mexico'
        MB_slab = 6.7;
        MB_if  = 7.3;
        reg = 3;
    case '4_Japan'
        MB_slab = 7.6;
        MB_if  = 8.1;
        reg = 4;
    case '5_NewZealand'
        MB_slab = 7.5;
        MB_if  = 8.3;
        reg = 5;
    case '6_SouthAmerica'
        MB_slab = 7;
        MB_if  = 8.5;   
        reg = 6;
    case '7_Taiwan'
        MB_slab = 7.1;
        MB_if  = 7.1;
        reg = 7;
        
    case '0_global'
        MB_slab = 7.6;
        MB_if  = 7.9;
        reg = 100;
    otherwise
        disp('warning! invalid region string')
end


        
switch mechanism
    case 'interface'
        Fs = 0;
    case 'intraslab'
        Fs = 1;
end      


% regionalized coeff
theta1_if = mean(p.mu_theta1_if(1) + p.delta_theta1_if(reg));
theta1_slab = mean( p.mu_theta1_slab(1) + p.delta_theta1_slab(reg) );
theta7 = mean( p.mu_theta7(1) + p.delta_theta7(reg) );
theta_62 = mean( p.mu_theta6_B(1) + p.delta_theta62(reg) );


%% global adjustment
%if reg == 100
% theta1_if = mean(p.mu_theta1_if(1) + p.delta_theta1_if([4,6,7]));
% theta1_slab = mean( p.mu_theta1_slab(1) + p.delta_theta1_slab([4,6,7]) );
% theta7 = mean( p.mu_theta7(1) + p.delta_theta7([4,6,7]) );
%theta_62 = -0.0088;
    
    
%end


    
    

%% d  f
% global coeff
theta5 = p.theta5(1);
theta2_if = p.theta2_if(1);
theta2_slab = p.theta2_slab(1);
theta3 = p.theta3(1);
theta4_if = p.theta4_if(1);
theta4_slab = p.theta4_slab(1);
theta_nft1 = p.theta_nft1(1);
theta_nft2 = p.theta_nft2(1);
dZB_if = p.dZB_if(1);
dZB_slab = p.dZB_slab(1);
theta9_if = p.theta9_if(1);
theta9_slab = p.theta9_slab(1);


for i = 1:length(M)
    for j =1:length(Rrup)

intercept  = (1-Fs)*theta1_if + Fs*theta1_slab;

f_mag = (1-Fs)*lh(M(i),MB_if,theta4_if*(MB_if - 6),theta4_if,theta5,0.1)...
+ Fs*lh(M(i),MB_slab,theta4_slab *(MB_slab - 6),theta4_slab,theta5,0.1);


f_geom = (1-Fs)*(theta2_if + theta3*M(i)).* log(Rrup(j) + 10.^(theta_nft1 +theta_nft2*(M(i)-6)))...
+ Fs*(theta2_slab + theta3*M(i)).*log(Rrup(j) + 10.^(theta_nft1 +theta_nft2*(M(i)-6)));

f_depth = (1-Fs)*lh(ZTOR,ZB_if + dZB_if, theta9_if*(ZB_if +dZB_if -15),theta9_if,theta10,1)+...
Fs*lh(ZTOR,ZB_slab + dZB_slab, theta9_slab*(ZB_slab +dZB_slab -50),theta9_slab,theta10,1);


f_attn = theta_62*Rrup(j);

siterock = (theta7 + k2*n)*log(1100/k1);
PGArock  = exp(intercept + f_mag + f_geom + f_depth  + f_attn );

f_site = vsflag1* ( theta7 * log(vs30/k1) + k2* (log(PGArock + c*(vs30/k1)^n) - log(PGArock +c) ))...
+ vsflag2*(theta7 +k2*n)*log(vs30/k1);

y = intercept + f_mag + f_geom + f_depth + f_site + f_attn;

cav(j,i) = exp(y)*9.81;


tau(j,i) = p.betsigma(1);
phi(j,i) = p.sd_y(1);
sigma(j,i) = sqrt(tau(j,i)^2 + phi(j,i)^2);
    end
end





function res = lh (x,x0,a,b0,b1,delta)
 res = a + b0 * (x -x0) + (b1 - b0) * delta * log(1+ exp( (x-x0)/delta )  );
end





end