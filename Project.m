clc
clear 
close all

%% Creators
% Authors: Nima Anvari 401411237, Ahmadreza Ghaderi 401413065
% Fields and waves Project

%% Medium Specifiactions

% we have four mediums, then we'll specify each medium's property

% frequency 
frequency = 1e7;

% free space
eps0 = 8.854187817e-12; mu0 = 4*pi*1e-7;


% first medium
eps1 = 0.7*8.854187817e-12; mu1 = 0.8*4*pi*1e-7;
d1 = 0.05;

% second medium
eps2 = 0.4*8.854187817e-12; mu2 = 0.5*4*pi*1e-7;
d2 = 0.02;

% third medium
eps3 = 0.5*8.854187817e-12; mu3 = 0.2*4*pi*1e-7;
d3 = 0.1;

%fourth medium
etha_PEC = 0;

eps = [eps0,eps1,eps2,eps3];
mu = [mu0,mu1,mu2,mu3];
d = [d1,d2,d3];

%intrinstic impedance 
etha = sqrt(mu./eps);

%% Angle handling

% we'll assume incident angle is 45 degree
% Snell's law of reflection: incident and reflection angle are the same
% Snell's law of refraction: incident and reflcetion angle are related by
% refraction indices

theta_inc = 45; % degree
theta_trans = asind(sqrt(mu(2)*eps(2))/sqrt(mu(1)*eps(1))*sind(theta_inc));

 for i = 2:3
     theta_trans(i) = asind(sqrt(mu(i+1)*eps(i+1))/sqrt(mu(i)*eps(i))*sind(theta_trans(i-1)));
 end

 theta_inc(2) = theta_trans(1);
 theta_inc(3) = theta_trans(2);

fprintf('The transmission and incident angle in each medium is:\n');

for k = 1:numel(theta_inc)
    fprintf('Boundary %d: Tran = %4.2f°, Inc = %4.2f°\n', ...
        k, theta_trans(k), theta_inc(k));
end

%% Impedance calculation 

b = zeros(1,4);
 for i= 1:4
     beta(i) = 2*pi*frequency*sqrt(mu(i)*eps(i));
 end

 % impedance seen from medium 4 onward
 Z(3) = etha(4)*(+1i*etha(4)*tan(beta(4)*(d(3)))/etha(4));

 % impedance seen from medium 3 onward
 Z(2) = etha(3)*(Z(3)+1i*etha(3)*tan(beta(3)*(d(2)))/etha(3)+1i*Z(3)*tan(beta(3)*(d(2))));

 % impedance seen from medium 2 onward
 Z(1) = etha(2)*(Z(2)+1i*etha(2)*tan(beta(2)*(d(1)))/etha(2)+1i*Z(2)*tan(beta(2)*(d(1))));
 
 %% Reflection and transmission coeffs


 % if the incident wave is TM or parallel
%  for i=1:3
%      ref_coeff_TM(i)=(Z(i)*cosd(theta_trans(i))-etha(i)*cosd(theta_inc(i)))/ (Z(i)*cosd(theta_trans(i))+etha(i)*cosd(theta_inc(i)));
%  end
% 
%  for i=3
%      tran_coeff_TM(i)=(2*Z(i)*cosd(theta_trans(i)))/(Z(i)*cosd(theta_trans(i))+etha(i)*cosd(theta_inc(i)));
%  end

 % if the incident wave is TE or normal
 for i=1:3
     ref_coeff_TE(i)=(Z(i)*cosd(theta_inc(i))-etha(i)*cosd(theta_trans(i)))/ (Z(i)*cosd(theta_inc(i))+etha(i)*cosd(theta_trans(i)));
 end

 for i=3
     tran_coeff_TE(i)=(2*Z(i)*cosd(theta_inc(i)))/(Z(i)*cosd(theta_inc(i))+etha(i)*cosd(theta_trans(i)));
 end

%  fprintf('the reflection coeffs for TM incident');
%  ref12 = ref_coeff_TM(1);
%  ref23 = ref_coeff_TM(2);
%  ref34 = ref_coeff_TM(3);
% 
%   fprintf('the reflection coeffs for TE incident');
%  ref12 = ref_coeff_TE(1);
%  ref23 = ref_coeff_TE(2);
%  ref34 = ref_coeff_TE(3);


%% Genetic Algorithm to Minimize Reflection at First Medium (TM) for 3 Layers

% Bounds: [eps1, mu1, d1, eps2, mu2, d2, eps3, mu3, d3]
lb = [0.1*8.854e-12, 0.1*4*pi*1e-7, 0.01, ...
      0.1*8.854e-12, 0.1*4*pi*1e-7, 0.01, ...
      0.1*8.854e-12, 0.1*4*pi*1e-7, 0.01];  % lower bounds

ub = [2*8.854e-12, 2*4*pi*1e-7, 0.2, ...
      2*8.854e-12, 2*4*pi*1e-7, 0.2, ...
      2*8.854e-12, 2*4*pi*1e-7, 0.2];      % upper bounds

nvars = 9; % 3 layers, each with eps, mu, d

% GA options
options = optimoptions('ga','Display','iter', ...
    'PopulationSize',200,'MaxGenerations',400);

% Objective function
fitnessFunc = @(x) reflectionFitness3Layers(x, frequency, eps0, mu0, 45);

% Run GA
[x_opt,fval] = ga(fitnessFunc,nvars,[],[],[],[],lb,ub,[],options);

% Display optimized values
fprintf('Optimized parameters:\n');
for i = 1:3
    fprintf('Layer %d: eps = %.3e, mu = %.3e, d = %.3f m\n', ...
        i, x_opt(3*(i-1)+1), x_opt(3*(i-1)+2), x_opt(3*(i-1)+3));
end
fprintf('Reflection magnitude at first boundary = %.3e\n', fval);

%% --- Fitness function ---
function fitness = reflectionFitness3Layers(x, freq, eps0, mu0, theta_inc_deg)
    % x = [eps1, mu1, d1, eps2, mu2, d2, eps3, mu3, d3]

    eps = [eps0, x(1), x(4), x(7)];
    mu  = [mu0,  x(2), x(5), x(8)];
    d   = [x(3), x(6), x(9)];

    eta = sqrt(mu./eps);
    beta = 2*pi*freq*sqrt(mu.*eps);

    % PEC backing
    Z_next = 0; 
    for i = numel(d):-1:1
        Z(i) = eta(i) * (Z_next + 1i*eta(i)*tan(beta(i)*d(i))) / ...
                        (eta(i) + 1i*Z_next*tan(beta(i)*d(i)));
        Z_next = Z(i);
    end

    % Snell's law for first interface (TM)
    theta_trans_deg = asind(sqrt(mu(2)*eps(2))/sqrt(mu(1)*eps(1))*sind(theta_inc_deg));

    % Reflection coefficient at first boundary (TM)
    ref_TM = (Z(1)*cosd(theta_trans_deg) - eta(1)*cosd(theta_inc_deg)) / ...
             (Z(1)*cosd(theta_trans_deg) + eta(1)*cosd(theta_inc_deg));

    fitness = abs(ref_TM);  % minimize magnitude

end




