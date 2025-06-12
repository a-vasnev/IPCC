% v01 copied from .../Research/2021/Confidence interval combination/Code/interval_comb_v17d_ExampleD.m 
%     where it construct pictures for proposition 2 of Magnus and Vasnev (2023) 
%     "On the uncertainty of a combined forecast: The critical role of correlation" in version 17d
% v02 created for revision in January 2025
% v03 is a cleaned version for github with notation matching the paper
%     Magnus and Vasnev (2025) "The role of data and priors in estimating
%     climate sensitivity"

%% clear the memory

clear all;
close all;


%% input data

caseN = 14; % 14 for IPCC5 data, 4 for IPCC6; 
            % doesn't work for other choices
            % the number corresponds to the number of points and used later in the program

if caseN == 14
    % IPCC5 data without the first point
    b_0i = [0.68; 0.72; 0.72; 0.99; 0.97; 0.75; 0.89; 1.17; 0.93; 0.99; 1.09; 1.23; 1.15; 1.02];
    sigma_0i = [0.21; 0.33; 0.44; 0.18; 0.22; 0.52; 0.43; 0.08; 0.45; 0.40; 0.26; 0.17; 0.28; 0.47];
elseif caseN == 4
    % IPCC6 data
    b_0i = [1.22; 1.03; 1.2; 1.05];
    sigma_0i = [0.36; 0.39; 0.61; 0.36];
else
    error("no such case yet")
end


%% create variable needed in Section 5

x = b_0i;  % observations from Section 5
m = caseN; % number of points corresponds to m in Section 5
vones = ones(m,1); % vector of ones
tmp = m / (sigma_0i'*sigma_0i); % normalization constant for tr(V) = m restriction
V_0 = diag(tmp*(sigma_0i).^2);
v = diag(sqrt(V_0));

% create a grid of correlation values to investigate
r_set = [0:0.0001:0.99];
r_count = length(r_set);
% preallocate estimators
mu_hat_set = zeros(r_count,1); % to store \hat{\mu}
sigma2_hat_set = zeros(r_count,1); % to store \hat{\sigma}^2
tau2_hat_set = zeros(r_count,1); % to store \hat{\tau}
i_Vinv_i_set = zeros(r_count,1); % to store the denominator, \imath'V^{-1}\imath, that appears in the estimators; not relevant for the current paper


%% compute the estimators for each correlation value using the function for case E from Magnus and Vasnev (2023) 

for i = 1:r_count
    [mu_hat_set(i,1), sigma2_hat_set(i,1), tau2_hat_set(i,1), i_Vinv_i_set(i,1)] = ...
        mu_hat_sigma_hat(r_set(i),m,V_0,x); % case E
end

% compute mean and s2
fcst_mean =  mean(x);
err_hat = x - mean(x);
s2 = err_hat'*err_hat / m;

% compute limits using function for Proposition 2 from Magnus and Vasnev (2023)
lim_r = 1;
[lim_mu(1,1), lim_sigma2(1,1), lim_tau2(1,1)] = proposition2_lim(x, v); % case E


%% create figure to look at the behaviour of estimators when correlation changes

scrsz = get(0,'ScreenSize'); figure('PaperPositionMode','auto','Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/(2.5) scrsz(4)/(1.4)],'PaperOrientation','portrait'); hold on;
subplot(3,1,1); plot(r_set, mu_hat_set,'LineWidth',2); hold on;  plot(lim_r,lim_mu,'bo'); legend('Case E','limits','Location','southwest'); 
   title('Panel 1: Estimators of \mu'); 
   xlim([-0.3 1]); ylim([0 1.5]); xlabel('\rho');
subplot(3,1,2); plot(r_set, sigma2_hat_set,'LineWidth',2); hold on;   plot(lim_r,lim_sigma2,'bo'); legend('Case E','limits','Location','northwest'); 
   title('Panel 2: Estimators of \sigma^2'); 
   xlim([-0.3 1]); ylim([0 2.5]); xlabel('\rho');
subplot(3,1,3); plot(r_set, tau2_hat_set,'LineWidth',2); hold on; plot(lim_r,lim_tau2,'bo'); legend('Case E','limits','Location','northwest'); 
   title('Panel 3: Estimators of \tau^2'); 
   xlim([-0.3 1]); ylim([0 0.05]); xlabel('\rho');


% find \rho where \sigma_0 = \sigma_2, i.e., prior variance matches the posterior
if caseN == 14
    tmp = (sigma2_hat_set - 0.2809).^2; % sigma_2^2 =  0.2809 for IPCC5
elseif caseN == 4
    tmp = (sigma2_hat_set - 0.0729).^2; % sigma_2^2 =  0.0729 for IPCC6
else
    error("not implemented yet")
end
[M,I] = min(tmp);
r_tmp = r_set(I); % correlation of the threshold case where var of data = var of posterior for Section 5


%% create Table 5 for the paper

% get estimators for correlation values used in Table 5
r_0 = [0; 0.5; 0.6; 0.7; 0.8; r_tmp; 0.9; 0.95; 0.99];
mu_hat_0 = spline(r_set, mu_hat_set, r_0); % corresponds to b_0 in the table
sigma2_hat_0 = spline(r_set, sigma2_hat_set, r_0); % corresponds to \sigma_0 in the table
tau2_hat_0 = spline(r_set, tau2_hat_set, r_0); % corresponds to \tau in the table

% create Table 5 for the paper
my_table6 = [r_0 mu_hat_0 sqrt(sigma2_hat_0) sqrt(tau2_hat_0)];
my_table6 = [my_table6; ...
             lim_r lim_mu sqrt(lim_sigma2) sqrt(lim_tau2)];

% display table with 2 decimals for the paper
my_table6_rounded = round(my_table6, 2);
if caseN == 14
    disp("Table 6: IPCC5");
elseif caseN == 4
    disp("Table 6: IPCC6");
else
    error("not implemented yet")
end
disp(my_table6_rounded);

%% output produced by Matlab R2024b on MacBook and on PC

% Table 6: IPCC5
%       0    1.07    0.27    0.04
%    0.50    1.20    0.32    0.06
%    0.60    1.21    0.35    0.06
%    0.70    1.21    0.40    0.06
%    0.80    1.22    0.49    0.06
%    0.83    1.22    0.53    0.06
%    0.90    1.22    0.68    0.06
%    0.95    1.23    0.96    0.06
%    0.99    1.23    2.13    0.06
%    1.00    1.23     Inf    0.06

% Table 6: IPCC6
%        0    1.11    0.10    0.04
%     0.50    1.10    0.13    0.09
%     0.60    1.09    0.15    0.10
%     0.70    1.09    0.17    0.12
%     0.80    1.07    0.21    0.14
%     0.88    1.06    0.27    0.16
%     0.90    1.05    0.30    0.17
%     0.95    1.03    0.42    0.19
%     0.99    1.01    0.92    0.21
%     1.00    1.00     Inf    0.22

%% create Table 6 for the paper

if caseN == 14
    b_2 = 1.07;
    sigma_2 = 0.53;
elseif caseN == 4
    b_2 = 1.15;
    sigma_2 = 0.27;
else
    error("not implemented yet")
end

sigma_1 = sqrt((sigma2_hat_0(end-2:end)*sigma_2^2)./(sigma2_hat_0(end-2:end)-sigma_2^2));
b_1 = (sigma2_hat_0(end-2:end)*b_2-sigma_2^2*mu_hat_0(end-2:end))./(sigma2_hat_0(end-2:end)-sigma_2^2);

my_table7 = [r_0(end-2:end) b_1 sigma_1];
my_table7 = [my_table7; ...
             lim_r b_2 sigma_2];
my_table7_rounded = round(my_table7, 2);
if caseN == 14
    disp("Table 7: IPCC5 prior");
elseif caseN == 4
    disp("Table 7: IPCC6 prior");
else
    error("not implemented yet")
end
disp(my_table7_rounded);

%% output produced by Matlab R2024b on MacBook

% Table 7: IPCC5 prior
%     0.90    0.83    0.85
%     0.95    1.00    0.64
%     0.99    1.06    0.55
%     1.00    1.07    0.53

% Table 7: IPCC6 prior
%     0.90    1.66    0.67    
%     0.95    1.24    0.36    
%     0.99    1.16    0.28    
%     1.00    1.15    0.27    

%% functions

function [mu_hat, sigma2_hat, tau2_hat, i_Vinv_i] = mu_hat_sigma_hat(r,sizeN,V_0,fcst_point)
% this function computes the estimators for each correlation value r using case E from Magnus and Vasnev (2023) 

    A = ones(sizeN,sizeN)/sizeN;
    alpha = (sizeN - 1) * r +1;
    beta  = 1 - r;
    P = alpha * A + beta * (eye(sizeN) - A);
    P_inv = (1/alpha) * A + (1/beta) * (eye(sizeN) - A);
    P_det = alpha*beta^(sizeN-1);
    V_temp = diag(1./diag(sqrt(V_0))); %V_temp = inv(sqrtm(V_0));

    vones = ones(sizeN,1);
    mu_hat = (vones'*V_temp*P_inv*V_temp*fcst_point) / (vones'*V_temp*P_inv*V_temp*vones);
    err_hat = fcst_point - mu_hat;
    sigma2_hat = err_hat'*V_temp*P_inv*V_temp*err_hat / sizeN; % equivalent to equation (10) in the paper
    
    V_inv = V_temp * (A/alpha + (eye(sizeN) - A)/beta ) * V_temp;
    i_Vinv_i =  (vones' * V_inv * vones); 
    tau2_hat = sigma2_hat / (vones' * V_inv * vones); % equation (12) in the paper
end

function [mu_hat_lim, sigma2_hat_lim, tau2_hat_lim] = proposition2_lim(x, v)
% this function computes limiting values of the estimators when correlation = 1 
% using function for Proposition 2 from Magnus and Vasnev (2023)

   sizeN = length(x);
   Eps = 1.0e-5;
   vones = ones(sizeN,1);
   A = ones(sizeN,sizeN)/sizeN;
   V0_sqrt = diag(v);
   V0_sqrt_inv = diag(1./v);
   V0 = diag(v.^2);
   w = V0_sqrt_inv * vones;
   y = V0_sqrt_inv * x;
   
   % compute covariances
   C_wy = w' * (eye(sizeN) - A) * y / sizeN;
   C_ww = w' * (eye(sizeN) - A) * w / sizeN;
   C_xx = x' * (eye(sizeN) - A) * x / sizeN;
   C_vv = v' * (eye(sizeN) - A) * v / sizeN;
   C_yy = y' * (eye(sizeN) - A) * y / sizeN;
   
   % limits for mu_hat
   if sum(v == vones) == sizeN
       mu_hat_lim = mean(x);
   else
       mu_hat_lim = C_wy / C_ww;
   end
   
   % limits for sigma2_hat
   if (sum(v == vones) == sizeN) && ( sum( (x/mean(x)) == vones ) == sizeN )
       sigma2_hat_lim = 0;
   elseif (sum(v == vones) < sizeN) && (sum((x - [vones v]*([vones v]\x)).^2) < Eps)
       sigma2_hat_lim = C_xx / (sizeN * C_vv);
   else
       sigma2_hat_lim = Inf;
   end
   
   % limits for tau2_hat
   if ( sum( (x/mean(x)) == vones ) == sizeN )
       tau2_hat_lim = 0;
   elseif ((sum((x/mean(x) == vones)) < sizeN) && (sum(v == vones) < sizeN))
       tau2_hat_lim = (C_ww * C_yy - C_wy^2) / (sizeN * C_ww^2);
   else
       tau2_hat_lim = Inf;
   end
   
end