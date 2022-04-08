% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; R1U = 1/15;
kPL = 0.05; 
std_noise = 0.01;

input_condition = 1; % choose from various simulated starting conditions
switch input_condition
    case 1
        % gamma variate input function - most realistic
        Tarrival = 0;  Tbolus = 12;
        input_function = realistic_input_function(N, TR, Tarrival, Tbolus);
        Mz0 = [0,0,0]; 
    case 2
        % boxcar input function
        Tbolus = 12;  Tarrival = 0;
        Ibolus = [1:round(Tbolus/TR)] + round(Tarrival/TR);
        Rinj = 1/Tbolus; % normalize for a total magnetization input = 1
        Mz0 = [0,0,0]; input_function(Ibolus) =  Rinj*TR;
    case 3
        Mz0 = [1.5,0,1]; % no input function
    case 4
        Tin = 6; Mz0 = Tin; % no input function, delayed start
        input_function = [];
end

flips = repmat([20;35;20]*pi/180,[1 N]); 

t = [0:N-1]*TR + Tin;

% generate simulated data
noise_S = randn([3 N])*std_noise;  % same noise for all flip schedules
UreaPyr_ratio = rand(1)*2;
UreaPyr_ratio = 1/2;

[Mxy_PL, Mz_PL] = simulate_Nsite_model(Mz0(1:2), [R1P R1L], [kPL 0], flips(1:2,:), TR, input_function);
[Mxy_U, Mz_U] = simulate_Nsite_model(Mz0(3), [R1U], [], flips(3,:), TR, input_function);


Mxy = [Mxy_PL; Mxy_U * UreaPyr_ratio];
Mz = [Mz_PL; Mz_U * UreaPyr_ratio];
Sn = Mxy + noise_S;

figure
subplot(311), plot(t, flips*180/pi)
title('flip angles')
legend('pyruvate', 'lactate', 'urea')
subplot(312), plot(t, Mxy, t, Sn)
title('M_{XY} and signal')
legend('pyruvate', 'lactate', 'urea')
subplot(313), plot(t, Mz)
title('M_Z')
legend('pyruvate', 'lactate', 'urea')



% initial parameter guesses
R1P_est = 1/25; R1L_est = 1/25; R1U_est = 1/15;
kPL_est = .02;
plot_fits = 1;
fit_function = @fit_pyr_kinetics_Ktrans_Extended;

%% Test fitting - 2-site
disp('2-site model: pyruvate -> lactate')
disp('Fitting kPL with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit 
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1U = R1U_est;
params_est.kPL = kPL_est; 


% no noise - this should be perfect
    [params_fit Sfit] = fit_pyr_kinetics_Ktrans_Extended(Mxy, TR, flips, params_fixed, params_est, [], plot_fits);
pause
    [params_fit_nourea Sfit_nourea] = fit_pyr_kinetics(Mxy(1:2,:), TR, flips(1:2,:), params_fixed, params_est, [], plot_fits);
pause
    
    % add noise
    [params_fitn Snfit ] = fit_pyr_kinetics_Ktrans_Extended(Sn, TR, flips, params_fixed, params_est, [], plot_fits);
pause
    [params_fitn_nourea Snfit_nourea ] = fit_pyr_kinetics(Sn(1:2,:), TR, flips(1:2,:), params_fixed, params_est, [], plot_fits);
pause

disp('Press any key to continue')
disp(' ')

return

%% Test fitting - 2-site plus T1 lactate fitting
disp('2-site model: pyruvate -> lactate')
disp('Fitting kPL  and T1 of lactate')
disp('Fitting both parameters leads to increases in variability, which can be')
disp('alleviated to some extent by constraints on the parameter values')
disp('(Pyruvate T1 values have very small effect on inputless fitting approach)')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est;
params_est.kPL = kPL_est;  params_est.R1L = R1L_est;
% set constraints on lactate T1:
params_est.R1L_lb = 1/40;
params_est.R1L_ub = 1/15;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1,1:size(Mxy,2),  Iflips)] = fit_function(Mxy(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1,1:size(Mxy,2),  Iflips)] = fit_function(Sn(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1,1:size(Mxy,2),  Iflips)] = fit_function(abs(Sn(1:2,:,Iflips)), TR, flips(1:2,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL      R1L '])
disp(num2str([kPL R1L],2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      R1L '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      R1L '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      R1L '])
disp([num2str(k_fit,2), flip_description_array])

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)

disp('Press any key to continue')
disp(' ')

pause

%% Test fitting - fit kPX only
disp('4-site model: pyruvate -> lactate, bicarbonate, alanine')
disp('Fitting kPL, kPB, and kPA, with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:3,1:size(Mxy,2),  Iflips)] = fit_function(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:3,1:size(Mxy,2),  Iflips)] = fit_function(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
 %   [params_fitn_mag(:,Iflips) Snfit_mag(1:3,1:size(Mxy,2),  Iflips)] = fit_function(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL      KPB      KPA    '])
disp(num2str([kPL, kPB, kPA],2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:3];[1:3] + Nparams_fit;[1:3] + 2*Nparams_fit]); 
disp(['KPL      KPB      KPA    '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:3];[1:3] + Nparams_fit;[1:3] + 2*Nparams_fit]); 
disp(['KPL      KPB      KPA    '])
disp([num2str(k_fit,2), flip_description_array])
% disp('Noisy magnitude fit results:')
% k_fit = struct2array(params_fitn_mag);
% k_fit = k_fit([[1:3];[1:3] + Nparams_fit;[1:3] + 2*Nparams_fit]); 
% disp(['KPL      KPB      KPA    '])
% disp([num2str(k_fit,2), flip_description_array])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
%plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits')% (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)
subplot(223) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
%plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')
subplot(224) , plot(t, squeeze(Sn(4,:,:)))
hold on, plot(t, squeeze(Snfit_complex(3,:,:)),':')
%plot(t, squeeze(Snfit_mag(3,:,:)),'--'), hold off
title('Alanine signals and fits')

disp('Press any key to continue')
disp(' ')

pause
%% Test fitting - 3-site
disp('3-site model: pyruvate -> lactate, bicarbonate')
disp('Fitting kPL, and kPB, with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:2,1:size(Mxy,2),  Iflips)] = fit_function(Mxy(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:2,1:size(Mxy,2),  Iflips)] = fit_function(Sn(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:2,1:size(Mxy,2),  Iflips)] = fit_function(abs(Sn(1:3,:,Iflips)), TR, flips(1:3,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL      KPB   '])
disp(num2str([kPL, kPB],2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      KPB   '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      KPB   '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      KPB   '])
disp([num2str(k_fit,2), flip_description_array])

figure
subplot(131) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(132) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)
subplot(133) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')

return
disp('Press any key to continue')
disp(' ')

pause

