% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; R1U = 1/15;
kPL = 0.05; 

%% Testing parameters here identify two potential areas of value:

test_scenario = 1;

switch test_scenario
    
    case 1 % "normal"
        std_noise = 0.02; UreaPyr_ratio = 1;  % comparable SNRs
        flips = repmat([20;35;20]*pi/180,[1 N]);  % flip anlges for [pyruvate; lactate; urea]
    case 2
        % if pyruvate SNR is low, but urea is high, this seems to improve the kPL
        % estimates
        std_noise = 0.04; UreaPyr_ratio = 4;  % low SNR pyruvate, high SNR urea
        %std_noise = 0.02; UreaPyr_ratio = 1;  % comparable SNRs
        %std_noise = 0.02; UreaPyr_ratio = 1/4;  % high SNR pyruvate, low SNR urea
        flips = repmat([20;35;20]*pi/180,[1 N]);  % flip anlges for [pyruvate; lactate; urea]
        
    case 3
        % can potentially use very small or no flip angle on pyruvate, and still
        % estimate kPL
        flips = repmat([1;35;20]*pi/180,[1 N]); % can get away without sampling pyruvate
        std_noise = 0.015; UreaPyr_ratio = 1;
        
end

%%

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
        Mz0 = [1.5,0,1.5*UreaPyr_ratio]; % no input function
        input_function = [];

end

t = [0:N-1]*TR + Tin;

% generate simulated data
noise_S = randn([3 N])*std_noise;  % same noise for all flip schedules

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

disp('Press any key to continue')
disp(' ')
pause

return


%% Test fitting - 2-site plus T1 lactate fitting
% Not sure there is any benefit here
disp('2-site model: pyruvate -> lactate')
disp('Fitting kPL  and T1 of lactate')
disp('Fitting both parameters leads to increases in variability, which can be')
disp('alleviated to some extent by constraints on the parameter values')
disp('(Pyruvate/Urea T1 values have very small effect on inputless fitting approach)')
disp('')

clear params_fixed params_est
params_fixed.R1P = R1P_est; params_fixed.R1U = R1U_est;
params_est.R1L = R1L_est; 
params_est.kPL = kPL_est;

% set constraints on lactate T1:
params_est.R1L_lb = 1/40;
params_est.R1L_ub = 1/15;


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

