 %% =========================================================================
% JIIT B.Tech Major Project – Semester VII (Batch A17)
% Project: Performance Analysis of OTFS under High-Mobility Multipath
%          and Fading Channels
%
% Team Members:
%   Arnav Maindola – 23119021
%   Mukta Pande    – 23119019
% Supervisor:
%   Prof. Alok Joshi
%
% Files:
%   1. otfs_transmitter.m - QPSK modulation, DD grid mapping, ISFFT, Heisenberg
%   2. otfs_channel.m     - EVA multipath, Doppler (0-500 km/h), 5 fading models
%   3. otfs_receiver.m    - MMSE equalizer, QPSK demodulation, BER computation
%   4. main.m             - Master simulation script for Part A & Part B
% =========================================================================

clear; clc; close all;

fprintf('=======================================================================\n');
fprintf('  JIIT Major Project: OTFS Performance Analysis (Batch A17)\n');
fprintf('  Arnav Maindola (23119021) & Mukta Pande (23119019)\n');
fprintf('=======================================================================\n\n');

%% 1. Simulation Parameters
params.M          = 32;          % Subcarriers / delay bins
params.N          = 16;          % Time slots / Doppler bins
params.delta_f    = 15e3;        % Subcarrier spacing = 15 kHz
params.fc         = 4e9;         % Carrier frequency = 4 GHz
params.mod_type   = 'QPSK';      % Modulation
params.detector   = 'MMSE';      % Equalizer

% Number of Monte Carlo frames (increased to 500 for smooth curves)
num_frames = 100;

% Evaluation ranges
snr_range_db       = 0:5:30;     % SNR from 0 to 30 dB
velocity_range_kmh = 0:50:500;   % Velocity from 0 to 500 km/h
fixed_velocity_kmh = 120;        % Fixed speed for BER vs SNR
fixed_snr_db       = 15;         % Fixed SNR for BER vs Velocity (optimal error statistics)

fading_models = {'Rayleigh', 'Rician', 'Nakagami-m', 'alpha-mu', 'kappa-mu'};

fprintf('Simulation Settings: %d frames per point (%d bits/point)\n\n', ...
        num_frames, num_frames * params.M * params.N * 2);

%% =========================================================================
% PART A: HIGH-MOBILITY VEHICULAR PROPAGATION (EVA CHANNEL)
% =========================================================================
fprintf('>>> PART A: OTFS over EVA Multipath Channel...\n');

%% Part A - Plot 1: BER vs SNR (v = 120 km/h)
fprintf('  [Plot 1] Simulating BER vs SNR at v = %d km/h...\n', fixed_velocity_kmh);
params_A1 = params;
params_A1.velocity_kmh = fixed_velocity_kmh;
params_A1.fading_model = 'Rayleigh';

ber_otfs_snr = zeros(size(snr_range_db));

tStart = tic;
for s = 1:length(snr_range_db)
    params_A1.snr_db = snr_range_db(s);
    errors = 0;
    total_bits = 0;

    for f = 1:num_frames
        [tx_sig, ~, tx_bits] = otfs_transmitter(params_A1);
        [rx_grid, H_eff, noise_var] = otfs_channel(tx_sig, params_A1);
        [~, ~, ~, n_err] = otfs_receiver(rx_grid, H_eff, noise_var, tx_bits, params_A1);

        errors = errors + n_err;
        total_bits = total_bits + length(tx_bits);
    end
    ber_otfs_snr(s) = errors / total_bits;
    fprintf('    SNR = %2d dB | BER = %.4e | Errors = %d / %d\n', ...
            params_A1.snr_db, ber_otfs_snr(s), errors, total_bits);
end
fprintf('  Plot 1 completed in %.2f s.\n\n', toc(tStart));

% Note: OFDM comparison curve will be added here in future work.

%% Part A - Plot 2: BER vs Velocity (SNR = 10 dB)
fprintf('  [Plot 2] Simulating BER vs Velocity at SNR = %d dB...\n', fixed_snr_db);
params_A2 = params;
params_A2.snr_db       = fixed_snr_db;
params_A2.fading_model = 'Rayleigh';

ber_otfs_vel = zeros(size(velocity_range_kmh));

tStart = tic;
% In main.m (lines 92-108):

for v = 1:length(velocity_range_kmh)
    params_A2.velocity_kmh = velocity_range_kmh(v);
    errors = 0;
    total_bits = 0;

    for f = 1:num_frames
        % Common Random Numbers (CRN): Synchronize seed per frame f
        rng(2000 + f, 'twister');

        [tx_sig, ~, tx_bits] = otfs_transmitter(params_A2);
        [rx_grid, H_eff, noise_var] = otfs_channel(tx_sig, params_A2);
        [~, ~, ~, n_err] = otfs_receiver(rx_grid, H_eff, noise_var, tx_bits, params_A2);

        errors = errors + n_err;
        total_bits = total_bits + length(tx_bits);
    end
    ber_otfs_vel(v) = errors / total_bits;
    fprintf('    v = %3d km/h | BER = %.4e | Errors = %d / %d\n', ...
            params_A2.velocity_kmh, ber_otfs_vel(v), errors, total_bits);
end
fprintf('  Plot 2 completed in %.2f s.\n\n', toc(tStart));

% Note: OFDM Doppler degradation curve will be added here in future work.


%% =========================================================================
% PART B: OTFS UNDER 5 STATISTICAL FADING MODELS
% =========================================================================
fprintf('>>> PART B: OTFS under 5 Fading Models (Rayleigh, Rician, Nakagami, alpha-mu, kappa-mu)...\n');

%% Part B - Plot 3: BER vs SNR across 5 Fading Models (v = 120 km/h)
fprintf('  [Plot 3] Simulating BER vs SNR for 5 Fading Models...\n');
snr_b3 = 0:5:25; % Evaluated up to 25 dB
ber_fading_snr = zeros(5, length(snr_b3));

% Rayleigh from Part A can be directly mapped to save computation time
[common_snr, idx_A, idx_B] = intersect(snr_range_db, snr_b3);
ber_fading_snr(1, idx_B) = ber_otfs_snr(idx_A);

tStart = tic;
for m = 2:5
    current_model = fading_models{m};
    fprintf('    -> Model: %-12s\n', current_model);

    params_B = params;
    params_B.velocity_kmh = fixed_velocity_kmh;
    params_B.fading_model = current_model;

    for s = 1:length(snr_b3)
        params_B.snr_db = snr_b3(s);
        errors = 0;
        total_bits = 0;

        for f = 1:num_frames
            [tx_sig, ~, tx_bits] = otfs_transmitter(params_B);
            [rx_grid, H_eff, noise_var] = otfs_channel(tx_sig, params_B);
            [~, ~, ~, n_err] = otfs_receiver(rx_grid, H_eff, noise_var, tx_bits, params_B);

            errors = errors + n_err;
            total_bits = total_bits + length(tx_bits);
        end
        ber_fading_snr(m, s) = errors / total_bits;
    end
end
fprintf('  Plot 3 completed in %.2f s.\n\n', toc(tStart));

%% Part B - Plot 4: BER vs Velocity across 5 Fading Models (SNR = 10 dB)
fprintf('  [Plot 4] Simulating BER vs Velocity for 5 Fading Models...\n');
vel_b4 = 0:50:500;
ber_fading_vel = zeros(5, length(vel_b4));

% Rayleigh from Part A mapped directly
% In main.m (lines 160-188):
tStart = tic;
for m = 1:5
    current_model = fading_models{m};
    fprintf('    -> Model: %-12s\n', current_model);

    params_B = params;
    params_B.snr_db       = fixed_snr_db;
    params_B.fading_model = current_model;

    for v = 1:length(vel_b4)
        params_B.velocity_kmh = vel_b4(v);
        errors = 0;
        total_bits = 0;

        for f = 1:num_frames
            % Same frame seed ensures identical bits, scatterer angles, and noise
            rng(4000 + f, 'twister');

            [tx_sig, ~, tx_bits] = otfs_transmitter(params_B);
            [rx_grid, H_eff, noise_var] = otfs_channel(tx_sig, params_B);
            [~, ~, ~, n_err] = otfs_receiver(rx_grid, H_eff, noise_var, tx_bits, params_B);

            errors = errors + n_err;
            total_bits = total_bits + length(tx_bits);
        end
        ber_fading_vel(m, v) = errors / total_bits;
    end
end
fprintf('  Plot 4 completed in %.2f s.\n\n', toc(tStart));


%% =========================================================================
% PLOTTING (4 Core Result Plots as per Synopsis Section 7)
% =========================================================================
fprintf('>>> Generating Plots...\n');

colors = [
    0.00, 0.45, 0.74;  % Blue
    0.85, 0.33, 0.10;  % Orange/Red
    0.93, 0.69, 0.13;  % Yellow/Amber
    0.49, 0.18, 0.56;  % Purple
    0.47, 0.67, 0.19   % Green
];
markers = {'o-', 's--', 'd-.', '^:', 'v-'};

%% --- Plot 1: Part A - BER vs SNR ---
fig1 = figure('Name', 'Part A - BER vs SNR', 'Color', 'w');
semilogy(snr_range_db, max(ber_otfs_snr, 1e-6), 'bo-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'MarkerFaceColor', [0.3 0.6 1.0], 'DisplayName', sprintf('OTFS MMSE (v = %d km/h)', fixed_velocity_kmh));
grid on; grid minor;
xlim([0, 30]); ylim([1e-5, 1]);
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Part A: OTFS Performance over EVA Channel (v = %d km/h)', fixed_velocity_kmh), ...
      'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
saveas(fig1, 'plot1_partA_ber_vs_snr.png');

%% --- Plot 2: Part A - BER vs Velocity ---
fig2 = figure('Name', 'Part A - BER vs Velocity', 'Color', 'w');
semilogy(velocity_range_kmh, ber_otfs_vel, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'MarkerFaceColor', [1.0 0.4 0.4], 'DisplayName', sprintf('OTFS MMSE (SNR = %d dB)', fixed_snr_db));
grid on; grid minor;
xlim([0, 500]); ylim([1e-3, 1e-1]);
xlabel('Vehicle Velocity (km/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Part A: OTFS Mobility Resilience over EVA Channel (SNR = %d dB)', fixed_snr_db), ...
      'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
saveas(fig2, 'plot2_partA_ber_vs_velocity.png');

%% --- Plot 3: Part B - BER vs SNR under 5 Fading Models ---
fig3 = figure('Name', 'Part B - BER vs SNR (5 Fading Models)', 'Color', 'w');
hold on;
model_names_plot = {'Rayleigh', 'Rician (K = 6 dB)', 'Nakagami-m (m = 2)', '\alpha-\mu (\alpha=1.2, \mu=0.8)', '\kappa-\mu (\kappa=3, \mu=2)'};
for m = 1:5
    semilogy(snr_b3, max(ber_fading_snr(m, :), 1e-6), markers{m}, ...
             'Color', colors(m, :), 'LineWidth', 2, 'MarkerSize', 8, ...
             'DisplayName', model_names_plot{m});
end
grid on; grid minor;
xlim([0, 25]); ylim([1e-5, 1]);
xlabel('SNR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Part B: OTFS under 5 Fading Models (EVA Profile, v = %d km/h)', fixed_velocity_kmh), ...
      'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
saveas(fig3, 'plot3_partB_fading_vs_snr.png');

%% --- Plot 4: Part B - BER vs Velocity under 5 Fading Models ---
fig4 = figure('Name', 'Part B - BER vs Velocity (5 Fading Models)', 'Color', 'w');
hold on;
for m = 1:5
    semilogy(vel_b4, ber_fading_vel(m, :), markers{m}, ...
             'Color', colors(m, :), 'LineWidth', 2, 'MarkerSize', 8, ...
             'DisplayName', model_names_plot{m});
end
grid on; grid minor;
xlim([0, 500]); ylim([1e-3, 1e-1]);
xlabel('Vehicle Velocity (km/h)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Bit Error Rate (BER)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Part B: OTFS across Fading Models vs Velocity (SNR = %d dB)', fixed_snr_db), ...
      'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'southwest', 'FontSize', 11);
saveas(fig4, 'plot4_partB_fading_vs_velocity.png');

%% =========================================================================
% SUMMARY TABLES PRINTED TO COMMAND WINDOW
% =========================================================================
fprintf('\n=======================================================================\n');
fprintf('                   SIMULATION RESULTS SUMMARY                          \n');
fprintf('=======================================================================\n');
fprintf('Part A: OTFS BER vs SNR (v = %d km/h)\n', fixed_velocity_kmh);  
disp(table(snr_range_db', ber_otfs_snr', 'VariableNames', {'SNR_dB', 'OTFS_BER'}));

fprintf('\nPart A: OTFS BER vs Velocity (SNR = %d dB)\n', fixed_snr_db);
disp(table(velocity_range_kmh', ber_otfs_vel', 'VariableNames', {'Velocity_kmh', 'OTFS_BER'}));

fprintf('\nPart B: BER vs SNR across 5 Fading Models\n');
T_snr = array2table(ber_fading_snr', 'VariableNames', matlab.lang.makeValidName(fading_models));
T_snr = [table(snr_b3', 'VariableNames', {'SNR_dB'}), T_snr];
disp(T_snr);

fprintf('\nPart B: BER vs Velocity across 5 Fading Models (SNR = %d dB)\n', fixed_snr_db);
T_vel = array2table(ber_fading_vel', 'VariableNames', matlab.lang.makeValidName(fading_models));
T_vel = [table(vel_b4', 'VariableNames', {'Velocity_kmh'}), T_vel];
disp(T_vel);

fprintf('\nAll 4 figures saved successfully as PNG graphics in working directory:\n');
fprintf('  1. plot1_partA_ber_vs_snr.png\n');
fprintf('  2. plot2_partA_ber_vs_velocity.png\n');
fprintf('  3. plot3_partB_fading_vs_snr.png\n');
fprintf('  4. plot4_partB_fading_vs_velocity.png\n');
fprintf('=======================================================================\n');
