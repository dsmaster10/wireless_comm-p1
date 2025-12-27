%% NOMA vs OMA BER Comparison with Data Rate Analysis
%              Demonstrates NOMA's superior spectral efficiency and higher data rates
%              compared to OMA in a two-user scenario
% Demonstrates proper SIC implementation for NOMA superiority

%clear all; close all; clc;

%% System Parameters
clc; clear; close all;
num_users = 2;              % Number of users
num_bits = 1e5;             % Number of bits for reliable BER
SNR_dB = 0:2:30;           % SNR range in dB
num_iter = 200;             % Monte Carlo iterations

% Power allocation (adjusted for better SIC)
alpha1 = 0.3;               % Strong user power coefficient
alpha2 = 0.7;               % Weak user power coefficient

% Channel conditions (increased disparity)
h_mean = [2.0, 0.1];        % Mean channel gains [strong, weak]

% Data rate parameters
bandwidth = 1e6;            % 1 MHz bandwidth
symbol_time = 1e-6;         % 1 μs symbol duration

%% Initialize metrics
BER_NOMA = zeros(num_users, length(SNR_dB));
BER_OMA = zeros(num_users, length(SNR_dB));
data_rate_NOMA = zeros(num_users, length(SNR_dB));
data_rate_OMA = zeros(num_users, length(SNR_dB));

%% Main Simulation Loop
for snr_idx = 1:length(SNR_dB)
    SNR_linear = 10^(SNR_dB(snr_idx)/10);
    noise_var = 1/SNR_linear;  % Noise variance
    
    noma_err = zeros(1, num_users);
    oma_err = zeros(1, num_users);
    noma_success = zeros(1, num_users);
    oma_success = zeros(1, num_users);
    
    for iter = 1:num_iter
        %% Data and Modulation (BPSK)
        data = randi([0 1], num_users, num_bits);
        mod_signal = 2*data - 1;  % BPSK modulation
        
        %% Channel Realization
        %% Channel Realization
        h = h_mean' .* (randn(num_users,1) + 1i*randn(num_users,1)) / sqrt(2);
        
        %% ===== NOMA Transmission =====
        tx_noma = sqrt(alpha1)*mod_signal(1,:) + sqrt(alpha2)*mod_signal(2,:);
        
        % Received signals with noise
        noise = sqrt(noise_var/2)*(randn(num_users,num_bits) + 1i*randn(num_users,num_bits));
        rx_noma = h.*tx_noma + noise;
        
        %% NOMA Reception - Enhanced MMSE-SIC
        % Strong user (User 1) processing
        rx_user1 = rx_noma(1,:);
        
        % MMSE detection for weak user's signal (User 2)
        W = (sqrt(alpha2)*h(1)') / (alpha2*abs(h(1))^2 + alpha1*abs(h(1))^2 + noise_var);
        dec_user2_at_user1 = real(W * rx_user1) > 0;
        
        % Reconstruct and subtract User 2's signal
        reconstructed = sqrt(alpha2)*(2*dec_user2_at_user1 - 1);
        residual_signal = (rx_user1 - h(1)*reconstructed)/h(1);
        
        % Decode User 1's signal
        dec_user1 = real(residual_signal/sqrt(alpha1)) > 0;
        
        % Weak user (User 2) processing
        snr_eff = (alpha2*abs(h(2))^2)/(alpha1*abs(h(2))^2 + noise_var);
        dec_user2 = real(rx_noma(2,:)/(h(2)*sqrt(alpha2))) > 0;
        
        % Count errors
        noma_err = noma_err + [sum(dec_user1 ~= data(1,:)), sum(dec_user2 ~= data(2,:))];
        noma_success = noma_success + [sum(dec_user1 == data(1,:)), sum(dec_user2 == data(2,:))];
        
        %% ===== OMA Transmission =====
        half_len = floor(num_bits/2);
        tx_oma = [mod_signal(1,1:half_len), mod_signal(2,half_len+1:end)];
        
        % Received signals with noise
        noise_oma = sqrt(noise_var/2)*(randn(1,num_bits) + 1i*randn(1,num_bits));
        rx_oma = [h(1)*tx_oma(1:half_len), h(2)*tx_oma(half_len+1:end)] + noise_oma;
        
        % OMA Decoding
        dec_oma_user1 = real(rx_oma(1:half_len)/h(1)) > 0;
        dec_oma_user2 = real(rx_oma(half_len+1:end)/h(2)) > 0;
        
        oma_err = oma_err + [sum(dec_oma_user1 ~= data(1,1:half_len)), ...
                            sum(dec_oma_user2 ~= data(2,half_len+1:end))];
        oma_success = oma_success + [sum(dec_oma_user1 == data(1,1:half_len)), ...
                                   sum(dec_oma_user2 == data(2,half_len+1:end))];
    end
    
    %% Calculate Performance Metrics
    BER_NOMA(:,snr_idx) = noma_err'/(num_iter*num_bits);
    BER_OMA(1,snr_idx) = oma_err(1)/(num_iter*half_len);
    BER_OMA(2,snr_idx) = oma_err(2)/(num_iter*(num_bits-half_len));
    
    % Data rates in Mbps
    data_rate_NOMA(:,snr_idx) = bandwidth * noma_success'/(num_iter*num_bits) / 1e6;
    data_rate_OMA(1,snr_idx) = (bandwidth/2) * oma_success(1)/(num_iter*half_len) / 1e6;
    data_rate_OMA(2,snr_idx) = (bandwidth/2) * oma_success(2)/(num_iter*(num_bits-half_len)) / 1e6;
end

%% Plot Results
figure;

% BER Comparison
subplot(2,1,1);
semilogy(SNR_dB, BER_NOMA(1,:), 'b-o', 'LineWidth', 2, 'DisplayName', 'NOMA User 1 (Strong)');
hold on;
semilogy(SNR_dB, BER_NOMA(2,:), 'r-s', 'LineWidth', 2, 'DisplayName', 'NOMA User 2 (Weak)');
semilogy(SNR_dB, BER_OMA(1,:), 'b--x', 'LineWidth', 2, 'DisplayName', 'OMA User 1');
semilogy(SNR_dB, BER_OMA(2,:), 'r--+', 'LineWidth', 2, 'DisplayName', 'OMA User 2');
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance: NOMA vs OMA');
legend('Location', 'best');
ylim([1e-5 1]);

% Data Rate Comparison
subplot(2,1,2);
plot(SNR_dB, data_rate_NOMA(1,:), 'b-o', 'LineWidth', 2, 'DisplayName', 'NOMA User 1');
hold on;
plot(SNR_dB, data_rate_NOMA(2,:), 'r-s', 'LineWidth', 2, 'DisplayName', 'NOMA User 2');
plot(SNR_dB, data_rate_OMA(1,:), 'b--x', 'LineWidth', 2, 'DisplayName', 'OMA User 1');
plot(SNR_dB, data_rate_OMA(2,:), 'r--+', 'LineWidth', 2, 'DisplayName', 'OMA User 2');
plot(SNR_dB, sum(data_rate_NOMA), 'k-*', 'LineWidth', 2, 'DisplayName', 'NOMA Sum Rate');
plot(SNR_dB, sum(data_rate_OMA), 'k--d', 'LineWidth', 2, 'DisplayName', 'OMA Sum Rate');
grid on;
xlabel('SNR (dB)');
ylabel('Data Rate (Mbps)');
title('Data Rate Comparison: NOMA vs OMA');
legend('Location', 'best');

%% Expected Results
fprintf('\n=== Expected Improvements ===\n');
fprintf('1. NOMA User 1 BER < OMA User 1 BER (at high SNR)\n');
fprintf('2. NOMA Sum Rate ≈ 2x OMA Sum Rate\n');
fprintf('3. NOMA User 2 BER > OMA User 2 BER (acceptable trade-off)\n');