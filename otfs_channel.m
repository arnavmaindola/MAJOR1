function [rx_grid_dd, H_eff, noise_var] = otfs_channel(tx_sig, params)
% OTFS Channel
% Simulates EVA multipath, vehicle Doppler, statistical fading, and AWGN.
%
% Inputs:
%   tx_sig - Transmitted signal vector (MN x 1)
%   params - Struct with M, N, delta_f, fc, velocity_kmh, snr_db, fading_model
%
% Outputs:
%   rx_grid_dd - Received Delay-Doppler grid (M x N)
%   H_eff      - Effective Delay-Doppler channel matrix (MN x MN)
%   noise_var  - Noise variance sigma^2

    M = params.M;
    N = params.N;
    MN = M * N;
    delta_f = params.delta_f;
    fc = params.fc;
    v = params.velocity_kmh;
    snr_db = params.snr_db;
    fading = params.fading_model;

    Ts = 1 / (M * delta_f);

    % Extended Vehicular A (EVA) delay profile
    tau_sec = [0, 30, 150, 310, 700, 1090, 1730, 2510, 17900] * 1e-9;
    p_db = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
    p_lin = 10.^(p_db / 10);
    p_lin = p_lin / sum(p_lin);
    P = length(p_lin);

    % Doppler shift: nu = nu_max * cos(theta)
    nu_max = ((v / 3.6) * fc) / 3e8;
    theta = 2 * pi * rand(P, 1) - pi;
    nu = nu_max * cos(theta);

    % Fading envelope generation for the 5 models
    switch lower(fading)
        case 'rayleigh'
            % Non-Line-of-Sight Rayleigh fading
            h_gain = sqrt(p_lin / 2)' .* (randn(P, 1) + 1i * randn(P, 1));

        case 'rician'
            % Line-of-Sight on path 1 (K = 6 dB), scattered on remaining paths
            K_lin = 10^(6 / 10);
            h_gain = sqrt(p_lin / 2)' .* (randn(P, 1) + 1i * randn(P, 1));
            los = sqrt(K_lin / (K_lin + 1));
            nlos = sqrt(1 / (2 * (K_lin + 1))) * (randn() + 1i * randn());
            h_gain(1) = sqrt(p_lin(1)) * (los + nlos);

        case 'nakagami-m'
            % Nakagami-m fading (m = 2: light fading / multi-cluster)
            m = 2.0;
            R = sqrt(local_gamma(m, p_lin / m));
            h_gain = R' .* exp(1i * 2 * pi * rand(P, 1));

        case 'alpha-mu'
            % alpha-mu fading (alpha = 1.2, mu = 0.8: severe non-linear fading)
            alpha = 1.2; mu = 0.8;
            r_hat = sqrt(p_lin .* (mu^(2 / alpha)) * gamma(mu) / gamma(mu + 2 / alpha));
            Y = local_gamma(mu, 1 / mu, P);
            R = r_hat .* (Y.^(1 / alpha));
            h_gain = R' .* exp(1i * 2 * pi * rand(P, 1));

        case 'kappa-mu'
            % kappa-mu fading (kappa = 3.0, mu = 2.0: dominant specular clusters)
            kappa = 3.0; mu = 2.0;
            d_sq = kappa * p_lin / (1 + kappa);
            sigma_sq = p_lin / (2 * mu * (1 + kappa));
            G_comp = sqrt(d_sq / 4) + sqrt(sigma_sq) .* randn(4, P);
            R = sqrt(sum(G_comp.^2, 1));
            h_gain = R' .* exp(1i * 2 * pi * rand(P, 1));

        otherwise
            error('Unknown fading model: %s', fading);
    end

    % Precalculated sinc interpolation weights and transformation matrix
    persistent sinc_mat U
    if isempty(sinc_mat)
        max_delay = min(ceil(max(tau_sec) / Ts) + 3, M - 1);
        L = max_delay + 1;
        sinc_mat = zeros(L, P);
        for l = 0:max_delay
            for i = 1:P
                val = (l * Ts - tau_sec(i)) / Ts;
                if val == 0
                    sinc_mat(l+1, i) = 1;
                else
                    sinc_mat(l+1, i) = sin(pi * val) / (pi * val);
                end
            end
        end
        FN = dftmtx(N) / sqrt(N);
        U = kron(FN, eye(M));
    end

    % Time-domain channel matrix G
    t_samp = (0:MN-1)' * Ts;
    D = exp(1i * 2 * pi * t_samp * nu'); % MN x P
    G = zeros(MN, MN);
    max_delay = size(sinc_mat, 1) - 1;
    row_idx = (1:MN)';
    for l = 0:max_delay
        gl = (D .* h_gain.') * sinc_mat(l+1, :)';
        col_idx = mod(row_idx - 1 - l, MN) + 1;
        G(sub2ind([MN, MN], row_idx, col_idx)) = gl;
    end

    % Effective Delay-Doppler channel matrix
    H_eff = U * G * U';

    % Convert transmitted time signal to Delay-Doppler domain: x_dd = U * tx_sig
    x_dd = U * tx_sig(:);

    % Additive White Gaussian Noise
    noise_var = 10^(-snr_db / 10);
    noise_dd = sqrt(noise_var / 2) * (randn(MN, 1) + 1i * randn(MN, 1));

    % Received Delay-Doppler signal
    y_dd = H_eff * x_dd + noise_dd;
    rx_grid_dd = reshape(y_dd, M, N);
end

% Self-contained Gamma random generator (Marsaglia-Tsang method)
function g = local_gamma(a, b, P)
    if nargin < 3, P = length(b); end
    if isscalar(b), b = repmat(b, 1, P); end
    g = zeros(1, P);
    for k = 1:P
        g(k) = b(k) * gamma_single(a);
    end
end

function val = gamma_single(a)
    if a < 1
        val = gamma_single(a + 1) * (rand()^(1 / a));
        return;
    end
    d = a - 1 / 3;
    c = 1 / sqrt(9 * d);
    while true
        z = randn();
        v = (1 + c * z)^3;
        if v > 0
            u = rand();
            if u < 1 - 0.0331 * (z^4)
                val = d * v;
                return;
            end
            if log(u) < 0.5 * (z^2) + d - d * v + d * log(v)
                val = d * v;
                return;
            end
        end
    end
end
