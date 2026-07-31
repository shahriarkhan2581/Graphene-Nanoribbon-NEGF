% all grn iv curve.m
% Generates IV curves for a set of graphene nanoribbon (GNR) widths using
% a simple, tunable toy-model transmission (sum of Lorentzians) and the
% Landauer formula. This script is intended to produce reproducible
% IV-characteristics for multiple ribbon widths so you can visualize and
% compare trends. Replace the toy-model `compute_transmission` with a
% real NEGF transmission routine from this repository if available.
%
% Output:
% - Figure saved to `gnr_iv_curves.png`
% - Data saved to `gnr_iv_curves.mat` and `gnr_iv_curves.csv`
%
% Usage:
% - Run from MATLAB: `run('all grn iv curve.m')`
% - Edit `widths` and `biases` variables below to change the sweep.

clearvars; close all; clc;

% Physical constants
e = 1.602176634e-19;   % elementary charge, C
h = 6.62607015e-34;    % Planck constant, J*s
kB = 1.380649e-23;     % Boltzmann constant, J/K

% Simulation parameters
tempK = 300;                      % temperature in K
kT_eV = kB*tempK / e;             % kT in eV for FD functions

widths = [4, 8, 12, 16];          % ribbon "widths" (arbitrary units)
length_unit = 50;                 % length (affects resonant density)

biases = linspace(0, 0.6, 61);    % bias sweep in Volts (0 to 0.6 V)
energy_grid = linspace(-2.5, 2.5, 4001);  % energy grid in eV for integration

% Transmission toy-model parameters (tunable)
gamma0 = 0.05;    % broadening (eV) base
E_spacing0 = 0.2; % base spacing between resonances (eV)

results = struct();
results.widths = widths;
results.biases = biases;
results.energy_grid = energy_grid;
results.tempK = tempK;

% Precompute Fermi function helper
fermi = @(E, mu) 1 ./ (1 + exp((E - mu) ./ kT_eV));

G0 = 2 * e^2 / h; % conductance quantum (S)

for iw = 1:numel(widths)
    W = widths(iw);

    % Build a simple effective-channel model where larger widths -> more
    % channels and denser resonances. This is a phenomenological model.
    Nchannels = max(1, round(W/2));
    gamma = gamma0 * (1 + 0.3*(W-4)/max(1,W));      % slight broadening with width
    E_spacing = E_spacing0 / (1 + 0.05*(W-4));     % smaller spacing for larger W
    center_offset = 0;                              % place resonances near Fermi

    % Create resonance energies for the toy transmission
    En = center_offset + ((1:Nchannels) - (Nchannels+1)/2) * E_spacing;

    % Transmission function as sum of Lorentzians normalized to at most 1
    T_of_E = zeros(size(energy_grid));
    for n = 1:numel(En)
        T_of_E = T_of_E + (gamma^2) ./ ((energy_grid - En(n)).^2 + gamma^2);
    end
    % normalize so peak transmission <= Nchannels (so per-channel max ~1)
    T_of_E = T_of_E ./ max(T_of_E) * min(Nchannels, 1.0*Nchannels);
    % optionally cap transmission
    T_of_E(T_of_E > Nchannels) = Nchannels;

    I_vs_V = zeros(size(biases));

    for ib = 1:numel(biases)
        V = biases(ib);
        muL =  V/2;   % left chemical potential in eV
        muR = -V/2;   % right chemical potential in eV

        fL = fermi(energy_grid, muL);
        fR = fermi(energy_grid, muR);

        integrand = T_of_E .* (fL - fR);
        % Integrate over energy in eV. Landauer: I = (2e^2/h) * 
        % 
        I = G0 * trapz(energy_grid, integrand); % Ampere (since energy in eV and prefactor includes e)
        I_vs_V(ib) = I;
    end

    results(iw).width = W;
    results(iw).Nchannels = Nchannels;
    results(iw).energy = energy_grid;
    results(iw).T = T_of_E;
    results(iw).I = I_vs_V;
end

% Plot IV curves
figure('Color','w','Position',[200 200 700 500]); hold on; box on;
colors = lines(numel(widths));
for iw = 1:numel(widths)
    plot(biases, results(iw).I*1e6, '-o', 'Color', colors(mod(iw-1,size(colors,1))+1,:), 'LineWidth',1.5, 'MarkerSize',4);
end
xlabel('Bias (V)');
ylabel('Current (\muA)');
title('GNR IV curves (toy-model transmission)');
legend(arrayfun(@(w) sprintf('W=%d',w), widths, 'UniformOutput', false), 'Location','northwest');
grid on;

% Save outputs
saveas(gcf, 'gnr_iv_curves.png');
save('gnr_iv_curves.mat','results');

% Also export CSV with columns: bias, I(W1), I(W2), ... (in A)
csv_header = ['bias_V', arrayfun(@(w) sprintf('I_W%d_A',w), widths, 'UniformOutput', false)];
csv_data = zeros(numel(biases), numel(csv_header));
csv_data(:,1) = biases(:);
for iw = 1:numel(widths)
    csv_data(:, 1+iw) = results(iw).I(:);
end

% Write CSV with header
fid = fopen('gnr_iv_curves.csv', 'w');
fprintf(fid, '%s', csv_header{1});
for c = 2:numel(csv_header)
    fprintf(fid, ',%s', csv_header{c});
end
fprintf(fid, '\n');
fclose(fid);
% Append numeric data
dlmwrite('gnr_iv_curves.csv', csv_data, '-append');

fprintf('Done: generated IV curves for widths: %s\n', mat2str(widths));
fprintf('Saved: gnr_iv_curves.png, gnr_iv_curves.mat, gnr_iv_curves.csv\n');

% End of script
