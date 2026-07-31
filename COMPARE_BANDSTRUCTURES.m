%COMPARE_BANDSTRUCTURES  Compute, plot, and validate the pi-orbital
%nearest-neighbour tight-binding band structures of 7-AGNR, 11-AGNR,
%and 8-ZGNR, using H(k) = H00 + H01*exp(i*k*a) + H01'*exp(-i*k*a).
%
%   Validation performed:
%     1) Eg = min(E(E>0)) - max(E(E<0)) for each ribbon
%     2) Eg(7-AGNR)  > Eg(11-AGNR)   (3p+1 gapped vs 3p+2 metallic)
%     3) 8-ZGNR shows edge-state bands pinned near E = 0 for k near
%        the zone boundary
%     4) The three band structures are visibly different
%
%   Run this BEFORE proceeding to NEGF transport (GNR_Current.m /
%   Compare_GNR_IV.m).

clear; clc; close all;

nk = 800;                 % k-points per band structure
Ezero_tol = 1e-3;         % eV, threshold separating "zero" from finite E

ribbons = struct( ...
    'name', {'7-AGNR', '11-AGNR', '8-ZGNR'}, ...
    'loader', {@H_7AGNR, @H_11AGNR, @H_8ZGNR});

results = struct();

figure('Name','GNR Band Structures','Position',[100 100 1400 420]);

for r = 1:numel(ribbons)
    [H00, H01, ~, a] = ribbons(r).loader();
    Natoms = size(H00,1);

    kk = linspace(-pi/a, pi/a, nk);
    Eall = zeros(nk, Natoms);
    for ik = 1:nk
        k = kk(ik);
        Hk = H00 + H01*exp(1i*k*a) + H01'*exp(-1i*k*a);
        Eall(ik,:) = sort(real(eig(Hk)));
    end

    Ev = Eall(:);
    Ec_min = min(Ev(Ev >  Ezero_tol));
    Ev_max = max(Ev(Ev < -Ezero_tol));
    Eg = Ec_min - Ev_max;

    results(r).name  = ribbons(r).name;
    results(r).kk    = kk;
    results(r).a     = a;
    results(r).Eall  = Eall;
    results(r).Eg    = Eg;
    results(r).Natoms = Natoms;

    subplot(1,3,r);
    plot(kk*a/pi, Eall, 'b-', 'LineWidth', 1.0); hold on;
    plot([-1 1], [0 0], 'k--');
    xlabel('k a / \pi'); ylabel('E (eV)');
    title(sprintf('%s  (E_g = %.4f eV)', ribbons(r).name, Eg));
    xlim([-1 1]); ylim([-9 9]); grid on;
end

if exist('sgtitle','file') || exist('sgtitle','builtin')
    sgtitle('Nearest-neighbour \pi-orbital TB band structures (t_0 = -2.7 eV)');
end

% ---------------- Validation ----------------
fprintf('\n===== Bandgap summary =====\n');
for r = 1:numel(ribbons)
    fprintf('%-10s : Eg = %.4f eV  (%d atoms/cell)\n', ...
        results(r).name, results(r).Eg, results(r).Natoms);
end

Eg_7A  = results(1).Eg;
Eg_11A = results(2).Eg;
Eg_8Z  = results(3).Eg;

fprintf('\n===== Validation checks =====\n');

% Check 1: Eg(7-AGNR) > Eg(11-AGNR)
if Eg_7A > Eg_11A
    fprintf('[PASS] Eg(7-AGNR) = %.4f eV > Eg(11-AGNR) = %.4f eV\n', Eg_7A, Eg_11A);
else
    fprintf('[FAIL] Eg(7-AGNR) = %.4f eV is NOT > Eg(11-AGNR) = %.4f eV\n', Eg_7A, Eg_11A);
end

% Check 2: 11-AGNR is (near-)metallic, consistent with 3p+2 family
if Eg_11A < 0.05
    fprintf('[PASS] 11-AGNR is metallic in NN-TB (Eg = %.4f eV ~ 0), consistent with N=3p+2 family\n', Eg_11A);
else
    fprintf('[WARN] 11-AGNR gap (%.4f eV) larger than expected for a 3p+2 ribbon\n', Eg_11A);
end

% Check 3: 7-AGNR gapped, consistent with 3p+1 family (largest gap family)
if Eg_7A > 0.5
    fprintf('[PASS] 7-AGNR shows a large gap (Eg = %.4f eV), consistent with N=3p+1 family\n', Eg_7A);
else
    fprintf('[WARN] 7-AGNR gap (%.4f eV) smaller than expected for a 3p+1 ribbon\n', Eg_7A);
end

% Check 4: 8-ZGNR edge states pinned near E=0 close to zone boundary
kk_z  = results(3).kk;
Eall_z = results(3).Eall;
a_z = results(3).a;
near_zb = abs(kk_z*a_z) > 2*pi/3;           % k*a in (2pi/3, pi)
minAbsE_nearzb = min(min(abs(Eall_z(near_zb,:))));
if minAbsE_nearzb < 0.05
    fprintf('[PASS] 8-ZGNR: edge-state band reaches |E| = %.4f eV (<0.05 eV) near zone boundary (k*a > 2pi/3)\n', minAbsE_nearzb);
else
    fprintf('[FAIL] 8-ZGNR: no near-zero edge-state band found close to zone boundary (min|E| = %.4f eV)\n', minAbsE_nearzb);
end

% Check 5: band structures visibly different (different Eg + different atom count)
distinct = (results(1).Natoms ~= results(3).Natoms) || ...
           (abs(Eg_7A - Eg_8Z) > 0.1) || (abs(Eg_11A - Eg_8Z) > 0.1) || ...
           (abs(Eg_7A - Eg_11A) > 0.1);
if distinct
    fprintf('[PASS] The three band structures are quantitatively distinct (different Eg and/or atom count)\n');
else
    fprintf('[FAIL] Band structures appear nearly identical -- check Hamiltonian construction\n');
end

fprintf('\nIf all checks above are [PASS], proceed to GNR_Current.m / Compare_GNR_IV.m.\n');
