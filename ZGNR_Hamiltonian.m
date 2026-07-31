function [H00, H01, coords] = ZGNR_Hamiltonian(N)
%ZGNR_HAMILTONIAN  Coordinate-based nearest-neighbour pi-orbital
%tight-binding Hamiltonian for an N-ZGNR (zigzag-edge graphene
%nanoribbon), transport direction = ribbon axis (zigzag direction).
%
%   [H00,H01,coords] = ZGNR_Hamiltonian(N)
%
%   N      : number of zigzag chains across ribbon width (ZGNR index)
%   H00    : (2N x 2N) intra-unit-cell Hamiltonian block
%   H01    : (2N x 2N) inter-unit-cell hopping block, H(cell n -> cell n+1)
%   coords : (2N x 2) real-space atomic coordinates (Angstrom) in one cell
%
%   Convention: H(k) = H00 + H01*exp(1i*k*a) + H01'*exp(-1i*k*a)
%   with a = sqrt(3)*a_cc the ZGNR lattice constant (transport period).
%
%   Method: identical geometric/distance-based construction as
%   AGNR_Hamiltonian.m, but the roles of transport/width directions are
%   swapped (ZGNR transport period is sqrt(3)*a_cc, the "short" period;
%   width grows along y). This automatically produces the characteristic
%   flat zigzag edge-state bands near E = 0.

t0 = -2.7;      % nearest-neighbour hopping integral, eV
ac = 1.42;      % C-C bond length, Angstrom
tol = 1e-3;     % geometric window tolerance, Angstrom
distTol = 0.01*ac;

a1 = ac*[ sqrt(3)/2,  3/2];
a2 = ac*[-sqrt(3)/2,  3/2];
d1 = ac*[0, 1];

Ttrans = [sqrt(3)*ac, 0];        % ZGNR transport (periodic) vector
Xmax   = sqrt(3)*ac;             % one ZGNR unit cell along transport (x)
Ymin   = -0.5*ac;                 % width-window offset: starts the cut on
% a proper 2-fold-coordinated edge atom
% (a naive cut at y=0 lands on an atom
% with only 1 remaining bond, which is
% not a physical zigzag edge)
Ymax   = Ymin + 1.5*N*ac;         % ribbon width window (y-direction)

R = N + 6;
A_all = [];
B_all = [];
for n1 = -R:R
    for n2 = -R:R
        A = n1*a1 + n2*a2;
        B = A + d1;
        A_all = [A_all; A]; %#ok<AGROW>
        B_all = [B_all; B]; %#ok<AGROW>
    end
end
allpts = [A_all; B_all];

mask = allpts(:,1) >= -tol & allpts(:,1) < Xmax - tol & ...
    allpts(:,2) >= Ymin - tol & allpts(:,2) < Ymax - tol;
coords = unique(round(allpts(mask,:)*1e6)/1e6, 'rows');
coords = sortrows(coords, [2 1]);

Natoms = size(coords,1);
if Natoms ~= 2*N
    warning('ZGNR_Hamiltonian:atomCount', ...
        'Expected %d atoms for %d-ZGNR, generated %d.', 2*N, N, Natoms);
end

H00 = zeros(Natoms);
H01 = zeros(Natoms);

for i = 1:Natoms
    for j = 1:Natoms
        if i == j, continue; end
        d = norm(coords(i,:) - coords(j,:));
        if abs(d - ac) < distTol
            H00(i,j) = t0;
        end
    end
end

coords_next = coords + Ttrans;
for i = 1:Natoms
    for j = 1:Natoms
        d = norm(coords(i,:) - coords_next(j,:));
        if abs(d - ac) < distTol
            H01(i,j) = t0;
        end
    end
end

end
