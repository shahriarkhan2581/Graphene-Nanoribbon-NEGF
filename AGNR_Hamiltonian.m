function [H00, H01, coords] = AGNR_Hamiltonian(N)
%AGNR_HAMILTONIAN  Coordinate-based nearest-neighbour pi-orbital
%tight-binding Hamiltonian for an N-AGNR (armchair-edge graphene
%nanoribbon), transport direction = ribbon axis (armchair direction).
%
%   [H00,H01,coords] = AGNR_Hamiltonian(N)
%
%   N      : number of dimer lines across ribbon width (AGNR index)
%   H00    : (2N x 2N) intra-unit-cell Hamiltonian block
%   H01    : (2N x 2N) inter-unit-cell hopping block, H(cell n -> cell n+1)
%   coords : (2N x 2) real-space atomic coordinates (Angstrom) in one cell
%
%   Convention: H(k) = H00 + H01*exp(1i*k*a) + H01'*exp(-1i*k*a)
%   with a = 3*a_cc the AGNR lattice constant (transport period).
%
%   Method: atoms are generated directly from the honeycomb lattice
%   (A/B sublattice basis vectors), windowed to one AGNR unit cell.
%   Bonds (H00, H01) are then assigned purely by nearest-neighbour
%   distance (~a_cc), so the topology is guaranteed correct regardless
%   of manual index bookkeeping.

t0 = -2.7;      % nearest-neighbour hopping integral, eV
ac = 1.42;      % C-C bond length, Angstrom
tol = 1e-3;     % geometric window tolerance, Angstrom
distTol = 0.01*ac; % bond-length tolerance for nearest-neighbour test

% Graphene honeycomb lattice vectors (A sublattice), B = A + d1
a1 = ac*[ sqrt(3)/2,  3/2];
a2 = ac*[-sqrt(3)/2,  3/2];
d1 = ac*[0, 1];

Ttrans = [0, 3*ac];              % AGNR transport (periodic) vector
Wmax   = N*(sqrt(3)/2)*ac;       % ribbon width window (x-direction)
Ymax   = 3*ac;                   % one AGNR unit cell along transport (y)

% Enumerate candidate lattice points over a generous range
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

% Keep only points within the AGNR unit-cell window [0,Wmax) x [0,Ymax)
mask = allpts(:,1) >= -tol & allpts(:,1) < Wmax - tol & ...
       allpts(:,2) >= -tol & allpts(:,2) < Ymax - tol;
coords = unique(round(allpts(mask,:)*1e6)/1e6, 'rows');
coords = sortrows(coords, [2 1]);   % order by y then x

Natoms = size(coords,1);
if Natoms ~= 2*N
    warning('AGNR_Hamiltonian:atomCount', ...
        'Expected %d atoms for %d-AGNR, generated %d.', 2*N, N, Natoms);
end

H00 = zeros(Natoms);
H01 = zeros(Natoms);

% Intra-cell bonds
for i = 1:Natoms
    for j = 1:Natoms
        if i == j, continue; end
        d = norm(coords(i,:) - coords(j,:));
        if abs(d - ac) < distTol
            H00(i,j) = t0;
        end
    end
end

% Inter-cell bonds: cell 0 (coords) -> cell +1 (coords shifted by Ttrans)
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
