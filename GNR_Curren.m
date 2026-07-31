function [V, I, T_of_E, Egrid] = GNR_Current(H00, H01, varargin)
%GNR_CURRENT  Ballistic I-V characteristic of a graphene nanoribbon via
%NEGF + Landauer formalism, with semi-infinite leads made of the SAME
%ribbon material (identical H00,H01), so transmission T(E) reduces to
%the number of conducting sub-bands ("modes") at energy E, broadened by
%the finite imaginary part used in the retarded Green's function.
%
%   [V,I] = GNR_Current(H00,H01)
%   [V,I] = GNR_Current(H00,H01,'Name',value,...)
%
%   H00, H01 : unit-cell Hamiltonian blocks from AGNR_Hamiltonian.m /
%              ZGNR_Hamiltonian.m (H01 = hopping cell(n) -> cell(n+1))
%
%   Optional name-value pairs:
%     'Ef'      Fermi level, eV                  (default 0)
%     'Vmax'    maximum bias magnitude, V         (default 1.0)
%     'Nv'      number of bias points             (default 41)
%     'NE'      number of energy grid points      (default 601)
%     'Erange'  [Emin Emax] energy window, eV     (default [-3 3])
%     'eta'     small imaginary part, eV          (default 1e-4)
%     'T'       temperature, K                    (default 300)
%
%   Outputs:
%     V       : (Nv x 1) bias vector, V
%     I       : (Nv x 1) current (arbitrary units of 2e/h, i.e. I is in
%               units of [2e/h * eV] so this is a relative/comparative
%               current -- valid for comparing ribbons, see note below)
%     T_of_E  : (NE x Nv) transmission at each energy/bias (T at mu
%               window used); also see Egrid
%     Egrid   : (NE x 1) energy grid used for the E=0-bias transmission
%
%   Method:
%     1) Lopez-Sancho-Sancho-Rubio iterative decimation gives the
%        surface Green's function of each semi-infinite lead.
%     2) Self-energies: Sigma_L = H01' * gL * H01 ,
%                        Sigma_R = H01  * gR * H01'
%     3) Broadening: Gamma_L,R = i*(Sigma_L,R - Sigma_L,R')
%     4) Device Green's function (device = one unit cell, H00):
%          G(E) = [(E+i*eta)*I - H00 - Sigma_L(E) - Sigma_R(E)]^-1
%     5) Transmission: T(E) = Trace[Gamma_L * G * Gamma_R * G']
%     6) Landauer current:
%          I(V) = (2e/h) * Integral T(E) [f(E-muL)-f(E-muR)] dE
%        with muL = Ef+eV/2, muR = Ef-eV/2. Implemented with e=h=1
%        (i.e. I returned in units of 2e/h * eV); use physical prefactor
%        (2e/h = 7.748e-5 S) to convert to Amperes if needed.

p = inputParser;
addParameter(p,'Ef',0);
addParameter(p,'Vmax',1.0);
addParameter(p,'Nv',41);
addParameter(p,'NE',601);
addParameter(p,'Erange',[-3 3]);
addParameter(p,'eta',1e-4);
addParameter(p,'T',300);
parse(p,varargin{:});
opt = p.Results;

kB = 8.617333e-5;    % eV/K
Ef = opt.Ef;
Natoms = size(H00,1);
Imat = eye(Natoms);

V = linspace(-opt.Vmax, opt.Vmax, opt.Nv)';
Egrid = linspace(opt.Erange(1), opt.Erange(2), opt.NE)';

% ---- Precompute T(E) once on the fixed energy grid (bias-independent,
%      since the channel material itself is unbiased/rigid: standard
%      simplification for a ballistic homogeneous-lead 2-terminal wire).
Tvec = zeros(opt.NE,1);
for iE = 1:opt.NE
    E = Egrid(iE);
    % lopez_sancho_surface_gf(E,H00,X,eta) grows a semi-infinite chain in
    % the +X direction from the surface layer. The LEFT lead extends in
    % the -H01 direction (pass H01' to grow it); the RIGHT lead extends
    % in the +H01 direction (pass H01).
    gL = lopez_sancho_surface_gf(E, H00, H01', opt.eta);
    gR = lopez_sancho_surface_gf(E, H00, H01,  opt.eta);

    SigL = H01' * gL * H01;
    SigR = H01  * gR * H01';

    GamL = 1i*(SigL - SigL');
    GamR = 1i*(SigR - SigR');

    G = ((E + 1i*opt.eta)*Imat - H00 - SigL - SigR) \ Imat;

    Tvec(iE) = real(trace(GamL*G*GamR*G'));
end
Tvec(Tvec < 0) = 0;   % numerical noise guard

T_of_E = Tvec;

% ---- Landauer current integration for each bias point
I = zeros(opt.Nv,1);
for iv = 1:opt.Nv
    muL = Ef + V(iv)/2;
    muR = Ef - V(iv)/2;
    fL = 1./(1+exp((Egrid-muL)/(kB*opt.T)));
    fR = 1./(1+exp((Egrid-muR)/(kB*opt.T)));
    integrand = Tvec.*(fL-fR);
    I(iv) = trapz(Egrid, integrand);   % units: 2e/h * eV
end

end

% ================= local function =================
function gs = lopez_sancho_surface_gf(E, H00, H01, eta)
%LOPEZ_SANCHO_SURFACE_GF  Iterative (Sancho-Sancho-Rubio 1985) surface
%Green's function of a semi-infinite periodic lead described by onsite
%block H00 and hopping block H01 (H01 = hopping from a layer to the
%NEXT layer, i.e. the lead grows in the +H01 direction).
%
%   gs = lopez_sancho_surface_gf(E,H00,H01,eta)

maxIter = 200;
tolIter = 1e-10;

n = size(H00,1);
Imat = eye(n);
w = (E + 1i*eta)*Imat;

eps_i  = H00;
eps_s  = H00;
alpha  = H01;
beta   = H01';

for iter = 1:maxIter
    g = (w - eps_i) \ Imat;
    ag  = alpha*g;
    bg  = beta*g;

    eps_s_new = eps_s + ag*beta;
    eps_i_new = eps_i + ag*beta + bg*alpha;
    alpha_new = ag*alpha;
    beta_new  = bg*beta;

    if norm(alpha_new,'fro') < tolIter && norm(beta_new,'fro') < tolIter
        eps_s = eps_s_new;
        break;
    end

    eps_s = eps_s_new;
    eps_i = eps_i_new;
    alpha = alpha_new;
    beta  = beta_new;
end

gs = (w - eps_s) \ Imat;

end
