function [H00,H01,coords,a] = H_11AGNR()

% 11-AGNR Wrapper

ac = 1.42;
N = 11;

[H00,H01,coords] = AGNR_Hamiltonian(N);

a = 3*ac;

end
