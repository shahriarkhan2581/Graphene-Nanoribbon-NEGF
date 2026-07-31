function [H00,H01,coords,a] = H_7AGNR()

% 7-AGNR Wrapper

ac = 1.42;
N = 7;

[H00,H01,coords] = AGNR_Hamiltonian(N);

a = 3*ac;

end
