function [H00,H01,coords,a] = H_8ZGNR()

% 8-ZGNR Wrapper

ac = 1.42;
N = 8;

[H00,H01,coords] = ZGNR_Hamiltonian(N);

a = sqrt(3)*ac;

end
