function A = ellipsoidal_basis_gen(Jmax, t)
S = zeros(length(t), Jmax);
C = zeros(length(t), Jmax);
for j = 1:Jmax
    S(:,j) = sin(j*t);
    C(:,j) = cos(j*t);
end

% construct a basis matrix
A = [ones(length(t),1) C S];