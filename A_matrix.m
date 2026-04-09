function A = A_matrix(beta, L, lp, EI1, EI2, m01, It, Mt, beta_coeff2)
% Builds the 8x8 characteristic matrix for the segmented beam
beta2 = beta * beta_coeff2;
w1_sq = beta^4 * (EI1/m01);

A = zeros(8,8);
% x=0: displacement zero
A(1,1:4) = [1, 0, 1, 0];
% x=0: slope zero
A(2,1:4) = [0, beta, 0, beta];
% x=L: moment balance
A(3,5:8) = EI2*beta2^2*[-cos(beta2*L), -sin(beta2*L), cosh(beta2*L), sinh(beta2*L)] - ...
           w1_sq*It*beta2*[-sin(beta2*L), cos(beta2*L), sinh(beta2*L), cosh(beta2*L)];
% x=L: shear balance
A(4,5:8) = EI2*beta2^3*[sin(beta2*L), -cos(beta2*L), sinh(beta2*L), cosh(beta2*L)] + ...
           w1_sq*Mt*[cos(beta2*L), sin(beta2*L), cosh(beta2*L), sinh(beta2*L)];
% x=lp: displacement continuity
A(5,1:4) = [cos(beta*lp), sin(beta*lp), cosh(beta*lp), sinh(beta*lp)];
A(5,5:8) = -[cos(beta2*lp), sin(beta2*lp), cosh(beta2*lp), sinh(beta2*lp)];
% x=lp: slope continuity
A(6,1:4) = beta*[-sin(beta*lp), cos(beta*lp), sinh(beta*lp), cosh(beta*lp)];
A(6,5:8) = -beta2*[-sin(beta2*lp), cos(beta2*lp), sinh(beta2*lp), cosh(beta2*lp)];
% x=lp: moment continuity
A(7,1:4) = EI1*beta^2*[-cos(beta*lp), -sin(beta*lp), cosh(beta*lp), sinh(beta*lp)];
A(7,5:8) = -EI2*beta2^2*[-cos(beta2*lp), -sin(beta2*lp), cosh(beta2*lp), sinh(beta2*lp)];
% x=lp: shear continuity
A(8,1:4) = EI1*beta^3*[sin(beta*lp), -cos(beta*lp), sinh(beta*lp), cosh(beta*lp)];
A(8,5:8) = -EI2*beta2^3*[sin(beta2*lp), -cos(beta2*lp), sinh(beta2*lp), cosh(beta2*lp)];
end