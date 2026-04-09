function [alpha,beta,dmin]=correct_cuboidalmagnet_coeff(MA,W,H,L,d)
% CORRECT_CUBOIDALMAGNET_COEFF Computes magnetic force correction coefficients
% for cuboidal magnets with square cross-section by introducing α and β.
%-------------------------------------------------------------------------
%   [alpha,beta,dmin] = CORRECT_CUBOIDALMAGNET_COEFF(MA,W,H,L,d) calculates
%   correction coefficients α and β for magnetic force modeling between
%   identical cuboidal magnets. Both magnets share magnetization magnitude MA.
%-------------------------------------------------------------------------
% Inputs:
%   MA  : Magnetization magnitude (A/m)
%   W   : Half-width of magnet cross-section (m)
%   H   : Half-height of magnet (m)
%   L   : Half-length of magnet (m)
%   d   : Geometric spacing between magnets (m) [scalar or vector]
%-------------------------------------------------------------------------
% Outputs:
%   alpha, beta : Dimensionless correction coefficients
%   dmin        : Lower spacing limit for valid application range (m)
%-------------------------------------------------------------------------
% Notes:
%   - Magnet dimensions are defined using half-parameters following standard
%     cuboid magnet modeling conventions
%   - Input d can be a vector for batch processing
%-------------------------------------------------------------------------
% Reference: 
%   Yang, Y., Xiang, H. (2023). A simple and precise formula for magnetic
%   forces in nonlinear piezoelectric energy harvesting. Nonlinear Dynamics,
%   111, 6085-6110. https://doi.org/10.1007/s11071-023-08288-y
%-------------------------------------------------------------------------
% Code developed by:
%   Yi Yang, Ph.D. Candidate
%   Beijing Jiaotong University
%   Contact: 19115045@bjtu.edu.cn
%-------------------------------------------------------------------------
% Copyright (C) 2023-2025, Yi Yang
% Last updated: July 2025
%-------------------------------------------------------------------------

% Number of discretization points for force calculation
num = 12001;  
w0 = linspace(-60e-3, 60e-3, num);

% Calculate minimum valid spacing for the analytical model
dmin = 1.587*sqrt(H^2 - L^2) - 0.5577*2*H + 1.589*2*L;

% Preallocate output arrays
alpha = zeros(size(d));
beta = zeros(size(d));

% Precompute trigonometric constants (theta = π)
c60 = cos(pi); s60 = sin(pi);  

% Parallel processing over spacing values
parfor k = 1:length(d)
    % Magnet dimensions (fixed per iteration)
    a0 = W; b0 = H; c0 = L; 
    A0 = a0; B0 = b0; C0 = c0;
    
    % Effective spacing and magnetization
    d0 = d(k) - 2*L;
    M0 = MA;
    
    % Precompute doubled dimensions (used repeatedly)
    a_val = 2*a0; b_val = 2*b0; c_val = 2*c0;
    A_val = 2*A0; B_val = 2*B0; C_val = 2*C0;
    z01_val = d0 + 2*C0 + 2*c0;  % Fixed z-position
    
    % Initialize force array
    Fmy = zeros(1, num);
    
    % Calculate magnetic force at each y-displacement
    for m = 1:num
        y01 = 2*B0 + w0(m);  % Varying y-position
        
        % Compute force components using analytical model
        term1 = calc_term(0, A_val, y01, z01_val, c60, s60, 0, 0, B_val, C_val, M0);
        term2 = calc_term(-a_val, -a_val + A_val, y01, z01_val, c60, s60, 0, 0, B_val, C_val, M0);
        term3 = calc_term(-a_val, -a_val + A_val, y01, z01_val, c60, s60, b_val, 0, B_val, C_val, M0);
        term4 = calc_term(0, A_val, y01, z01_val, c60, s60, b_val, 0, B_val, C_val, M0);
        term5 = calc_term(-a_val, -a_val + A_val, y01, z01_val, c60, s60, 0, c_val, B_val, C_val, M0);
        term6 = calc_term(0, A_val, y01, z01_val, c60, s60, 0, c_val, B_val, C_val, M0);
        term7 = calc_term(0, A_val, y01, z01_val, c60, s60, b_val, c_val, B_val, C_val, M0);
        term8 = calc_term(-a_val, -a_val + A_val, y01, z01_val, c60, s60, b_val, c_val, B_val, C_val, M0);
        
        % Combine components to get total force
        Fmy(m) = term1 - term2 + term3 - term4 + term5 - term6 + term7 - term8;
    end
    
    % Numerical integration using trapezoidal rule
    Um0 = zeros(1, num);
    for m = 2:num
        Um0(m) = 0.5*(Fmy(m) + Fmy(m-1))*(w0(m) - w0(m-1)) + Um0(m-1);
    end

    % Find force maximum and its position
    [~, max_idx] = max(Fmy);
    w0_max = w0(max_idx);

    % Calculate dimensionless parameters
    A_sq = A0^2 - C0^2;
    miu = A_sq / d(k)^2;          % Dimensionless area parameter
    eta = abs(w0_max / d(k));     % Normalized displacement
    
    % Compute correction coefficients
    beta(k) = 3.166 - 4.637*miu - 12.32*eta - 10.46*miu^2 + ...
              4.232*miu*eta + 26.07*eta^2 - 22.36*eta^3 + 16.04*miu*eta^2;
    alpha(k) = 0.3065 - 22.94*miu - 4.872*eta - 25.12*miu^2 + ...
               33.56*miu*eta + 34.36*eta^2 - 44.23*eta^3 + 27.88*miu*eta^2;
end

end

% Helper function: Calculates force component term
function val = calc_term(v, w, y01, z01, ctheta, stheta, b, c, B, C, M)
% Compute force component using analytical expression
val = 1e-7 * M^2 * (f2(v, w, y01, z01, ctheta, stheta, b, c, B, C) - ...
                  f2(v, w, y01, z01, ctheta, stheta, b, c, B, 0));
end

% Helper function: Intermediate force calculation
function res = f2(v, w, y01, z01, ctheta, stheta, b, c, B, z0)
% Calculate force difference components
res = f3(w, y01, z01, ctheta, stheta, b, c, B, z0) - ...
      f3(v, y01, z01, ctheta, stheta, b, c, B, z0) - ...
      f3(w, y01, z01, ctheta, stheta, b, c, 0, z0) + ...
      f3(v, y01, z01, ctheta, stheta, b, c, 0, z0);
end

% Core magnetic field calculation function
function res = f3(u, y01, z01, ctheta, stheta, b, c, y0, z0)
% Compute fundamental magnetic field components

% 1. Calculate spatial transformation terms
f6 = y01*ctheta + z01*stheta - b*ctheta - c*stheta + y0;
f5 = -y01*stheta + z01*ctheta + b*stheta - c*ctheta + z0;

% 2. Calculate distance term with numerical safety
f4_sq = u^2 + f5^2 + f6^2;
if f4_sq < 0
    f4 = 0;  % Prevent complex numbers (physical safeguard)
else
    f4 = sqrt(f4_sq);
end

% 3. Logarithmic term calculations with domain checks
log_arg1 = -u + f4;
if log_arg1 <= 0
    part1 = 0;  % Handle invalid log domain
else
    part1 = u * f6 * log(log_arg1);
end

log_arg2 = f4 + f6;
if log_arg2 <= 0
    part2 = 0;
    part6 = 0;
else
    part2 = -u^2 * log(log_arg2);
    part6 = 0.5*(u^2 + f5^2) * log(log_arg2);
end

% 4. Arctangent term with singularity protection
denom = f5 * f6;
if abs(denom) < 1e-10
    part3 = 0;  % Avoid division by zero
else
    atan_arg = (-f5^2 - u^2 + u*f4) / denom;
    part3 = u * f5 * atan(atan_arg);
end

% 5. Special function terms
part4 = 0.5 * u * pi * abs(f5) * sign(f6);  % Sign-dependent term
part5 = 0.5 * f6 * f4;  % Direct spatial term

% 6. Combine all components
res = part1 - u*f6 + part2 + part3 + part4 + part5 + part6;
end