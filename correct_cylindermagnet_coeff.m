function [alpha,beta,dmin]=correct_cylindermagnet_coeff(MA,R,T,d)
% CORRECT_CYLINDERMAGNET_COEFF Computes magnetic force correction coefficients
% for cylindrical magnets by introducing α and β.
%
%   [ALPHA,BETA,DMIN] = CORRECT_CYLINDERMAGNET_COEFF(MA,R,T,d) calculates
%   correction coefficients α and β for magnetic force modeling between
%   identical cylindrical magnets. Both magnets share magnetization magnitude MA.
%
% Inputs:
%   MA  : Magnetization magnitude (A/m)
%   R   : Magnet radius (m)
%   T   : Half-thickness of magnet (m)
%   d   : Geometric spacing between magnets (m) [scalar or vector]
%
% Outputs:
%   alpha, beta : Dimensionless correction coefficients
%   dmin        : Lower spacing limit for valid application range (m)
%
% Reference: 
%   Yang, Y., Xiang, H. (2023). A simple and precise formula for magnetic
%   forces in nonlinear piezoelectric energy harvesting. Nonlinear Dynamics,
%   111, 6085-6110. https://doi.org/10.1007/s11071-023-08288-y
%
% Code developed by:
%   Yi Yang, Ph.D. Candidate
%   Beijing Jiaotong University
%   Contact: 19115045@bjtu.edu.cn
%
% Copyright (C) 2023-2025, Yi Yang
% Last updated: July 2025

% ========== OPTIMIZED SETUP ==========
num = 12001;  % Reduced number of points while maintaining accuracy
r = linspace(-60e-3, 60e-3, num);  % Displacement range
miu0 = 4*pi*1e-7;                  % Vacuum permeability
t = 2*T;                           % Full thickness of magnet
M = MA;                            % Magnetization
A = (1/4)*R^2 - (1/3)*T^2;         % Dimensionless area parameter

% Preallocate output arrays
alpha = zeros(size(d));
beta = zeros(size(d));
tic;
% ========== EFFICIENT PARALLEL COMPUTATION ==========
parfor k = 1:length(d)
    s = d(k) - 2*T;  % Effective spacing
    
    % Skip invalid spacing values
    if s <= 0
        alpha(k) = NaN;
        beta(k) = NaN;
        continue;
    end
    
    % Dynamically determine integration range - adaptive based on spacing
    max_x = min(500, max(100, 50/(s/R + 0.1)));  % Adaptive integration upper limit
    
    % Precompute constants
    pi_miu0_M2_R2 = pi * miu0 * M^2 * R^2;
    t_over_R = t/R;
    s_over_R = s/R;
    
    % Dynamic tolerance control
    if s < 1e-3  % Small spacing requires higher precision
        rel_tol = 1e-8;
        abs_tol = 1e-10;
    else
        rel_tol = 1e-6;
        abs_tol = 1e-8;
    end
    
    % Calculate force curve - only critical regions
    Fy = zeros(1, num);
    
    % Stage 1: Sparse sampling to locate approximate peak position
    sparse_indices = round(linspace(1, num, min(200, num)));
    for i = sparse_indices
        fun2 = @(x) besselj(1, r(i)*x/R)./x .* (besselj(1,x)).^2 .* ...
                    (1-exp(-x*t_over_R)).^2 .* exp(-x*s_over_R);
        Fy(i) = pi_miu0_M2_R2 * integral(fun2, 0, max_x, ...
                 'RelTol', rel_tol, 'AbsTol', abs_tol);
    end
    
    % Locate peak region
    [~, max_idx] = max(Fy(sparse_indices));
    max_idx = sparse_indices(max_idx);
    peak_region = max(1, max_idx-200):min(num, max_idx+200);
    
    % Stage 2: Compute only the peak region
    for i = peak_region
        if Fy(i) == 0  % Only compute uncomputed points
            fun2 = @(x) besselj(1, r(i)*x/R)./x .* (besselj(1,x)).^2 .* ...
                        (1-exp(-x*t_over_R)).^2 .* exp(-x*s_over_R);
            Fy(i) = pi_miu0_M2_R2 * integral(fun2, 0, max_x, ...
                     'RelTol', rel_tol, 'AbsTol', abs_tol);
        end
    end
    
    % ========== ACCURATE PEAK DETECTION ==========
    % Find the maximum point in the peak region
    [~, local_max_idx] = max(Fy(peak_region));
    max_idx = peak_region(1) + local_max_idx - 1;
    
    % Ensure enough points for interpolation
    idx_low = max(1, max_idx-2);
    idx_high = min(num, max_idx+2);
    r_segment = r(idx_low:idx_high);
    Fy_segment = Fy(idx_low:idx_high);
    
    % Use quadratic interpolation for peak refinement
    if length(r_segment) >= 3
        p = polyfit(r_segment, Fy_segment, 2);
        w_tip_max = -p(2)/(2*p(1));  % Vertex of parabola
        
        % Clamp to segment range
        w_tip_max = max(min(w_tip_max, max(r_segment)), min(r_segment));
    else
        w_tip_max = r(max_idx);
    end
    
    % Dimensionless parameters
    miu_val = A / d(k)^2;
    eta = abs(w_tip_max / d(k));
    
    % Compute correction coefficients
    beta(k) = 3.309 + 10.09*miu_val - 9.176*eta + 21.74*miu_val^2 - ...
              23.47*miu_val*eta + 8.336*eta^2;
    alpha(k) = 0.04007 - 39.13*miu_val + 2.1*eta - 48.77*miu_val^2 + ...
               68.64*miu_val*eta + 5.643*eta^2 - 12.06*eta^3;
end

% Calculate minimum spacing
R1 = sqrt(max(0, R^2 - (4/3)*T^2));
dmin = 1.408*R1 - 0.6363*R + 0.7451*t;
toc;
end