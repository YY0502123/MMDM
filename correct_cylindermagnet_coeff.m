function [alpha,beta,dmin]=correct_cylindermagnet_coeff(MA,R,T,d)
% CORRECT_CYLINDERMAGNET_COEFF Computes magnetic force correction coefficients
% for cylindrical magnets by introducing dimensionless α and β.
%
%   [ALPHA,BETA,DMIN] = CORRECT_CYLINDERMAGNET_COEFF(MA,R,T,d) calculates
%   correction coefficients α and β for magnetic force modeling between
%   two identical face-to-face cylindrical magnets sharing magnetization MA.
%
%--------------------------------------------------------------------------
% Workflow:
%   1. Polynomial region  (dmin < d ≤ 150 mm):
%        α, β from least-squares polynomial fit over calibration dataset.
%        Exact Hankel-transform magnetostatics integral used internally to
%        locate the transverse force peak (→ dimensionless parameters μ, η).
%
%   2. Exponential extrapolation (d > 150 mm):
%        α(d) = 1 + (α₁₅₀−1)·exp(−λ_α·(d−150 mm))
%        λ fitted by origin-constrained log-space regression using reference
%        distances [30,60,90,120,150] mm (points below dmin auto-skipped).
%        No hard upper cutoff; naturally approaches 1 as d → ∞.
%
%   3. d ≤ dmin: outputs NaN (model geometrically inapplicable).
%--------------------------------------------------------------------------
% Applicable geometry ranges (calibration dataset):
%   R   : magnet radius,          2~10 mm
%   2T  : total magnet thickness, 2~10 mm
%   d   : inter-magnet spacing,   dmin ~ ∞ (extrapolation beyond 150 mm)
%--------------------------------------------------------------------------
% Exact force model used internally (Hankel-transform integral):
%   Fy(r) = π μ₀ M² R² ∫₀^∞ J₁(ξr/R)/ξ · J₁²(ξ) ·
%                        (1−e^{−ξt/R})² · e^{−ξs/R} dξ
%   where s = d − 2T is the surface-to-surface gap.
%   Computed numerically using MATLAB adaptive quadrature (integral).
%--------------------------------------------------------------------------
% Inputs:
%   MA : Magnetization magnitude (A/m)
%   R  : Magnet radius (m)
%   T  : Magnet half-thickness (m)
%   d  : Center-to-center spacing (m) [scalar or vector]
%
% Outputs:
%   alpha : Dimensionless correction coefficient (→ 1 as d → ∞)
%   beta  : Dimensionless correction coefficient (→ 1 as d → ∞)
%   dmin  : Lower bound of valid spacing range (m)
%--------------------------------------------------------------------------
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
% Copyright (C) 2023-2026, Yi Yang
% Last updated: May 2026
%--------------------------------------------------------------------------

% ===== CONSTANTS =====
miu0    = 4*pi*1e-7;                % Vacuum permeability (H/m)
t       = 2*T;                      % Full magnet thickness (m)
M       = MA;                       % Magnetization (A/m)
A       = (1/4)*R^2 - (1/3)*T^2;   % Geometric parameter for dimensionless miu
d_trans = 150e-3;                   % Polynomial upper limit (m); extrapolation beyond
tol_fit = 1e-5;                     % Min deviation threshold for λ regression

% ===== MINIMUM VALID SPACING =====
% Derived from magnet geometry; below this the model is inapplicable.
R1   = sqrt(max(0, R^2 - (4/3)*T^2));
dmin = 1.408*R1 - 0.6363*R + 0.7451*t;

% ===== PRE-COMPUTE REFERENCE VALUES AND FIT λ =====
% Evaluate polynomial at multiple reference distances to robustly fit the
% exponential decay rate λ for extrapolation beyond d_trans.
% Reference distances: [30,60,90,120,150] mm; points ≤ dmin are auto-skipped.
d_refs_all = [30, 60, 90, 120, 150] * 1e-3;
valid_ref  = (d_refs_all > dmin) & (d_refs_all <= d_trans);

if sum(valid_ref) < 1
    % Fallback: no valid reference points (very large dmin); set λ = 0
    warning('correct_cylindermagnet_coeff: no valid reference points for lambda fitting.');
    lam_a = 0;
    lam_b = 0;
    [alpha_trans, beta_trans] = eval_cyl_poly(d_trans, R, T, M, miu0, A);

else
    d_refs     = d_refs_all(valid_ref);
    alpha_refs = zeros(1, length(d_refs));
    beta_refs  = zeros(1, length(d_refs));

    % Sequential evaluation at reference distances (pre-parfor)
    for ki = 1:length(d_refs)
        [alpha_refs(ki), beta_refs(ki)] = eval_cyl_poly(d_refs(ki), R, T, M, miu0, A);
    end

    % Anchor values at d_trans (guarantees C⁰ continuity of extrapolation)
    if abs(d_refs(end) - d_trans) < 1e-10
        alpha_trans = alpha_refs(end);
        beta_trans  = beta_refs(end);
    else
        [alpha_trans, beta_trans] = eval_cyl_poly(d_trans, R, T, M, miu0, A);
        d_refs     = [d_refs,     d_trans     ];
        alpha_refs = [alpha_refs, alpha_trans ];
        beta_refs  = [beta_refs,  beta_trans  ];
    end

    %----------------------------------------------------------------------
    % Fit λ via origin-constrained log-space linear regression (inlined)
    %
    % Model: f(d) = 1 + (f_trans−1)·exp(−λ·(d−d_trans))
    % At each reference point d_i < d_trans:
    %   log(|f_i−1| / |f_trans−1|) = −λ·(d_i − d_trans)
    % Let x_i = d_i−d_trans (<0),  y_i = log(|dev_i/dev_trans|) (>0):
    %   λ = −Σ(x_i·y_i) / Σ(x_i²)   [origin-constrained OLS]
    % Points at d_i = d_trans contribute 0/0 and are excluded automatically.
    %----------------------------------------------------------------------

    % ---- λ for alpha ----
    dev_trans_a = alpha_trans - 1;
    if abs(dev_trans_a) < tol_fit
        lam_a = 0;   % Deviation already negligible; no decay needed
    else
        dev_a  = alpha_refs - 1;
        keep_a = (abs(dev_a) > tol_fit) & (abs(d_refs - d_trans) > 1e-10);
        if sum(keep_a) < 1
            lam_a = 0;
        else
            xa    = d_refs(keep_a) - d_trans;
            ya    = log(abs(dev_a(keep_a) ./ dev_trans_a));
            lam_a = max(0, -sum(xa .* ya) / sum(xa.^2));
        end
    end

    % ---- λ for beta ----
    dev_trans_b = beta_trans - 1;
    if abs(dev_trans_b) < tol_fit
        lam_b = 0;
    else
        dev_b  = beta_refs - 1;
        keep_b = (abs(dev_b) > tol_fit) & (abs(d_refs - d_trans) > 1e-10);
        if sum(keep_b) < 1
            lam_b = 0;
        else
            xb    = d_refs(keep_b) - d_trans;
            yb    = log(abs(dev_b(keep_b) ./ dev_trans_b));
            lam_b = max(0, -sum(xb .* yb) / sum(xb.^2));
        end
    end

end % if sum(valid_ref)

% ===== MAIN COMPUTATION =====
alpha = zeros(size(d));
beta  = zeros(size(d));

tic;
parfor k = 1:length(d)

    if d(k) <= dmin
        alpha(k) = NaN;
        beta(k)  = NaN;
    elseif d(k) > d_trans
        %------------------------------------------------------------------
        % Exponential extrapolation region (d > 150 mm)
        % α(d) = 1 + (α_trans−1)·exp(−λ_α·(d−d_trans))
        % Naturally approaches 1 as d → ∞ with no hard upper cutoff.
        %------------------------------------------------------------------
        if lam_a > 0
            alpha(k) = 1 + (alpha_trans - 1) * exp(-lam_a * (d(k) - d_trans));
        else
            alpha(k) = 1;
        end
        if lam_b > 0
            beta(k)  = 1 + (beta_trans  - 1) * exp(-lam_b * (d(k) - d_trans));
        else
            beta(k) = 1;
        end

    elseif d(k) - 2*T <= 0
        %------------------------------------------------------------------
        % d ≤ dmin: magnets overlap or touch; model inapplicable
        %------------------------------------------------------------------
        alpha(k) = NaN;
        beta(k)  = NaN;

    else
        %------------------------------------------------------------------
        % Polynomial region (dmin < d ≤ 150 mm)
        % Exact force integral → peak location → miu, η → polynomial
        %------------------------------------------------------------------
        [alpha(k), beta(k)] = eval_cyl_poly(d(k), R, T, M, miu0, A);

    end

end % parfor
toc;

end % main function


%==========================================================================
% LOCAL FUNCTION: eval_cyl_poly
%   Evaluates polynomial correction at a single spacing dk.
%   Computes exact Hankel-transform Fy(r), locates peak, derives miu & η.
%==========================================================================
function [alpha_raw, beta_raw] = eval_cyl_poly(dk, R, T, M, miu0, A)

num  = 12001;
r    = linspace(-60e-3, 60e-3, num);
t    = 2*T;
s    = dk - 2*T;   % Surface-to-surface gap (m)

if s <= 0
    alpha_raw = NaN;  beta_raw = NaN;  return;
end

% Adaptive integration upper limit: larger gap → faster decay → smaller limit
max_x         = min(500, max(100, 50/(s/R + 0.1)));
pi_miu0_M2_R2 = pi * miu0 * M^2 * R^2;   % Integral prefactor (N)
t_over_R      = t / R;
s_over_R      = s / R;

% Tighter tolerance for near-contact (steep integrand)
if s < 1e-3
    rel_tol = 1e-8;  abs_tol = 1e-10;
else
    rel_tol = 1e-6;  abs_tol = 1e-8;
end

% ---- Exact transverse force Fy(r) via Hankel-transform integral ----
%   Kernel: J₁(ξr/R)/ξ · J₁²(ξ) · (1−e^{−ξt/R})² · e^{−ξs/R}
Fy = zeros(1, num);

% Stage 1: Coarse probe (~200 pts) to approximately locate the force peak
sp = round(linspace(1, num, min(200, num)));
for i = sp
    fun2  = @(x) besselj(1, r(i)*x/R)./x .* (besselj(1,x)).^2 .* ...
                 (1-exp(-x*t_over_R)).^2 .* exp(-x*s_over_R);
    Fy(i) = pi_miu0_M2_R2 * integral(fun2, 0, max_x, ...
             'RelTol',rel_tol, 'AbsTol',abs_tol);
end

% Stage 2: Dense evaluation in ±200-point neighbourhood of identified peak
[~, mx]   = max(Fy(sp));
mx        = sp(mx);
pk_region = max(1, mx-200) : min(num, mx+200);
for i = pk_region
    if Fy(i) == 0   % Skip points already evaluated in Stage 1
        fun2  = @(x) besselj(1, r(i)*x/R)./x .* (besselj(1,x)).^2 .* ...
                     (1-exp(-x*t_over_R)).^2 .* exp(-x*s_over_R);
        Fy(i) = pi_miu0_M2_R2 * integral(fun2, 0, max_x, ...
                 'RelTol',rel_tol, 'AbsTol',abs_tol);
    end
end

% Sub-grid peak refinement via local quadratic (parabolic) fit
[~, lmx]  = max(Fy(pk_region));
mx        = pk_region(1) + lmx - 1;
il        = max(1,   mx-2);
ih        = min(num, mx+2);
rs        = r(il:ih);
Fs        = Fy(il:ih);
if length(rs) >= 3
    pp        = polyfit(rs, Fs, 2);       % Fit ax²+bx+c
    w_tip_max = -pp(2) / (2*pp(1));       % Vertex: x* = −b/(2a)
    w_tip_max = max(min(w_tip_max, max(rs)), min(rs));
else
    w_tip_max = r(mx);
end

% Dimensionless parameters entering the correction polynomials
miu_val = A / dk^2;          % Geometric ratio (→ 0 as d → ∞)
eta     = abs(w_tip_max / dk); % Normalised peak displacement (→ 0 as d → ∞)

% Polynomial regression (calibrated over dmin ≤ d ≤ 150 mm)
beta_raw  =  3.309  + 10.09*miu_val  - 9.176*eta  + 21.74*miu_val^2 ...
            - 23.47*miu_val*eta + 8.336*eta^2;
alpha_raw =  0.04007 - 39.13*miu_val + 2.1*eta    - 48.77*miu_val^2 ...
            + 68.64*miu_val*eta + 5.643*eta^2 - 12.06*eta^3;

end % eval_cyl_poly