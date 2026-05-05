function [alpha,beta,dmin]=correct_cuboidalmagnet_coeff(MA,W,H,L,d)
% CORRECT_CUBOIDALMAGNET_COEFF Computes magnetic force correction coefficients
% for cuboidal magnets with square cross-section by introducing α and β.
%
%   [ALPHA,BETA,DMIN] = CORRECT_CUBOIDALMAGNET_COEFF(MA,W,H,L,d) calculates
%   correction coefficients α and β for magnetic force modeling between
%   two identical face-to-face cuboidal magnets sharing magnetization MA.
%
%--------------------------------------------------------------------------
% Workflow:
%   1. Polynomial region  (dmin < d ≤ 150 mm):
%        α, β from least-squares polynomial fit over calibration dataset.
%        Exact analytical Coulombian force model used internally to locate
%        the transverse force peak (→ dimensionless parameters μ, η).
%        Special case: when W = L, A_sq = W²−L² = 0 → μ ≡ 0. The polynomial
%        then depends only on η, and its d→∞ asymptote deviates from 1.
%        Exponential extrapolation is especially important in this case.
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
%   2W=2H : cross-section side length, 5~20 mm (square cross-section)
%   2L    : total magnet thickness,    2~20 mm
%   d     : inter-magnet spacing,      dmin ~ ∞
%--------------------------------------------------------------------------
% Exact force model used internally (analytical Coulombian charge model):
%   Fmy = μ₀/(4π)·M²·Σ(±) f(corner coordinates)
%   Assembled from logarithmic and arctangent primitives via triple
%   integration of the magnetic scalar potential (8-corner superposition).
%--------------------------------------------------------------------------
% Inputs:
%   MA : Magnetization magnitude (A/m)
%   W  : Magnet half-width  (m) [square cross-section: W = H]
%   H  : Magnet half-height (m)
%   L  : Magnet half-length (m) [half-thickness along magnetization axis]
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
d_trans = 150e-3;   % Polynomial upper limit (m); extrapolation beyond
tol_fit = 1e-5;     % Min deviation threshold for λ regression

% ===== MINIMUM VALID SPACING =====
dmin = 1.587*sqrt(max(0, H^2-L^2)) - 0.5577*2*H + 1.589*2*L;

% ===== PRE-COMPUTE REFERENCE VALUES AND FIT λ =====
d_refs_all = [30, 60, 90, 120, 150] * 1e-3;
valid_ref  = (d_refs_all > dmin) & (d_refs_all <= d_trans);

if sum(valid_ref) < 1
    warning('correct_cuboidalmagnet_coeff: no valid reference points for lambda fitting.');
    lam_a = 0;
    lam_b = 0;
    [alpha_trans, beta_trans] = eval_cub_poly(d_trans, W, H, L, MA);

else
    d_refs     = d_refs_all(valid_ref);
    alpha_refs = zeros(1, length(d_refs));
    beta_refs  = zeros(1, length(d_refs));

    for ki = 1:length(d_refs)
        [alpha_refs(ki), beta_refs(ki)] = eval_cub_poly(d_refs(ki), W, H, L, MA);
    end

    % Anchor values at d_trans
    if abs(d_refs(end) - d_trans) < 1e-10
        alpha_trans = alpha_refs(end);
        beta_trans  = beta_refs(end);
    else
        [alpha_trans, beta_trans] = eval_cub_poly(d_trans, W, H, L, MA);
        d_refs     = [d_refs,     d_trans     ];
        alpha_refs = [alpha_refs, alpha_trans ];
        beta_refs  = [beta_refs,  beta_trans  ];
    end

    %----------------------------------------------------------------------
    % Fit λ via origin-constrained log-space linear regression (inlined)
    %
    % Model: f(d) = 1 + (f_trans−1)·exp(−λ·(d−d_trans))
    % Let x_i = d_i−d_trans (<0),  y_i = log(|dev_i/dev_trans|) (>0):
    %   λ = −Σ(x_i·y_i) / Σ(x_i²)   [origin-constrained OLS]
    %----------------------------------------------------------------------

    % ---- λ for alpha ----
    dev_trans_a = alpha_trans - 1;
    if abs(dev_trans_a) < tol_fit
        lam_a = 0;
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

    elseif d(k) - 2*L <= 0
        %------------------------------------------------------------------
        % d ≤ dmin: model inapplicable
        %------------------------------------------------------------------
        alpha(k) = NaN;
        beta(k)  = NaN;

    else
        %------------------------------------------------------------------
        % Polynomial region (dmin < d ≤ 150 mm)
        %------------------------------------------------------------------
        [alpha(k), beta(k)] = eval_cub_poly(d(k), W, H, L, MA);

    end

end % parfor
toc;

end % main function


%==========================================================================
% LOCAL FUNCTION: eval_cub_poly
%   Evaluates polynomial correction at a single spacing dk.
%   Computes exact analytical Coulombian Fmy(w), locates peak, derives μ, η.
%==========================================================================
function [alpha_raw, beta_raw] = eval_cub_poly(dk, W, H, L, MA)

num  = 12001;
w0   = linspace(-60e-3, 60e-3, num);
c60  = cos(pi);   % = −1  (face-to-face configuration: rotation angle θ = π)
s60  = sin(pi);   % =  0

A0 = W;  B0 = H;  C0 = L;
d0 = dk - 2*L;    % Surface-to-surface gap (m)

a_val = 2*W;   b_val = 2*H;   c_val = 2*L;
A_val = 2*A0;  B_val = 2*B0;  C_val = 2*C0;
z01_val = d0 + 2*C0 + 2*L;

% ---- Exact transverse force Fmy(w) via 8-corner Coulombian model ----
% Force assembled from 8 corner contributions via inclusion-exclusion:
%   Fmy = μ₀/(4π)·M²·(t1−t2+t3−t4+t5−t6+t7−t8)
Fmy = zeros(1, num);
for m = 1:num
    y01 = 2*B0 + w0(m);
    t1 = calc_term(0,      A_val,         y01,z01_val,c60,s60,0,    0,    B_val,C_val,MA);
    t2 = calc_term(-a_val,-a_val+A_val,   y01,z01_val,c60,s60,0,    0,    B_val,C_val,MA);
    t3 = calc_term(-a_val,-a_val+A_val,   y01,z01_val,c60,s60,b_val,0,    B_val,C_val,MA);
    t4 = calc_term(0,      A_val,         y01,z01_val,c60,s60,b_val,0,    B_val,C_val,MA);
    t5 = calc_term(-a_val,-a_val+A_val,   y01,z01_val,c60,s60,0,    c_val,B_val,C_val,MA);
    t6 = calc_term(0,      A_val,         y01,z01_val,c60,s60,0,    c_val,B_val,C_val,MA);
    t7 = calc_term(0,      A_val,         y01,z01_val,c60,s60,b_val,c_val,B_val,C_val,MA);
    t8 = calc_term(-a_val,-a_val+A_val,   y01,z01_val,c60,s60,b_val,c_val,B_val,C_val,MA);
    Fmy(m) = t1-t2+t3-t4+t5-t6+t7-t8;
end

[~, mx] = max(Fmy);
w0_max  = w0(mx);

% Dimensionless parameters
% Note: A_sq = W²−L² = 0 when W = L (degenerate geometry, miu ≡ 0)
A_sq  = A0^2 - C0^2;
miu   = A_sq / dk^2;
eta   = abs(w0_max / dk);

% Polynomial regression (calibrated over dmin ≤ d ≤ 150 mm)
beta_raw  =  3.166  - 4.637*miu  - 12.32*eta  - 10.46*miu^2 ...
            + 4.232*miu*eta + 26.07*eta^2 - 22.36*eta^3 + 16.04*miu*eta^2;
alpha_raw =  0.3065 - 22.94*miu  - 4.872*eta  - 25.12*miu^2 ...
            + 33.56*miu*eta + 34.36*eta^2 - 44.23*eta^3 + 27.88*miu*eta^2;

end % eval_cub_poly


%==========================================================================
% HELPER: calc_term
%   One signed corner contribution: μ₀/(4π)·M² prefactor × double-difference
%==========================================================================
function val = calc_term(v,w,y01,z01,ctheta,stheta,b,c,B,C,M)
val = 1e-7*M^2*(f2(v,w,y01,z01,ctheta,stheta,b,c,B,C) - ...
                f2(v,w,y01,z01,ctheta,stheta,b,c,B,0));
end

%==========================================================================
% HELPER: f2
%   Double difference over B-dimension (inclusion-exclusion over magnet 2
%   y-boundaries).
%==========================================================================
function res = f2(v,w,y01,z01,ctheta,stheta,b,c,B,z0)
res = f3(w,y01,z01,ctheta,stheta,b,c,B,z0) - ...
      f3(v,y01,z01,ctheta,stheta,b,c,B,z0) - ...
      f3(w,y01,z01,ctheta,stheta,b,c,0,z0) + ...
      f3(v,y01,z01,ctheta,stheta,b,c,0,z0);
end

%==========================================================================
% HELPER: f3
%   Core analytical primitive of the Coulombian force integral.
%   Evaluates closed-form antiderivative at corner (u,y0,z0) after rotation
%   by angle θ. Combines logarithmic, arctangent, and linear terms from
%   triple integration of the magnetic scalar potential.
%==========================================================================
function res = f3(u,y01,z01,ctheta,stheta,b,c,y0,z0)
% Rotated coordinate offsets
f6 =  y01*ctheta + z01*stheta - b*ctheta - c*stheta + y0;
f5 = -y01*stheta + z01*ctheta + b*stheta - c*ctheta + z0;

% Distance (safeguard against negative radicand at degenerate geometry)
f4_sq = u^2 + f5^2 + f6^2;
if f4_sq < 0; f4 = 0; else; f4 = sqrt(f4_sq); end

% Logarithmic term 1: u·f6·ln(−u+r)
log_arg1 = -u + f4;
if log_arg1 <= 0; part1 = 0;
else;             part1 = u*f6*log(log_arg1); end

% Logarithmic terms 2 & 6: involve ln(r+f6)
log_arg2 = f4 + f6;
if log_arg2 <= 0
    part2 = 0;  part6 = 0;
else
    part2 = -u^2*log(log_arg2);
    part6 =  0.5*(u^2+f5^2)*log(log_arg2);
end

% Arctangent term (solid-angle contribution; singularity-protected)
denom = f5*f6;
if abs(denom) < 1e-10
    part3 = 0;
else
    part3 = u*f5*atan((-f5^2 - u^2 + u*f4) / denom);
end

% Sign-dependent and direct spatial terms
part4 = 0.5*u*pi*abs(f5)*sign(f6);
part5 = 0.5*f6*f4;

res = part1 - u*f6 + part2 + part3 + part4 + part5 + part6;
end % f3