function varargout = aluffi_pentini_2D(varargin)
%ALUFFI_PENTINI_2D  Aluffi-Pentini 2D test problem (heterogeneous WORK-space wrapper).
%
% INPUT SPACE (illustrative heterogeneity):
%
%   x1   in [-10  ,  10 ]   (range:  20)
%   x2   in [-100 , 100 ]   (range: 200)
%
% Effective contrast ratio (max range / min range): 10
%
% Known global minimum:
%   x* = [-1.0465; 0]
%   f* = -0.3523
%
% USAGE:
%   f        = aluffi_pentini_2D(x)   % Evaluate at a 2D point x
%   [lb, ub] = aluffi_pentini_2D(n)   % Get bounds (n must equal 2)
%   info     = aluffi_pentini_2D()    % Get problem metadata
%
% References:
%   Aluffi-Pentini test function:
%     Ali, M. M., Khompatraporn, C., Zabinsky, Z. B. (2005).
%     A numerical evaluation of several stochastic algorithms on
%     selected continuous global optimization test problems.
%     Journal of Global Optimization, 31:635-672.
%
%   Heterogeneous-bounds wrapper formulation:
%     J. F. A. Madeira,
%     "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%      Engineering Design Optimization",
%     Journal of Computational Design and Engineering, 2026.
%     Manuscript ID JCDE-2026-065.
%
% Copyright (C) 2026 J. F. A. Madeira.
% SPDX-License-Identifier: LGPL-3.0-or-later

nloc           = 2;
lb_orig        = [-10; -10];
ub_orig        = [ 10;  10];
lb_work        = [-10; -100];
ub_work        = [ 10;  100];
scale_factors  = [1; 10];
contrast_ratio = 10;

if nargin == 0
    info.name              = mfilename;
    info.problem           = 'aluffi_pentini';
    info.dimension         = nloc;
    info.strategy          = 'illustrative';
    info.kappa             = contrast_ratio;
    info.lb_orig           = lb_orig;       info.ub_orig = ub_orig;
    info.lb_work           = lb_work;       info.ub_work = ub_work;
    info.scale_factors     = scale_factors;
    info.contrast_ratio    = contrast_ratio;
    info.global_min_known  = true;
    info.f_global_min      = -0.3523;
    info.x_global_min_orig = [-1.0465; 0];
    info.x_global_min_work = [-1.0465; 0];
    info.global_min_note   = ['Aluffi-Pentini (2D): single global minimum at ', ...
                              'x*=(-1.0465, 0), f*=-0.3523. Two local minima ', ...
                              'are also present. Ref: Ali et al. (2005).'];
    info.mapping           = 'x_orig = lb_orig + clip01((x-lb_work)/(ub_work-lb_work)).*(ub_orig-lb_orig)';
    varargout{1}           = info;
    return
end

arg1 = varargin{1};
if isscalar(arg1) && arg1 == round(arg1)
    if arg1 ~= nloc, error('This instance is 2D only.'); end
    varargout{1} = lb_work;
    if nargout >= 2, varargout{2} = ub_work; end
    return
end

x = arg1(:);
if numel(x) ~= nloc
    error('Input x must have 2 components.');
end
range            = ub_work - lb_work;
range(range==0)  = 1;
t                = (x - lb_work) ./ range;
t                = max(0, min(1, t));
x_orig           = lb_orig + t .* (ub_orig - lb_orig);
varargout{1}     = aluffi_pentini_orig(x_orig);
return

% -------------------------------------------------------------------------
% Embedded original objective (verbatim).
% -------------------------------------------------------------------------
function f = aluffi_pentini_orig(x)
%
%    Aluffi-Pentini function, n = 2.
%    f(x) = 0.25*x1^4 - 0.5*x1^2 + 0.1*x1 + 0.5*x2^2
%
%    Minimum global value : -0.3523
%    Global minimum        : x = (-1.0465, 0)
%    Local minima          : two
%    Search domain         : -10 <= x_i <= 10
%
%    Reference: Ali et al. (2005).
%
% Written by A. L. Custodio and J. F. A. Madeira.
% Version June 2012.
%
f = 0.25*x(1)^4 - 0.5*x(1)^2 + 0.1*x(1) + 0.5*x(2)^2;
%
% End of aluffi_pentini.
