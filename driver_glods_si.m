% driver_glods_si.m
% =========================================================================
% Example driver for GLODS_SI (Scale-Invariant GLODS).
%
% This script runs GLODS_SI on the Aluffi-Pentini 2D test problem with
% bounds [-10, 10] x [-100, 100], which gives a 10x range contrast
% between the two coordinates. The aluffi_pentini_2D wrapper handles
% the y -> x mapping internally; the solver works in normalized
% coordinates y in [0, 1]^n and reports back in the work-space
% coordinates passed to it.
%
% Companion code to:
%   J. F. A. Madeira,
%   "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%    Engineering Design Optimization",
%   Journal of Computational Design and Engineering, 2026.
%   Manuscript ID JCDE-2026-065.
%
% Copyright (C) 2026 J. F. A. Madeira.
% SPDX-License-Identifier: LGPL-3.0-or-later
% =========================================================================

clear; clc;
format compact;

% -------------------------------------------------------------------------
% Bounds for the Aluffi-Pentini wrapper (work space).
% -------------------------------------------------------------------------
lbound = [-10; -100];
ubound = [ 10;  100];

% -------------------------------------------------------------------------
% Run GLODS_SI.
% -------------------------------------------------------------------------
[glods_profile, Plist, flist, alfa, radius, func_eval] = ...
    glods_si(@aluffi_pentini_2D, [], [], lbound, ubound);

% -------------------------------------------------------------------------
% Display basic results.
% -------------------------------------------------------------------------
fprintf('\nGLODS_SI finished.\n');
fprintf('Number of function evaluations: %d\n', func_eval);

if ~isempty(flist)
    [fbest, ibest] = min(flist);
    xbest          = Plist(:, ibest);

    fprintf('Best objective value: %.16e\n', fbest);
    fprintf('Best point (work-space coordinates):\n');
    disp(xbest.');
end

% -------------------------------------------------------------------------
% Optional: convergence profile (best-so-far value vs. function evaluations).
% -------------------------------------------------------------------------
if ~isempty(glods_profile)
    figure;
    plot(glods_profile, 'LineWidth', 1.5);
    grid on;
    xlabel('Function evaluations');
    ylabel('Best objective value');
    title('GLODS\_SI convergence profile (Aluffi-Pentini 2D)');
end

% =========================================================================
% Using GLODS_SI with your own objective function
% =========================================================================
%
% The first argument to glods_si is a function handle returning a scalar
% f(x), where x is a column vector of length n in the original variable
% space [lbound, ubound]. Three common usage patterns:
%
%   (a) Anonymous function:
%
%         f      = @(x) sum(x.^2);
%         lbound = [-5; -5; -5];
%         ubound = [ 5;  5;  5];
%         [profile, Plist, flist, alfa, r, fevals] = ...
%             glods_si(f, [], [], lbound, ubound);
%
%   (b) Function defined in its own file my_objective.m:
%
%         [profile, Plist, flist, alfa, r, fevals] = ...
%             glods_si(@my_objective, [], [], lbound, ubound);
%
%   (c) Self-contained wrapper with an internal y -> x mapping (as in
%       aluffi_pentini_2D.m used above), recommended for problems whose
%       variables span heterogeneous physical scales.
%
% Algorithmic options (initial step size, search strategy, evaluation
% budget, etc.) are set in parameters_glods_si.m.
%
% =========================================================================
% End of driver_glods_si.m
% =========================================================================
