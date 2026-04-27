% parameters_glods_si.m script file
%
% Purpose:
%
% File parameters_glods_si sets the algorithmic strategies, parameters,
% constants, and tolerances values to be used by the function glods_si,
% which can be user modified.
%
% These are the parameters used to produce the numerical experiments
% reported in the accompanying paper (JCDE-2026-065).
%
% -------------------------------------------------------------------------
% Provenance and attribution:
%   Based on the parameters file of the GLODS reference implementation:
%     GLODS Version 0.2
%     Copyright (C) 2014 A. L. Custodio and J. F. A. Madeira.
%
%   Modified for the scale-invariant (two-space) formulation introduced in:
%     J. F. A. Madeira,
%     "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%      Engineering Design Optimization",
%     Journal of Computational Design and Engineering, 2026.
%     Manuscript ID JCDE-2026-065.
%
%   Modifications copyright (C) 2026 J. F. A. Madeira.
%   SPDX-License-Identifier: LGPL-3.0-or-later
% -------------------------------------------------------------------------
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 3.0 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, see <https://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
output = 2; % 0-2 variable: 0 if only a final report is displayed in the
            % screen; 1 if at each iteration output is displayed in the
            % screen and recorded in a text file stored at the current directory;
            % 2 level of output similar to 1 but additionally, at each
            % iteration, current approximations to local minimizers are
            % recorded in a file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stopping Criteria.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
stop_alfa = 1;     % 0-1 variable: 1 if the stopping criterion is based
                   % on the step size parameter; 0 otherwise.
tol_stop  = 10^-8; % Lowest value allowed for the step size parameter.
%
stop_feval = 1;     % 0-1 variable: 1 if the stopping criterion is based
                    % on a maximum number of function evaluations; 0
                    % otherwise.
max_fevals = 20000; % Maximum number of function evaluations allowed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algorithmic Options.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
suf_decrease = 0; % 0-1 variable: 1 if the algorithm uses a globalization
                  % strategy based on imposing a sufficient decrease
                  % condition; 0 otherwise.
%
% Cache Use.
%
cache     = 0;        % 0-1 variable: 0 if point evaluation is always done;
                      % 1 if a cache is maintained.
tol_match = tol_stop; % Tolerance used in point comparisons, when
                      % considering a cache.
%
% Initialization.
%
list = 6; % 0-6 variable: 0 if the algorithm initializes the list of
          % points with a single one; 1 if a latin hypercube sampling
          % strategy is considered for initialization; 2 if
          % random sampling is used; 3 if points are considered
          % equally spaced in a line segment, joining the
          % variable upper and lower bounds, jointly with
          % the central point; 4 if the algorithm is
          % initialized with a list provided by the user;
          % 5 if Halton sequences; 6 if Sobol numbers.
%
user_list_size = 1;  % 0-1 variable: 1 if the user sets below the size of
                     % the initial list of points; 0 if the initial size
                     % equals the problem dimension.
nPini          = 30; % Number of points to be considered in the initialization.
%
% Search step.
%
search_option = 5; % 0-5 variable: 0 if no search step should be performed;
                   % 1 if a latin hypercube sampling strategy is considered;
                   % 2 if random sampling is used; 3 if points are uniformly
                   % spaced in a grid (2^n-Centers strategy);
                   % 4 if points are generated using Halton sequences;
                   % 5 if generation is based on Sobol numbers.
%
search_size   = 0; % Number of points to be generated at each search step.
                   % If search_size is set equal to 0 (which should be the
                   % case when search_option is set equal to 3) then
                   % search_size is equal to problem dimension.
%
search_freq_type = 0; % 0-1 variable: 0 if the execution of the search step
                      % is based on the frequency of unsuccessful iterations;
                      % 1 if the execution depends on the number of current
                      % active points not yet identified as local minimizers.
%
search_freq = 10; % Number of consecutive unsuccessful iterations, after
                  % which the search step should be performed. Only used
                  % when search_freq_type is set equal to 0.
%
min_active_points = 1; % Number of active points for which a search step is
                       % performed. Only used when search_freq_type is set
                       % equal to 1.
%
tol_active_points = tol_stop; % Tolerance value used for identifying a point
                              % as a good candidate to a local minimizer.
                              % Only used when search_freq_type is set equal
                              % to 1.
%
% Directions and step size.
%
dir_dense = 0; % 0-1 variable: 1 if an asymptotically dense set of
               % directions should be considered for polling; 0 otherwise.
%
% NOTE on alfa_ini and radius_ini (GLODS-SI):
%   In normalized coordinates, alfa_ini and radius_ini are dimensionless
%   and represent fractions of each variable range. The defaults below
%   match the convention used to generate Table 3 of the accompanying
%   paper: an initial step of 10% and a comparison radius of 20% of
%   each variable range. The condition r_0 >= d_max * alfa_ini is
%   satisfied with d_max = 1 for the positive basis [I -I], so nearby
%   initial points are not merged immediately.
%
alfa_ini   = 0.1;             % Initial step size in normalized space.
radius_ini = 2 * alfa_ini;    % Initial comparison radius.
%
beta_par = 0.5; % Coefficient for step size contraction.
%
gamma_par = 2;  % Coefficient for step size expansion. In case of suf_decrease, we
                % recommend setting gamma_par equal to 1.
%
% End of parameters_glods_si.
