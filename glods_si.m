function [glods_profile,Plist,flist,alfa,radius,func_eval] = glods_si(func_f,file_ini,x_ini,lbound,ubound)
%
% GLODS_SI - GLODS with Scale-Invariant (two-space) formulation
%
% This solver is a scale-invariant adaptation of the GLODS framework:
% all geometric operations (polling, distances, merging) are performed in
% normalized coordinates y in [0,1]^n, while objective values are evaluated
% in the original variables x in [lbound,ubound] via an affine map.
%
% Purpose:
%   Solve the bound constrained problem:
%        min f(x)  s.t.  lbound <= x <= ubound,
%   where x is a real vector of dimension n. Derivatives are not used.
%
% Input:
%   func_f   Objective f(x) defined in original space
%   file_ini Initialization file name (only used when list==4)
%   x_ini    Initial point in original space (only used when list==0)
%   lbound   Lower bounds (original space)
%   ubound   Upper bounds (original space)
%
% Output:
%   glods_profile Best-so-far record vs function evaluations
%   Plist         Approximations to local minimizers (returned in original space)
%   flist         Corresponding function values
%   alfa          Corresponding step size parameters
%   radius        Corresponding comparison radii
%   func_eval     Total number of NEW function evaluations performed
%
% -------------------------------------------------------------------------
% Provenance and attribution:
%   Based on the MATLAB reference implementation of GLODS:
%     A. L. Custodio and J. F. A. Madeira,
%     "GLODS: Global and Local Optimization using Direct Search",
%     Journal of Global Optimization, 62 (2015), 1-28.
%
%   Scale-invariant (two-space) reformulation introduced in:
%     J. F. A. Madeira,
%     "GLODS-SI: Scale-Invariant Global-Local Direct Search for
%      Engineering Design Optimization",
%     Journal of Computational Design and Engineering, 2026.
%     Manuscript ID JCDE-2026-065.
%
%   Copyright (C) 2026 J. F. A. Madeira.
%   SPDX-License-Identifier: LGPL-3.0-or-later
% -------------------------------------------------------------------------

tStart = tic;
format long e;
% Note: no global "warning off all" — silence specific warnings only if needed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Load parameters first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters_glods_si;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define SI map (x<->y) and evaluation wrapper f(y)=f_orig(y2x(y))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure column vectors
lbound_orig = lbound(:);
ubound_orig = ubound(:);

n_vars = length(lbound_orig);
if length(ubound_orig) ~= n_vars
    error('Error: lbound and ubound must have the same dimension.');
end

scale = ubound_orig - lbound_orig;
if any(scale == 0)
    error('Error: some variables have zero range (ubound == lbound). Handle fixed vars explicitly.');
end

% Normalized bounds in y-space
lbound_norm = zeros(n_vars,1);
ubound_norm = ones(n_vars,1);

% Map y -> x
y2x = @(y) lbound_orig + y .* scale;

% Objective wrapper: ALWAYS takes y, evaluates in x
func_f_orig = func_f;
func_f_eval = @(y) feval(func_f_orig, y2x(y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Initialization (Pini always in y)
%   - list==0 and list==4 normalization happens inside init_glods_si
%   - list in {1,2,3,5,6} generated in y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Pini, f_ini, alfa_list, radius_list, n, nPini, has_f_ini] = init_glods_si( ...
    list, file_ini, x_ini, lbound_norm, ubound_norm, ...
    user_list_size, nPini, lbound_orig, ubound_orig, scale);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Cache initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

func_eval = 0;

% Cache containers (used only when cache ~= 0; harmless when cache == 0).
CacheP     = [];
CachenormP = [];
CacheF     = [];

% Seed for stochastic init strategies (only affects list==1 LHS and list==2 random).
% rng('default') gives reproducible runs; switch to rng('shuffle') for non-reproducible.
if (list == 1) || (list == 2)
    rng('default');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Build initial lists (evaluate+merge unless list==4 provides f_ini)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Plist  = [];
flist  = [];
alfa   = [];
radius = [];
active = [];

% Initialize profile (pre-allocated to the evaluation budget).
% Falls back to 10000 if max_fevals is not set or non-positive.
if exist('max_fevals','var') && isnumeric(max_fevals) && max_fevals > 0
    glods_profile = zeros(1, max_fevals);
else
    glods_profile = zeros(1, 10000);
end

% We mimic original behaviour: iterate through Pini and MERGE sequentially.
% - If list==4: use f_ini/alfa_list/radius_list (no evaluations)
% - Otherwise: evaluate points via evaluate_with_cache (new evaluations)
%
% The while-loop remains to handle cases where all points are infeasible or f=Inf,
% especially for random/lds lists (1/2/5/6).

while isempty(flist)

    % For random/lds lists, Pini may be regenerated externally if needed.
    % (init_glods_si already creates Pini; but if everything fails, we can regenerate here)
    % NOTE: list==3 is deterministic in y; list==0 deterministic given x_ini.

    for i = 1:size(Pini,2)

        y_ini = Pini(:,i);

        % Feasibility in y
        feasible = is_feasible_y(y_ini, lbound_norm, ubound_norm);

        if feasible

            if has_f_ini
                % list==4: f_ini already known
                ftemp = f_ini(i);
            else
                % Evaluate with cache (in y)
                [ftemp, func_eval, CacheP, CachenormP, CacheF, ~] = evaluate_with_cache( ...
                    y_ini, func_f_eval, cache, CacheP, CachenormP, CacheF, tol_match, func_eval);
            end

            if isfinite(ftemp)

                if isempty(flist)
                    % Start lists
                    flist = ftemp;
                    Plist = y_ini;

                    if has_f_ini
                        alfa   = alfa_list(i);
                        radius = radius_list(i);
                    else
                        alfa   = alfa_ini;
                        radius = radius_ini;
                    end

                    active = true(1,1);

                    % Profile update (only meaningful when we have evaluations)
                    if func_eval >= 1
                        glods_profile(func_eval) = ftemp;
                    end

                else
                    % Merge subsequent points
                    if has_f_ini
                        alfa_aux   = alfa_list(i);
                        radius_aux = radius_list(i);
                    else
                        alfa_aux   = alfa_ini;
                        radius_aux = radius_ini;
                    end

                    [~,Plist,flist,alfa,radius,active,~] = merge( ...
                        y_ini, ftemp, alfa_aux, radius_aux, ...
                        Plist, flist, alfa, radius, active, ...
                        suf_decrease, 0, []);

                    if func_eval >= 1
                        if func_eval > length(glods_profile)
                            glods_profile(end+1:func_eval) = min(flist);
                        else
                            glods_profile(func_eval) = min(flist);
                        end
                    end
                end
            end
        end
    end

    % If still empty, handle failure cases
    if isempty(flist) && (list~=1) && (list~=2) && (list~=5) && (list~=6)
        fprintf('Error: The optimizer did not generate a feasible point\n');
        fprintf('or the initial point provided is not feasible.\n');
        fprintf('Please try list=1 or list=2 in parameters file.\n\n');
        return
    end

    if isempty(flist) && stop_feval && (func_eval >= max_fevals)
        fprintf('Error: The optimizer did not generate a feasible point,\n');
        fprintf('considering the budget of functions evaluations provided.\n\n');
        return
    end

    % If random/lds list and everything failed, regenerate Pini and retry
    if isempty(flist) && ((list==1)||(list==2)||(list==5)||(list==6))
        Pini = generate_Pini_glods_si(list, n, nPini, lbound_norm, ubound_norm);
    end
end

% Finalize profile if list==4 (no new evaluations)
if has_f_ini && func_eval == 0
    glods_profile = min(flist); % scalar baseline
else
    % truncate to current func_eval (later again at end)
    if func_eval > 0
        glods_profile = glods_profile(1:func_eval);
    else
        glods_profile = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: Set seed for poll directions and initialize counters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (dir_dense == 1)
    rng(1);
end

if search_size == 0
    search_size = n;
end

halt            = 0;
iter            = 0;
iter_suc        = 0;
unsuc_consec    = 0;
grid_size       = 1;
label_grid_size = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 7: Print header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if output
    fprintf('Iteration Report (GLODS_SI): \n\n');
    fprintf('| iter  | success | #active points |    min fvalue   |    min alpha    |    max alpha    |\n');
    print_format = ['| %5d |   %2s    |    %5d       | %+13.8e | %+13.8e | %+13.8e |\n'];

    active_indices = find(active);
    active_indices = active_indices(active_indices <= length(flist));
    active_indices = active_indices(active_indices <= length(alfa));

    if ~isempty(active_indices)
        min_f = min(flist(active_indices));
        min_a = min(alfa(active_indices));
        max_a = max(alfa(active_indices));
    else
        min_f = Inf; min_a = 0; max_a = 0;
    end

    fprintf(print_format, iter, '--', length(active_indices), min_f, min_a, max_a);

    fresult = fopen('glods_report.txt','w');
    fprintf(fresult,'Iteration Report (GLODS_SI): \n\n');
    fprintf(fresult,'| iter  | success | #active points |    min fvalue   |    min alpha    |    max alpha    |\n');
    fprintf(fresult,print_format, iter, '--', length(active_indices), min_f, min_a, max_a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP (Search + Poll)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (~halt)

    func_iter   = 0;
    aux_success = 0;
    success     = 0;
    poll        = 1;
    search      = 0;
    changes     = zeros(1,size(flist,2));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Search Step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if iter ~= 0 && search_option ~= 0
        if search_freq_type == 0
            if search_freq == 0
                search = 1;
            else
                if unsuc_consec == search_freq
                    search = 1;
                end
            end
        else
            index = find(alfa >= tol_active_points);
            aux_active = active(index);
            index = find(aux_active);
            if length(index) <= min_active_points
                search = 1;
            end
        end
    end

    finite = 0;
    if search && ~halt
        while ~finite
            unsuc_consec = 0;

            [Psearch,grid_size,label_grid_size] = search_step(search_option, ...
                search_size, lbound_norm, ubound_norm, grid_size, label_grid_size);

            if ~isempty(Psearch)
                for i = 1:size(Psearch,2)

                    ytemp = Psearch(:,i);

                    if is_feasible_y(ytemp, lbound_norm, ubound_norm)

                        [ftemp, func_eval, CacheP, CachenormP, CacheF, ~] = evaluate_with_cache( ...
                            ytemp, func_f_eval, cache, CacheP, CachenormP, CacheF, tol_match, func_eval);

                        func_iter = func_iter + 1;

                        if isfinite(ftemp)
                            finite = 1;

                            [success,Plist,flist,alfa,radius,active,changes] = merge( ...
                                ytemp, ftemp, alfa_ini, radius_ini, ...
                                Plist, flist, alfa, radius, active, ...
                                suf_decrease, 0, changes);

                            aux_success = aux_success + success;

                            % Profile
                            if func_eval > length(glods_profile)
                                glods_profile(end+1:func_eval) = min(flist);
                            else
                                glods_profile(func_eval) = min(flist);
                            end
                        end
                    end
                end

                if stop_feval && (func_eval >= max_fevals)
                    halt = 1;
                end

                if aux_success > 0
                    success = 1;
                    poll    = 0;
                else
                    success = 0;
                    poll    = 1;
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Poll Step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if poll && ~halt

        % Generate positive basis
        if (~dir_dense)
            D = [eye(n) -eye(n)];
        else
            v     = 2*rand(n,1)-1;
            [Q,R] = qr(v);
            if ( R(1) > 1 )
                D = Q * [ eye(n) -eye(n) ];
            else
                D = Q * [ -eye(n) eye(n) ];
            end
        end

        % Reorder points: by f then by alpha stop threshold then active first
        [flist,index] = sort(flist,'ascend');
        Plist         = Plist(:,index);
        alfa          = alfa(:,index);
        radius        = radius(:,index);
        active        = active(:,index);

        index1 = find(alfa >= tol_stop);
        index2 = find(alfa <  tol_stop);
        index1 = index1(index1 <= size(Plist,2));
        index2 = index2(index2 <= size(Plist,2));

        Plist  = [Plist(:,index1),Plist(:,index2)];
        flist  = [flist(index1),flist(index2)];
        alfa   = [alfa(index1),alfa(index2)];
        radius = [radius(index1),radius(index2)];
        active = [active(index1),active(index2)];

        active_indices   = find(active);
        inactive_indices = find(~active);

        active_indices   = active_indices(active_indices <= size(Plist,2));
        inactive_indices = inactive_indices(inactive_indices <= size(Plist,2));

        Plist  = [Plist(:,active_indices),Plist(:,inactive_indices)];
        flist  = [flist(active_indices),flist(inactive_indices)];
        alfa   = [alfa(active_indices),alfa(inactive_indices)];
        radius = [radius(active_indices),radius(inactive_indices)];
        active = [true(1,length(active_indices)), false(1,length(inactive_indices))];

        % Poll loop
        nd      = size(D,2);
        count_d = 1;
        changes = [1,zeros(1,length(flist)-1)];

        while ~success && (count_d <= nd)

            ytemp = Plist(:,1) + alfa(1) * D(:,count_d);

            if is_feasible_y(ytemp, lbound_norm, ubound_norm)

                [ftemp, func_eval, CacheP, CachenormP, CacheF, ~] = evaluate_with_cache( ...
                    ytemp, func_f_eval, cache, CacheP, CachenormP, CacheF, tol_match, func_eval);

                if isfinite(ftemp)
                    [success,Plist,flist,alfa,radius,active,changes] = merge( ...
                        ytemp, ftemp, alfa_ini, radius_ini, ...
                        Plist, flist, alfa, radius, active, ...
                        suf_decrease, 1, changes);

                    % Profile
                    if func_eval > length(glods_profile)
                        glods_profile(end+1:func_eval) = min(flist);
                    else
                        glods_profile(func_eval) = min(flist);
                    end
                end
            end

            count_d = count_d + 1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update step sizes and stopping criteria
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if success
        iter_suc     = iter_suc + 1;
        unsuc_consec = 0;

        update_indices = find(changes & active);
        update_indices = update_indices(update_indices <= length(alfa));
        if ~isempty(update_indices)
            alfa(update_indices)   = alfa(update_indices) * gamma_par;
            radius(update_indices) = max(radius(update_indices), alfa(update_indices));
        end
    else
        unsuc_consec = unsuc_consec + 1;

        update_indices = find(changes & active);
        update_indices = update_indices(update_indices <= length(alfa));
        if ~isempty(update_indices)
            alfa(update_indices) = alfa(update_indices) * beta_par;
        end
    end

    active_indices = find(active);
    active_indices = active_indices(active_indices <= length(alfa));
    if stop_alfa && ~isempty(active_indices) && (sum(alfa(active_indices) >= tol_stop) == 0)
        halt = 1;
    end
    if stop_feval && (func_eval >= max_fevals)
        halt = 1;
    end

    iter = iter + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print iteration report
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if output
        print_format = ['| %5d |    %1d    |    %5d       | %+13.8e | %+13.8e | %+13.8e |\n'];

        active_indices = find(active);
        active_indices = active_indices(active_indices <= length(flist));
        active_indices = active_indices(active_indices <= length(alfa));

        if ~isempty(active_indices)
            min_f = min(flist(active_indices));
            min_a = min(alfa(active_indices));
            max_a = max(alfa(active_indices));
        else
            min_f = Inf; min_a = 0; max_a = 0;
        end

        fprintf(print_format,iter,success,length(active_indices),min_f,min_a,max_a);
        fprintf(fresult,print_format,iter,success,length(active_indices),min_f,min_a,max_a);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print current active points (output == 2) in ORIGINAL space x
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (output == 2)
        fglods = fopen('glods_partial_results.txt','w');
        n_dim  = size(Plist,1);

        active_indices = find(active);
        active_indices = active_indices(active_indices <= size(Plist,2));
        active_indices = active_indices(active_indices <= length(flist));
        active_indices = active_indices(active_indices <= length(alfa));
        active_indices = active_indices(active_indices <= length(radius));

        m = length(active_indices);

        if m > 0
            fprintf(fglods,'%d %d\n\n', n_dim, m);

            % Build format string
            format_str = repmat(' %+21.16e', 1, m);
            format_str = [format_str, '\n'];

            % Denormalize active points to x
            Plist_orig_partial = zeros(n_dim, m);
            for k = 1:m
                Plist_orig_partial(:,k) = lbound_orig + Plist(:,active_indices(k)) .* scale;
            end

            % Write x points
            for j = 1:n_dim
                fprintf(fglods, format_str, Plist_orig_partial(j,:));
            end
            fprintf(fglods, '\n');
            fprintf(fglods, format_str, flist(active_indices));
            fprintf(fglods, '\n');
            fprintf(fglods, format_str, alfa(active_indices));
            fprintf(fglods, '\n');
            fprintf(fglods, format_str, radius(active_indices));
        else
            fprintf(fglods, '%d %d\n\n', n_dim, 0);
        end

        fclose(fglods);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL: Return active points in ORIGINAL space x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

active_indices_final = find(active);
active_indices_final = active_indices_final(active_indices_final <= size(Plist,2));
active_indices_final = active_indices_final(active_indices_final <= length(flist));
active_indices_final = active_indices_final(active_indices_final <= length(alfa));
active_indices_final = active_indices_final(active_indices_final <= length(radius));

if ~isempty(active_indices_final)
    Plist_y = Plist(:, active_indices_final);

    % Denormalize to original space x
    Plist = zeros(size(Plist_y));
    for k = 1:size(Plist_y,2)
        Plist(:,k) = lbound_orig + Plist_y(:,k) .* scale;
    end

    flist  = flist(active_indices_final);
    alfa   = alfa(active_indices_final);
    radius = radius(active_indices_final);
else
    Plist  = [];
    flist  = [];
    alfa   = [];
    radius = [];
end

% Truncate profile to actual number of evaluations
if func_eval > 0
    if length(glods_profile) >= func_eval
        glods_profile = glods_profile(1:func_eval);
    else
        glods_profile(end+1:func_eval) = min(flist);
    end
else
    % list==4 or no evaluation happened
    if ~isempty(flist)
        glods_profile = min(flist);
    else
        glods_profile = [];
    end
end

time = toc(tStart);

fprintf('\n Final Report (GLODS_SI - Scale-Invariant): \n\n');
fprintf('Elapsed Time = %10.3e \n\n', time);
fprintf('| #iter | #isuc | #active points | #fevals |    min fvalue   |\n');
fprintf('| %5d | %5d |    %5d       |  %5d  | %+13.8e |\n\n', ...
    iter, iter_suc, size(Plist,2), func_eval, min(flist));

if output
    fprintf(fresult,'\n Final Report (GLODS_SI): \n\n');
    fprintf(fresult,'Elapsed Time = %10.3e \n\n', time);
    fprintf(fresult,'| #iter | #isuc | #active points | #fevals |   min fvalue    |\n');
    fprintf(fresult,'| %5d | %5d |    %5d       |  %5d  | %+13.8e |\n\n', ...
        iter, iter_suc, size(Plist,2), func_eval, min(flist));
    fclose(fresult);
end

end % end of glods_si


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local helpers (can also be placed in separate .m files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ftemp, func_eval, CacheP, CachenormP, CacheF, match] = evaluate_with_cache( ...
    y, func_f_eval, cache, CacheP, CachenormP, CacheF, tol_match, func_eval)
% Evaluate f(y) using cache when enabled.

    match = 0;

    if cache ~= 0
        y_norm1 = norm(y,1);
        if ~isempty(CacheP)
            [match, y, ftemp] = match_point(y, y_norm1, CacheP, CacheF, CachenormP, tol_match);
        end
    end

    if ~match
        ftemp     = feval(func_f_eval, y);
        func_eval = func_eval + 1;

        if cache ~= 0
            CacheP     = [CacheP, y];
            CachenormP = [CachenormP, norm(y,1)];
            CacheF     = [CacheF, ftemp];
        end
    end
end


function feasible = is_feasible_y(y, lbound_norm, ubound_norm)
% Bound feasibility test in normalized space.
    y = y(:);
    lbound_norm = lbound_norm(:);
    ubound_norm = ubound_norm(:);

    if length(y) ~= length(lbound_norm) || length(y) ~= length(ubound_norm)
        error('Error: dimension mismatch in is_feasible_y.');
    end

    bound = [y - ubound_norm; -y + lbound_norm];
    feasible = (sum(bound <= 0) == 2*length(y));
end


function [Pini, f_ini, alfa_list, radius_list, n, nPini, has_f_ini] = init_glods_si( ...
    list, file_ini, x_ini, lbound_norm, ubound_norm, ...
    user_list_size, nPini, lbound_orig, ubound_orig, scale) %#ok<INUSD>
% Outputs:
%   Pini        - n x nPini matrix of initial points in normalized space [0,1]^n
%   f_ini       - 1 x nPini vector of function values (only for list==4)
%   alfa_list   - 1 x nPini vector of alfa values (only for list==4)
%   radius_list - 1 x nPini vector of radius values (only for list==4)
%   n           - problem dimension
%   nPini       - number of initial points (may differ from input nPini)
%   has_f_ini   - boolean indicating if f_ini is provided
%
% INIT_GLODS_SI
% Build initial list Pini in normalized space y in [0,1]^n.
%
% - list==0: x_ini is in ORIGINAL x (normalize x->y) or center (y). has_f_ini=false
% - list==4: file provides x in ORIGINAL space + f_ini/alfa_list/radius_list (normalize x->y). has_f_ini=true
% - list in {1,2,3,5,6}: generate Pini directly in y. has_f_ini=false
%
% Notes:
% - ubound_orig is used to recompute scale = ubound_orig - lbound_orig (robust).
% - Input "scale" is kept in signature for backward compatibility; it is not used.

    % Defaults
    f_ini       = [];
    alfa_list   = [];
    radius_list = [];
    has_f_ini   = false;

    % Dimension (trust normalized bounds for n)
    n = size(lbound_norm,1);

    % Defensive shaping (column vectors)
    lbound_orig = lbound_orig(:);
    ubound_orig = ubound_orig(:);

    if numel(lbound_orig) ~= n || numel(ubound_orig) ~= n
        error('Error: lbound_orig and ubound_orig must have %d elements.', n);
    end

    % Recompute scale for consistency
    scale = ubound_orig - lbound_orig;
    if any(scale == 0)
        error('Error: some variables have zero range (ubound == lbound). Handle fixed vars explicitly.');
    end

    % Decide nPini if user_list_size==0
    if (user_list_size == 0)
        nPini = n;
    end

    % -------------------------
    % list == 0 (single point)
    % -------------------------
    if (list == 0)
        if ~isempty(x_ini)
            if size(x_ini,1) ~= n
                error('Error: x_ini must have %d rows (one per variable).', n);
            end
            % Normalize ORIGINAL x_ini -> y
            Pini = (x_ini - lbound_orig) ./ scale;
            % Clamp to [0,1] (avoid tiny numerical violations)
            Pini = min(max(Pini, 0), 1);
        else
            Pini = (lbound_norm + ubound_norm)/2;
        end
        nPini = size(Pini,2);
        return;
    end

    % -------------------------
    % list == 4 (file in x + f/alfa/radius)
    % -------------------------
    if (list == 4)
        fpoints = fopen(file_ini,'r');
        if fpoints < 0
            error('Error: could not open initialization file: %s', file_ini);
        end
        c = onCleanup(@() fclose(fpoints)); %#ok<NASGU>

        aux = str2num(fgetl(fpoints)); %#ok<ST2NM>
        n_file = aux(1);
        m      = aux(2);

        if n_file ~= n
            error('Error: file dimension n=%d does not match bounds dimension n=%d.', n_file, n);
        end

        str2num(fgetl(fpoints)); %#ok<ST2NM> % keep file format (unused header line)

        x_file = zeros(n,m);
        for i = 1:n
            line = fgetl(fpoints);
            if ~ischar(line)
                error('Error: unexpected end-of-file while reading x_file (row %d).', i);
            end
            aux = str2num(line); %#ok<ST2NM>
            if numel(aux) < m
                error('Error: not enough entries in x_file row %d (expected %d).', i, m);
            end
            x_file(i,:) = aux(1:m);
        end

        % Normalize x -> y (implicit expansion in recent MATLAB)
        Pini = (x_file - lbound_orig) ./ scale;
        % Clamp to [0,1]
        Pini = min(max(Pini, 0), 1);

        % Read f_ini
        line = fgetl(fpoints);
        if ~ischar(line), error('Error: unexpected EOF before f_ini header.'); end %#ok<NASGU>
        aux  = str2num(fgetl(fpoints)); %#ok<ST2NM>
        if numel(aux) < m
            error('Error: not enough entries for f_ini (expected %d).', m);
        end
        f_ini = aux(1:m);

        % Read alfa_list
        line = fgetl(fpoints);
        if ~ischar(line), error('Error: unexpected EOF before alfa_list header.'); end %#ok<NASGU>
        aux  = str2num(fgetl(fpoints)); %#ok<ST2NM>
        if numel(aux) < m
            error('Error: not enough entries for alfa_list (expected %d).', m);
        end
        alfa_list = aux(1:m);

        % Read radius_list
        line = fgetl(fpoints);
        if ~ischar(line), error('Error: unexpected EOF before radius_list header.'); end %#ok<NASGU>
        aux  = str2num(fgetl(fpoints)); %#ok<ST2NM>
        if numel(aux) < m
            error('Error: not enough entries for radius_list (expected %d).', m);
        end
        radius_list = aux(1:m);

        nPini     = size(Pini,2);
        has_f_ini = true;
        return;
    end

    % ------------------------------------
    % list in {1,2,3,5,6}: generate in y
    % ------------------------------------
    Pini  = generate_Pini_glods_si(list, n, nPini, lbound_norm, ubound_norm);
    nPini = size(Pini,2); % important for list==3 (adds center)
end



function Pini = generate_Pini_glods_si(list, n, nPini, lbound_norm, ubound_norm)
% GENERATE_PINI_GLODS_SI - Generate initial points in normalized space y

    if (list == 1)
        % Latin Hypercube
        Pini = repmat(lbound_norm,1,nPini) + lhsdesign(nPini,n)' .* ...
               repmat((ubound_norm-lbound_norm),1,nPini);

    elseif (list == 2)
        % Random
        Pini = repmat(lbound_norm,1,nPini) + rand(n,nPini) .* ...
               repmat((ubound_norm-lbound_norm),1,nPini);

    elseif (list == 3)
        % Equally spaced points in the line segment + central point
        center = (lbound_norm + ubound_norm)/2;
        if (nPini <= 1)
            Pini = center;
        else
            t = linspace(0, 1, nPini);
            Pini = repmat(lbound_norm,1,nPini) + repmat(t,n,1) .* ...
                   repmat((ubound_norm-lbound_norm),1,nPini);
            Pini = [Pini, center];
        end

    elseif (list == 5)
        % Halton
        Lhalton = haltonset(n,'Skip',n+1);
        Pini = repmat(lbound_norm,1,nPini) + repmat((ubound_norm-lbound_norm),1,nPini) .* ...
               Lhalton(1:nPini,:)';

    elseif (list == 6)
        % Sobol
        Lsobol = sobolset(n,'Leap',2^n);
        Pini = repmat(lbound_norm,1,nPini) + repmat((ubound_norm-lbound_norm),1,nPini) .* ...
               Lsobol(1:nPini,:)';
    else
        Pini = [];
    end
end


function [success,Plist,flist,alfa,radius,active,changes] = merge(x,f,...
         alfa_ini,radius_ini,Plist,flist,alfa,radius,...
         active,suf_decrease,poll,changes)
%
% Purpose:
%
%    Function merge compares a new evaluated point with the current list of
%    points, deciding if it should be added to it and updating the list.
%
% Input:  
%
%         x (Point to be compared.)
%
%         f (Corresponding objective function value.)
%
%         alfa_ini (Initial step size.)
%
%         radius_ini (Initial radius of comparison.)
%
%         Plist (Current list of points.)
%
%         flist (Corresponding objective function values.)
%
%         alfa (Corresponding step sizes.)
%
%         radius (Corresponding radius of comparison.)
%
%         active (Corresponding point status.) 
%
%         suf_decrease (0-1 variable: 1 if the algorithm uses a
%                      globalization strategy based in imposing a 
%                      sufficient decrease condition; 0 otherwise.)
% 
%         poll (0-1 variable: 1 if merging is performed inside the poll
%              step, 0 otherwise.)
%
%         changes (Record of points under analysis.)
%
% Output: 
%
%         success (0-1 variable: 1 if a better point was found; 0 otherwise.)
%
%         Plist (Updated list of points.)
%
%         flist (Corresponding objective function values.)
%
%         alfa (Corresponding step sizes.)
%
%         radius (Corresponding radius of comparison.)
%
%         active (Corresponding point status.)
%
%         changes (Record of points analysed.)
%
% Functions called: forcing (Provided by the optimizer.)
%
% GLODS Version 0.2
%
% Copyright (C) 2014 A. L. Custódio and J. F. A. Madeira.
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
% License along with this library; if not, write to the Free Software
% Foundation, Inc.,51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
%
%
success   = 0;
dist_list = sqrt(sum((Plist-repmat(x,1,size(Plist,2))).^2,1));
if min(dist_list-radius)>0
   success = 1;
   Plist   = [Plist,x];
   flist   = [flist,f];
   alfa    = [alfa,alfa_ini];
   radius  = [radius,radius_ini];
   active  = [active,1];
   changes = [changes,1];  
else
    if min(dist_list) ~= 0
       index      = find(dist_list-radius<=0);
       m_index    = size(index,2);
       active_new = 0;
       alfa_new   = 0;
       radius_new = 0;
       idom       = 0;
       pdom       = 0;
       icomp      = 0;
       for i=1:m_index
           if f < flist(index(i)) - forcing(alfa(index(i)),suf_decrease)
              icomp            = 1;
              idom             = idom + active(index(i));
              active(index(i)) = 0;
              if alfa(index(i)) > alfa_new
                 alfa_new   = alfa(index(i));
                 radius_new = radius(index(i));
              end   
           else
              if flist(index(i)) <= f - forcing(alfa(index(i)),suf_decrease)
                 pdom = 1;
              end
           end
       end
       if pdom == 0
           active_new = 1;
           success    = 1;
       end
       if (idom > 0) || (suf_decrease == 0 && pdom == 0 && icomp == 1)
          Plist   = [Plist,x];
          flist   = [flist,f];
          active  = [active,active_new];
          changes = [changes,1];
          if poll
             alfa   = [alfa,alfa(1)];
             radius = [radius,alfa(1)];  
          else
             alfa   = [alfa,alfa_new];
             radius = [radius,radius_new];    
          end
       end
    end
end
%
% End of merge.
end

function [rho] = forcing(alfa,suf_decrease)
%
% Purpose:
%
%    Function forcing implements a globalization strategy based on
%    imposing a sufficient decrease condition.
%
% Input:  
%
%         alfa (Step size, given by the optimizer.)
%
%         suf_decrease (0-1 variable: 1 if the algorithm uses a
%                      globalization strategy based in imposing a 
%                      sufficient decrease condition; 0 otherwise.)
%
% Output: 
%
%         rho (Forcing function value at the given step size.)        
%
% GLODS Version 0.1
%
% Copyright (C) 2013 A. L. Custódio and J. F. A. Madeira.
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
% License along with this library; if not, write to the Free Software
% Foundation, Inc.,51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
%
%
if suf_decrease
   rho = alfa^2;
else
   rho = 0;
end
%
% End of forcing.
end

function [match,x,f] = match_point(x,xnorm,CacheP,CacheF,CachenormP,tol_match);
%
% Purpose:
%
%    Function match_point scans a list of previously evaluated points,
%    trying to match a point provided by the optimizer. When matching 
%    is successful, the function returns the corresponding objective
%    function value.
%
% Input: 
%
%         x (Point to be checked.)
%
%         xnorm (1-norm of the point to be checked.)
%
%         CacheP (Matrix of points to be used in comparisons,
%                storing the points columnwise.)
%
%         CacheF (Matrix storing the objective function values
%                corresponding to the points in CacheP.)
%
%         CachenormP (Vector storing the 1-norm of the points in CacheP.)
%
%         tol_match (Tolerance value within which two points are
%                   considered as equal.)
%
% Output: 
%
%         match (0-1 variable: 1 if x was previously evaluated; 0
%         otherwise.)
%
%         x (Vector storing the matched point coordinates.)
%
%         f (Objective function value at the matched point.)
%
% GLODS Version 0.1
%
% Copyright (C) 2013 A. L. Custódio and J. F. A. Madeira.
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
% License along with this library; if not, write to the Free Software
% Foundation, Inc.,51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
%
%
% Prune the points that do not satisfy a 1-norm criterion.
%
f     = [];
index = find(abs(CachenormP - xnorm) <= tol_match);
if ~isempty(index)
   CacheP = CacheP(:,index);
   CacheF = CacheF(:,index);
%
% Finish search.
%
   nCacheP = size(CacheP,2);
   index   = find(max(abs(CacheP-repmat(x,1,nCacheP)),[],1) <= tol_match);
end
match = ~isempty(index);
%
% Retrieve the point coordinates and the objective function values. 
%
if match
   x = CacheP(:,index(1));
   f = CacheF(:,index(1));
end
%
% End of match_point.

function  [M_new] = replicate(v,M)
%
% Purpose:
%
%    Auxiliary function that computes a new matrix by duplicating a 
%    provided matrix and duplicating a provided vector.
%
% Input:  
%
%         v (Vector to be duplicated.)
%
%         M (Matrix to be duplicated.)
%
% Output: 
%
%         M_new (New matrix.)        
%
% GLODS Version 0.1
%
% Copyright (C) 2013 A. L. Custódio and J. F. A. Madeira.
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
% License along with this library; if not, write to the Free Software
% Foundation, Inc.,51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
%
%
M_new = [];
   for i = 1:numel(v)
    M_new = [M_new; M repmat(v(i),size(M,1),1)];
   end
end
%
% End of replicate.
end


function [Psearch,grid_size,label_grid_size] = search_step(search_option,...
                    search_size,lbound,ubound,grid_size,label_grid_size)
%     
% Purpose:
%   
%    Function search_step generates a new set of points, to be used in the
%    search step, which should be asymptotically dense in the feasible region.
%    
% Input:  
%
%         search_option(0-5 variable: 0 if no search step should be performed;
%                      1 if a latin hypercube sampling strategy is considered; 
%                      2 if random sampling is used; 3 if points are uniformly 
%                      spaced in a grid (2^n-Centers strategy); 4 if points 
%                      are generated using Halton sequences; 5 if generation 
%                      is based on Sobol numbers.)
%
%         search_size (Number of points to be generated at each search
%                     step, when search_option is not set equal to 3. If 
%                     search_option is set equal to 3 then search_size is
%                     equal to 2^(problem dimension).) 
%
%         lbound (Lower bound on the problem variables.)
%
%         ubound (Upper bound on the problem variables.)
%
%         grid_size (Level of the current grid when search_option is set 
%                   equal to 3.)
%
%         label_grid_size (Label that indicates the current point to be
%                         used in the current grid when search_option is
%                         set equal to 3.)
%
% Output: 
%
%         Psearch (New list of points to be evaluated in the search step.)
%
%         grid_size (Level of the last used grid when search_option is set 
%                   equal to 3.)
%
%         label_grid_size (Label that indicates the next point to be used 
%                         in the grid when search_option is set equal to 3.)
%
% Functions called: replicate (Provided by the optimizer).
%
%
% GLODS Version 0.1
%
% Copyright (C) 2013 A. L. Custódio and J. F. A. Madeira.
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
% License along with this library; if not, write to the Free Software
% Foundation, Inc.,51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
%
%
Psearch = [];
n       = size(lbound,1);
%
if (search_option == 1)
   Psearch = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
             lhsdesign(search_size,n)';
end
%
if (search_option == 2)
   Psearch = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
             rand(n,search_size);
end
%
if (search_option == 3)
    vaux       = [1/2^grid_size:1/2^(grid_size-1):1-1/2^grid_size]';
    Psearch    = vaux;
    for j = 1:n-1
       Psearch = replicate(vaux,Psearch);
    end
    if grid_size <= 2
       if grid_size == 1
          pindex = [1]; 
       else
          pindex = [1:2^n];
       end
       grid_size       = grid_size + 1;
       label_grid_size = 1;
    else
      index = [1]; 
      level = 0;
      count = size(index);
      while count < (2^n)^(grid_size-2)
         index_old = index;
         for i = 1:(2^(grid_size-2)-1)
           index_aux = (2^(grid_size-1))^level*i+index_old;
           index     = [index, index_aux];
         end
         count = size(index);
         level = level +1;
      end
      pindex = [index(label_grid_size) index(label_grid_size)+2^(grid_size-2)];
      level  = 1;
      while size(pindex) < 2^n  
        pindex = [pindex pindex+(2^(grid_size-1))^level*2^(grid_size-2)];
        level  = level +1;
      end
      if label_grid_size == (2^n)^(grid_size-2)
         grid_size       = grid_size + 1;
         label_grid_size = 1;
      else
         label_grid_size = label_grid_size + 1; 
      end
    end
    Psearch = Psearch(pindex,:)';
    if (grid_size == 2)
        vaux        = [0 1]';
        Psearch_aux = vaux;
        for j = 1:n-1
           Psearch_aux = replicate(vaux,Psearch_aux);
        end
        Psearch = [Psearch, Psearch_aux'];
    end
    Psearch = Psearch.*repmat((ubound-lbound),1,size(Psearch,2))+...
              repmat(lbound,1,size(Psearch,2));
end
%
if (search_option == 4)
    Lhalton   = haltonset(n,'Skip',n+1);
    Lhalton   = scramble(Lhalton,'RR2');
    Psearch   = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
                Lhalton(grid_size:grid_size + search_size - 1,:)';
    grid_size = grid_size + search_size;
end
%
if (search_option == 5)
    Lsobol    = sobolset(n,'Leap',2^n);
    Psearch   = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
                Lsobol(grid_size:grid_size + search_size - 1,:)';
    grid_size = grid_size + search_size;
end
%
% End of search_step.
end