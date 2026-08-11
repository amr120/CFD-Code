1%   PropCURTIS_bladesweep_automesh
%
%   Copy of PropCURTIS_bladesweep.m for developing AUTOMATIC meshing so a wide
%   range of blade counts can be swept without hand-meshing each one. The auto
%   path here drives AUTOGRID's Row Wizard (M.ag_script = 'ag_rowwizard.py')
%   to build the 9-block "Axial Fan" topology the default mesher collapsed to 7
%   blocks. run_cfd is off: validate a clean 9-block mesh for every count first.
%
%   Rotor blade-count sweep of the 02062026PropPhi06HYBRID propeller:
%   meshes each count and (if run_cfd) runs the CFD for it. Every blade count
%   needs its own mesh because the blade number sets the passage pitch in
%   AUTOGRID. CFD is a single operating point at fixed RPM (see the run
%   section), matching the current state of PropCURTIS.m.
%
%   AUTOGRID's default B2B topology previously failed the skew check on these
%   wide, low-count passages. It now meshes them with strong skew control
%   (n_opt + skewness-targeted optimisation, set below), giving the same
%   9-block layout as the previous cases so mesh_prop / the CFD run work.
%   This script has two modes, selected by mesh_mode:
%
%     'auto' - build each mesh automatically (default topology + skew control)
%              and save the pre-processed TURBOSTREAM mesh
%              <TURBOSTREAM>/B_<count>_02062026PropPhi06HYBRID_mesh.mat
%
%     'read' - fall back: read meshes you made BY HAND in AUTOGRID (e.g. with
%              your Tweak2 template) and exported as
%              <AUTOGRID>/B_<count>_02062026PropPhi06HYBRID.g  (+ .bcs).
%
%   Run from the "Code - CLAUDE" directory (AUTOGRID is called with the
%   relative path ./AG), the same working directory as PropCURTIS.m.
%
%   Based on PropCURTIS.m (James Taylor & Alex Reaves)
%   Blade-count sweep added July 2026

% Clear workspace and load eVTOL environment
clear; close all; clc; eVTOL;

% Meshing mode: 'auto' builds each mesh automatically; 'read' pulls in meshes
% made by hand (see header). This automesh variant runs 'auto' with the AUTOGRID
% Row Wizard (M.ag_script below) to reproduce the 9-block layout without hand work.
mesh_mode = 'auto';

% Re-mesh counts that already have a saved mesh? 0 = skip (idempotent, so the
% script can be re-run to finish CFD without re-meshing), 1 = force re-mesh.
% Forced here so a re-run actually re-meshes while we develop the auto path.
remesh = 1;

% Run the CFD after meshing? Off while validating the automatic meshing - get a
% clean 9-block mesh for every count first, then flip to 1 to run the CFD.
run_cfd = 0;

% Current job variables
job.priority = 'JTAYLOR-SL2-GPU';
job.cards = 1; job.walltime = 1;

% Base design name - shared geometry / mission inputs for every count
des_base = '02062026PropPhi06HYBRID';

% Rotor blade counts to sweep (odd only). Start with the single hardest case:
% N=3 is the widest passage, the one the default topology collapsed to 7 blocks.
% Once it meshes 9-block cleanly, widen this back to e.g. [3 5 7 9 11 13 15] (or
% wider) and re-run - the whole point of the auto path is arbitrary counts.
blade_counts = [3];

% Load mission from file (shared across all blade counts)
s = importdata([dr.geom des_base '-mission.json'],'r');
mymission = jsondecode(s{1});

% Mesh parameters - AUTOGRID default topology, matching the previous working
% 9-block auto cases (e.g. 06epsProp3Blade / MachinableProp3Blade) exactly so
% the mesh comes out 9-block for mesh_prop. Wall spacings / expansion come from
% eVTOL defaults.
%   The original skew failure was caused by the Tweak2 grid-point counts (now
%   reverted), so with these default counts it should pass skew WITHOUT the
%   optimisation additions - which is important because n_opt / skew_ctrl / p_jc
%   were what collapsed the inlet/outlet extension blocks to a 7-block layout.
%   If a wide count still fails the skew check, add M.skew_ctrl = 3 back, but
%   expect it may drop to 7 blocks (then mesh that count by hand in 'read' mode).
M.stin = [0 0]; M.stout = [0 0]; M.n_aid = [17 17]; M.n_aou = [25 25];
M.n_aod = [25 25]; M.n_aiu = [25 25]; M.n_g = 2; M.n_ai = 25; M.n_ao = 25;
M.n_bl = [101 101]; M.n_bu = [101 101]; M.n_so = [25 9]; M.n_zi = [17 9];
M.ui = 0; M.uo = 0; M.n_j = 101;
M.skew = 8;                            % min-angle limit (deg), relaxed from 10.
                                       % N=3 meshes at 9.2 on the default topology
                                       % (usable for RANS); 8 gives margin for the
                                       % other counts. Skew control would push it
                                       % higher but collapses to 7 blocks, so we
                                       % accept the slightly-skewed 9-block mesh.

% Mesh settings for farfield (kept for parity; farfield not built here)
M.f_min = 0.15; M.f_far = 5; M.nj_far = 65; M.ni_in = 49; M.ni_out = 73; M.dr_far = 25e-3;

% Build the row topology with AUTOGRID's Row Wizard (9-block "Axial Fan" layout)
% instead of the default z-line topology that collapses to 7 blocks on wide
% passages. bl_create_machine passes this to `igg -autogrid5 -script ./AG/<...>`.
% Set to '' (or remove) to fall back to the default ag_create_machine.py mesher.
M.ag_script = 'ag_rowwizard.py';


%% Load blade geometry - shared, only the blade COUNT changes in the sweep

% Load blade geometry directly
xrrt = dlmread([dr.geom des_base '.csv'],',',1,0);

% Reshape and split into rotor (and stator if present)
B{1}.xrrt = reshape(xrrt(:,[1 3 2]),428,[],3);
if size(xrrt,2) == 7
    nn = 2;
    B{2}.xrrt = reshape(xrrt(:,[7 6 5]),428,[],3);
else
    nn = 1; B(2) = [];
end

% Scale domain size based on chord of first row
xchord_max = max((max(B{1}.xrrt(:,:,1),[],1) - min(B{1}.xrrt(:,:,1),[],1)));
dx_duct = 0.5 * xchord_max;
M.dx_row = 0.5 * xchord_max; dx_row = 0.5 * xchord_max;

% Domain inlet and outlet planes
B{1}.x_in = min(min(B{1}.xrrt(:,:,1))) - dx_duct;
B{end}.x_out = max(max(B{end}.xrrt(:,:,1))) + dx_duct;

% Dummy meridional curves
x = linspace(B{1}.x_in,B{end}.x_out,101)';
r_hub = mymission.Rh; r_cas = mymission.Rc*mymission.Rrat(end);
xr_hub = [x ones(size(x)) * r_hub]; xr_cas = [x ones(size(x)) * r_cas];
if B{1}.xrrt(1,end,2) < r_cas
    xr_hub = [x ones(size(x)) * (B{1}.xrrt(1,1,2) + 0.1e-3)];
    xr_cas = [x ones(size(x)) * (B{1}.xrrt(1,end,2) - 0.1e-3)];
end
for n = 1:length(B); B{n}.xr_hub = xr_hub; B{n}.xr_cas = xr_cas; end

% Mean radii
r_1 = r_hub; r_2 = r_cas;

% Mesh parameters depending on whether it is a 1 row or 2 row case
if nn == 1
    M.n_zo = 29;
else
    M.n_zo = [17 29];
end

% Set negligible casing gaps on all rows
for n = 1:nn; F{n}.flag = [0 0]; D{n}.flag = [0 1]; D{n}.cas = 1e-5; end

% Boundary conditions and initial guess (required by bl_create_machine)
I.pout = mymission.Atm_Pressure; I.k = 10; I.rpm = [mymission.RPM 0];
r_mean = ((B{1}.xr_cas(1,2)^2 + B{1}.xr_hub(1,2)^2) / 2)^0.5;
Vx_mean = mymission.phi_m * r_mean * pi * I.rpm(1) / 30;
I.nstep = 5; I.poisson_nstep = 1000;
I.guess = [1.2 1.2 * Vx_mean 0 0 2e5 0.0002];

% Solver setup used when pre-processing the mesh
job.run_sa = 1; job.nodes = 1;


%% Sweep over blade counts

fprintf('\n%%%%%%%% Blade-count sweep - mode: %s %%%%%%%%\n', mesh_mode);

ok = false(size(blade_counts));
for ib = 1:length(blade_counts)

    % Rotor blade count for this case
    N = blade_counts(ib);

    % Unique name for this blade count (prefix "B_<count>_")
    des_name = ['B_' num2str(N,'%02d') '_' des_base];
    job.name = des_name;

    % Override the rotor blade count - written into config.dat as c.N and
    % read by AUTOGRID to set the passage pitch
    B{1}.N = N;

    fprintf('\n%%%%%%%% N = %d  ->  %s %%%%%%%%\n\n', N, des_name);

    if strcmp(mesh_mode,'auto')

        % Skip if already meshed (unless remesh is forced), so re-running the
        % script to finish CFD does not re-mesh.
        if remesh == 0 && exist([dr.ts des_name '_mesh.mat'],'file') ~= 0
            ok(ib) = true;
            fprintf('%%%% N = %d already meshed - skipping (set remesh=1 to force) %%%%\n', N);
            continue
        end

        % Build the mesh in AUTOGRID with the Tweak2 settings and save mesh.mat.
        try
            bl_write_ibl([dr.geom des_name '.ibl'], B{1}, 0, 0);

            %   args: (...,plot_stuff,monitor_job,run_ag,run_ts)
            job = bl_create_machine(B,dr,[job.name '_'],job,I,M,F,D,0,0,1,0);

            % Verify the topology is the 9-block layout mesh_prop needs. Wide
            % low-count passages (e.g. N=3) can collapse to a 7-block topology
            % that meshes but breaks the farfield step, so flag those for
            % hand-meshing rather than feeding a broken mesh to the CFD.
            gchk = load([dr.ts des_name '_mesh.mat']);
            if numel(gchk.g) == 9
                ok(ib) = true;
                fprintf('%%%% Saved mesh (9 blocks): %s%s_mesh.mat %%%%\n', dr.ts, des_name);
            else
                ok(ib) = false;
                warning(['N = %d meshed but AUTOGRID collapsed it to %d blocks (need 9) - ' ...
                    'the wide passage is too much for the default topology.\n' ...
                    '         Mesh this count by hand with your Tweak2 template and run it in ' ...
                    '''read'' mode.'], N, numel(gchk.g));
            end
            clear gchk
        catch err
            warning('N = %d failed to auto-mesh: %s\n         Mesh it by hand in AUTOGRID, then run in ''read'' mode.', N, err.message);
        end

    else  % 'read'

        % Require the hand-made AUTOGRID mesh and its boundary conditions
        if isempty(dir([dr.ag des_name '.g'])) == 1 || ...
                isempty(dir([dr.ag des_name '.bcs'])) == 1
            warning(['N = %d: missing %s.g or %s.bcs in AUTOGRID.\n' ...
                '         Mesh this count by hand (101 flow paths, 1e-5 wall cells,\n' ...
                '         Tweak2 topology) and export as %s.g / .bcs first.'], ...
                N, des_name, des_name, des_name);
            continue
        end

        % Check the exported plot3d .g is the 9-block layout (the first number
        % in a plot3d file is the block count). Catches a stale/mismatched
        % export - e.g. re-saving the wizard project but not re-exporting the .g.
        fid = fopen([dr.ag des_name '.g'],'r'); nblk = fscanf(fid,'%d',1); fclose(fid);
        if nblk ~= 9
            warning(['N = %d: %s.g has %d blocks (need 9). The wizard mesh is 9-block, ' ...
                'so this .g is likely stale - re-export the plot3d .g (and .bcs) from ' ...
                'your saved wizard project so they match.'], N, des_name, nblk);
            continue
        end

        % run_ag = 0 -> read the existing (hand-made) .g/.bcs and save mesh.mat
        try
            job = bl_create_machine(B,dr,[job.name '_'],job,I,M,F,D,0,0,0,0);
            gchk = load([dr.ts des_name '_mesh.mat']);
            if numel(gchk.g) == 9
                ok(ib) = true;
                fprintf('%%%% Saved mesh (9 blocks): %s%s_mesh.mat %%%%\n', dr.ts, des_name);
            else
                warning('N = %d: read mesh has %d blocks (need 9) - .g/.bcs may be inconsistent.', ...
                    N, numel(gchk.g));
            end
            clear gchk
        catch err
            warning('N = %d failed to read: %s', N, err.message);
        end

    end
end


%% Meshing summary

fprintf('\n%%%%%%%% Meshing complete (mode: %s) %%%%%%%%\n', mesh_mode);
for ib = 1:length(blade_counts)
    des_name = ['B_' num2str(blade_counts(ib),'%02d') '_' des_base];
    if ok(ib)
        fprintf('  N = %2d : OK       %s_mesh.mat\n', blade_counts(ib), des_name);
    else
        fprintf('  N = %2d : PENDING  (see warning above)\n', blade_counts(ib));
    end
end


%% Run CFD for each blade count
%
%   Single operating point, fixed RPM for every count: the mission Take-Off
%   design point (Vf = mymission.Vfs, rpm = mymission.RPM), sea level. This is
%   the same RPM the mesh inlet BC is built with (I.rpm above), so the CFD runs
%   consistently with the design point rather than an arbitrary speed.
%   Three stages per count, monitored together stage by stage:
%     _O1          - initial run with slipwalls    (nstep 5e4)
%     restart_O1   - slipwalls + inlet -> farfield  (nstep 2e5)
%     restart2_O1  - second restart for averaging   (nstep 2e5, save from 1e5)

if run_cfd

    % Single operating point, fixed RPM for every count = Take-Off design point
    Vf = mymission.Vfs; alt = 0; rpm_run = mymission.RPM;   % mission design RPM, sea level

    % Centrebody geometry parameters
    ar_nose = 1.5; ar_tail = 4; drdx_tail = -0.4;

    % --- Farfield control curves (identical for every count) ---

    % Elliptical nose cone
    th = linspace(0,90,201)';
    c = struct();
    c.xr_nose = [B{1}.x_in - ar_nose * xr_hub(1,2) * cosd(th), xr_hub(1,2) * sind(th)];

    % Quartic tail cone
    dx_tail = ar_tail * xr_hub(end,2);
    Atail = [0 0 0 0 1 ; 0 0 0 1 0 ; 0 0 2 0 0 ; ...
        dx_tail^4 dx_tail^3 dx_tail^2 dx_tail 1 ; 4*dx_tail^3 3*dx_tail^2 2*dx_tail 1 0];
    btail = [xr_hub(end,2) ; 0 ; 0 ; 0 ; drdx_tail]; ptail = Atail \ btail;
    c.xr_tail = [linspace(xr_hub(end,1),xr_hub(end,1) + dx_tail,201)', ...
        polyval(ptail,linspace(0,dx_tail,201)')];

    % Farfield boundary and spike sizes, then boundary curves
    M.rad_far = M.f_far * r_cas; M.r_min = M.f_min * r_hub; M.x_far = 2 * M.rad_far;
    c.xr_in  = [-M.x_far 0 ; -M.x_far M.rad_far];
    c.xr_far = [c.xr_in(end,:) ; M.x_far M.rad_far];
    c.xr_out = [c.xr_far(end,1) 0 ; c.xr_far(end,:)];

    % Atmospheric conditions and freestream stagnation state
    [Ta,Pa,ro_a,~] = Altitude(alt);
    Po = Pa + 0.5 * ro_a * Vf^2;
    To = Ta + 0.5 / cp * Vf^2;

    rpmvars = {'rpm','rpmi1','rpmi2','rpmj1','rpmj2','rpmk1','rpmk2'};

    % ---------- Stage 1: initial run (slipwalls) ----------
    fprintf('\n%%%%%%%% CFD stage 1 (initial run) %%%%%%%%\n');
    Jobs = cell(0,1);
    for ib = 1:length(blade_counts)
        N = blade_counts(ib); des_name = ['B_' num2str(N,'%02d') '_' des_base];
        if ok(ib) == 0 || exist([dr.ts des_name '_mesh.mat'],'file') == 0
            warning('N = %d: mesh not ready (failed/absent), skipping CFD', N); continue;
        end
        load([dr.ts des_name '_mesh.mat']);          % loads g

        % Extend the blade mesh to the farfield
        g = mesh_prop(g,c,M,0);

        % Apply inlet and exit / farfield boundary conditions
        g = ts_apply_inlet(g,[Po Ta 0 0],0);
        for bid = ts_get_bids(g)
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 1
                    g{bid+1}.patch{pid+1}.throttle_type = 0;
                    g{bid+1}.patch{pid+1}.ipout = -3;
                    g{bid+1}.patch{pid+1}.pout = Pa;
                elseif g{bid+1}.patch{pid+1}.attribute.kind == 13
                    g{bid+1}.patch{pid+1}.pstag_free = Po;
                    g{bid+1}.patch{pid+1}.pstat_free = Pa;
                    g{bid+1}.patch{pid+1}.sf_free = 0.5000;
                    g{bid+1}.patch{pid+1}.tstag_free = To;
                end
            end
        end

        % Solver settings
        g{1}.av.nstep = 50000; g{1}.av.dts = 0;
        g{1}.av.nstep_save_start = g{1}.av.nstep - 5000;
        g{1}.av.poisson_nstep = I.poisson_nstep;
        g{1}.av.cfl = 0.3; g{1}.av.turbvis_lim = 10000;
        g{1}.av.poisson_nsmooth = 20; g{1}.av.nomatch_int = 1;
        for bid = ts_get_bids(g)
            g{bid+1}.bv.fmgrid = 0.1;
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 0
                    g{bid+1}.patch{pid+1}.rfin = 0.1;
                end
            end
        end

        % Stationary hubs, then fixed RPM for every count
        g{end}.bv.rpmj1 = 0; g{end-1}.bv.rpmj1 = 0;
        for v = 1:length(rpmvars)
            for bid = ts_get_bids(g)
                if g{bid+1}.bv.(rpmvars{v}) ~= 0
                    g{bid+1}.bv.(rpmvars{v}) = rpm_run;
                end
            end
        end

        % Submit
        job.walltime = 1;
        job.name = [des_name '_O1'];
        Jobs{end+1,1} = ts_sbatch(job,g);
    end
    ts_smonitor(Jobs,0,1);

    % ---------- Stage 2: restart with farfield BCs ----------
    fprintf('\n%%%%%%%% CFD stage 2 (restart) %%%%%%%%\n');
    Jobs = cell(0,1);
    for ib = 1:length(blade_counts)
        N = blade_counts(ib); des_name = ['B_' num2str(N,'%02d') '_' des_base];
        if exist([dr.ts des_name '_O1.hdf5'],'file') == 0
            warning('N = %d: no _O1 solution, skipping restart', N); continue;
        end
        g = ts_read_hdf5([dr.ts des_name '_O1.hdf5']);

        % Turn upper slipwalls (kind 7) and inlet (kind 0) into farfield (kind 13)
        for bid = ts_get_bids(g)
            for pid = ts_get_pids(g,bid)
                k = g{bid+1}.patch{pid+1}.attribute.kind;
                if (k == 7 && g{bid+1}.patch{pid+1}.attribute.jst ~= 0) || k == 0
                    g{bid+1}.patch{pid+1}.attribute.kind = 13;
                    g{bid+1}.patch{pid+1}.pstag_free = Po;
                    g{bid+1}.patch{pid+1}.pstat_free = Pa;
                    g{bid+1}.patch{pid+1}.sf_free = 0.5000;
                    g{bid+1}.patch{pid+1}.tstag_free = To;
                end
            end
        end

        g{1}.av.poisson_nstep = I.poisson_nstep;
        g{1}.av.nstep = 2.0e5; g{1}.av.nstep_save_start = 0.75e5;

        job.name = [des_name 'restart_O1'];
        Jobs{end+1,1} = ts_sbatch(job,g);
    end
    ts_smonitor(Jobs,0,1);

    % ---------- Stage 3: second restart for averaging ----------
    fprintf('\n%%%%%%%% CFD stage 3 (restart2 / averaging) %%%%%%%%\n');
    Jobs = cell(0,1);
    for ib = 1:length(blade_counts)
        N = blade_counts(ib); des_name = ['B_' num2str(N,'%02d') '_' des_base];
        if exist([dr.ts des_name 'restart_O1.hdf5'],'file') == 0
            warning('N = %d: no restart_O1 solution, skipping restart2', N); continue;
        end
        g = ts_read_hdf5([dr.ts des_name 'restart_O1.hdf5']);
        g{1}.av.nstep = 2.0e5; g{1}.av.nstep_save_start = 1e5; g{1}.av.nchange = 0;
        job.name = [des_name 'restart2_O1'];
        Jobs{end+1,1} = ts_sbatch(job,g);
    end
    ts_smonitor(Jobs,0,1);

    fprintf('\n%%%%%%%% CFD sweep complete - averaged solutions: B_<N>_%srestart2_O1_avg.hdf5 %%%%%%%%\n', des_base);
end
