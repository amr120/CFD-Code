%   DuctedFanCURTIS
%                               
%   Generate meshes and CFD solutions for ALL designs
%
%   Written by James Taylor & Alex Reaves             
%   April 2023

% Clear workspace and load eVTOL environment
clear; close all; clc; eVTOL;

% Current job variables
job.priority = 'JTAYLOR-SL2-GPU';
job.cards = 1; job.walltime = 4;

% Script operation
run_case = 1;
restart_case = 0;

% Input design name
des_name = '3OCT2024'; 

% Load mission from file
s = importdata([dr.geom des_name '-mission.json'],'r');
 mymission = jsondecode(s{1});

 ScaleFactor = 3 %mymission.ScaleFactor; %scales velocites so that W1tip =100 m/s

% Flight operating points for mesh_far
% Vf = [1 10 20 30]; alt = zeros(size(Vf)); leg = {'1' '10' '20' '30'};
Vf = mymission.Vfs*ScaleFactor; alt = 0; leg = {'20'};

% Mesh parameters for all rows - high staggered topology
% M.stin = [1 0]; M.stout = [0 0]; M.n_aid = [33 25]; M.n_aou = [33 25];
% M.n_aod = [17 25]; M.n_aiu = [17 25]; M.n_g = 2;
% M.n_bl = [101 101]; M.n_bu = [101 101]; M.n_so = [17 9]; M.n_zi = [17 9]; 
% M.ui = 0; M.uo = 0; M.skew = 5; M.n_j = 101;

% Mesh parameters for all rows - normal O4H topology
M.stin = [0 0]; M.stout = [0 0]; M.n_aid = [17 17]; M.n_aou = [17 17];
M.n_aod = [17 17]; M.n_aiu = [17 17]; M.n_g = 2; M.n_ai = 17; M.n_ao = 17;
M.n_bl = [101 101]; M.n_bu = [101 101]; M.n_so = [17 9]; M.n_zi = [17 9]; 
M.ui = 0; M.uo = 0; M.skew = 10; M.n_j = 101;

% Mesh settings for farfield
M.f_min = 0.15; M.f_far = 5; M.nj_far = 65; M.ni_in = 33; M.ni_out = 53; M.dr_far = 25e-3;

% Centrebody geometry parameters
ar_nose = 1.5; ar_tail = 4; drdx_tail = -0.4;


%% Create mesh for the blades in AUTOGRID, extend to the farfield in MATLAB and run in TURBOSTREAM
if run_case == 1

    % Load blade geometry directly
    xrrt = dlmread([dr.geom des_name '.csv'],',',1,0);

    % Reshape and split into rotor and stator
    B{1}.xrrt = reshape(xrrt(:,[1 3 2]),428,[],3);
    if size(xrrt,2) == 7
        nn = 2;
        B{2}.xrrt = reshape(xrrt(:,[7 6 5]),428,[],3); %was 399 427
    else
        nn = 1; B(2) = [];
    end

    % Plot polar coordinates
    figure(); hold on; grid on; box on; axis equal; cols = lines(length(B));
    for n = 1:length(B)
        plot3(B{n}.xrrt(:,:,1),B{n}.xrrt(:,:,2),B{n}.xrrt(:,:,3),'-','color',cols(n,:))
    end

    % Scale domain size based on chord of first row
    xchord_max = (max(B{1}.xrrt(:,1,1)) - min(B{1}.xrrt(:,1,1)));
    dx_duct = 1 * xchord_max; M.dx_row = 0.2 * xchord_max; %dx_duct = 2 * xchord_max; M.dx_row = 0.2 * xchord_max;

    % Dummy domain size
    B{1}.x_in = min(min(B{1}.xrrt(:,:,1))) - dx_duct; 
    B{end}.x_out = max(max(B{end}.xrrt(:,:,1))) + dx_duct;

    % Dummy merdional curves
    x = linspace(B{1}.x_in,B{end}.x_out,101)'; 
    r_hub = mymission.Rh; r_cas = mymission.Rc;
    %r_hub = 
    xr_hub = [x ones(size(x)) * r_hub]; xr_cas = [x ones(size(x)) * r_cas];    
    for n = 1:length(B); B{n}.xr_hub = xr_hub; B{n}.xr_cas = xr_cas; end;
    
    % Mean radii
    r_1 = r_hub; r_2 = r_cas;

    % Blade counts
    B{1}.N = mymission.NB_r; if nn == 2; B{2}.N = mymission.NB_s; end;
    
    % Mesh parameters depending on whether it is a 1 row or 2 row case
    if nn == 1
        M.n_zo = 29;
    else
        M.n_zo = [17 29];
    end

    % Limit the blade counts and correct the pitch later
    %if B{1}.N < 15; N_temp = B{1}.N; B{1}.N = 15; end;
    
    % Set negligible casing gaps on both rows
    for n = 1:nn; F{n}.flag = [0 0]; D{n}.flag = [0 1]; D{n}.cas = 5e-5; end;

    % Boundary conditions
    I.pout = mymission.Atm_Pressure; I.k = 10; I.rpm = [mymission.RPM*ScaleFactor 0];
    r_mean = ((B{1}.xr_cas(1,2)^2 + B{1}.xr_hub(1,2)^2) / 2)^0.5;
    Vx_mean = mymission.phi_m * r_mean * pi * I.rpm(1) / 30; 

    % Set number of steps and initial guess
    I.nstep = 5000; I.poisson_nstep = 10000; 
    I.guess = [1.2 1.2 * Vx_mean 0 0 2e5 0.0002];

    % Mesh and run the design point in TURBOSTREAM
    job.name = des_name; job.run_sa = 1;
    job = bl_create_machine(B,dr,[job.name '_'],job,I,M,F,D,1,1,1,0);
    if isempty(dir([dr.ts des_name '_avg.*'])) == 1
        job = bl_create_machine(B,dr,[job.name '_'],job,I,M,F,D,1,1,1,1);
    end
    
    % Load in blade rows mesh
    load([dr.ts des_name '_mesh.mat']);

    % Figure window for geometry comparisons
    h = figure(); hold on; grid on; box on; axis equal; cols = lines(7); 
    
    % Plot hub and casing lines
    plot(xr_hub(:,1),xr_hub(:,2),'-','color',cols(3,:));
    plot(xr_cas(:,1),xr_cas(:,2),'-','color',cols(3,:));
        
    % Plot new blade designs
    for n = 1:length(B)
        plot(B{n}.xrrt([1 301],:,1)',B{n}.xrrt([1 301],:,2)','color',cols(1,:))
    end
    
    % Create elliptical nose cone
    th = linspace(0,90,201)';
    c.xr_nose = [B{1}.x_in - ar_nose * xr_hub(1,2) * cosd(th) xr_hub(1,2) * sind(th)];
    
    % Create quartic tail cone
    dx_tail = ar_tail * xr_hub(end,2);
    A = [0 0 0 0 1 ; 0 0 0 1 0 ; 0 0 2 0 0 ; dx_tail^4 dx_tail^3 dx_tail^2 dx_tail^1 1 ; ...
        4*dx_tail^3 3*dx_tail^2 2*dx_tail 1 0];
    b = [xr_hub(end,2) ; 0 ; 0 ; 0 ; drdx_tail]; p = A \ b;
    c.xr_tail = [linspace(xr_hub(end,1),xr_hub(end,1) + dx_tail,201)' ...
        polyval(p,linspace(0,dx_tail,201)')];
    
    % Calculate boundary and spike sizes
    M.rad_far = M.f_far * r_cas; M.r_min = M.f_min * r_hub;
    
    % Create far field boundary curves
    c.xr_in = [-M.rad_far 0 ; -M.rad_far M.rad_far];
    c.xr_far = [c.xr_in(end,:) ; M.rad_far M.rad_far];
    c.xr_out = [c.xr_far(end,1) 0 ; c.xr_far(end,:)];
    
    % Plot all external mesh control curves
    varnames = fieldnames(c);
    for v = 1:length(varnames)
        plot(c.(varnames{v})(:,1),c.(varnames{v})(:,2),'-','color',cols(2,:));
    end 
       
    % Create the mesh
    g = mesh_prop(g,c,M,1);

    % Loop over all operating points and run cases
    Jobs = cell(0,1);
    for n = 1:length(Vf)
    
        % Air conditions at altitude
        [Ta,Pa,ro,~] = Altitude(alt(n));
        
        % Calculate stagnation pressure from flight speed
        %Po = Pa + 0.5 * ro * Vf(n)^2;
        Po = Pa + 0.5 * ro * 100^2;
        

        % Apply boundary conditions for different flight speeds
        g = ts_apply_inlet(g,[Po Ta 0 0],0);
        for bid = ts_get_bids(g)
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 1
                    g{bid+1}.patch{pid+1}.throttle_type = 0;
                    %g{bid+1}.patch{pid+1}.ipout = -3;
                    g{bid+1}.patch{pid+1}.ipout = 0;
                    g{bid+1}.patch{pid+1}.pout = Pa;
                    %g{bid+1}.patch{pid+1}.pout = 100000;
                elseif g{bid+1}.patch{pid+1}.attribute.kind == 13
                    g{bid+1}.patch{pid+1}.pstag_free = Po;
                    g{bid+1}.patch{pid+1}.pstat_free = Pa;
                    g{bid+1}.patch{pid+1}.sf_free = 0.5000;
                    g{bid+1}.patch{pid+1}.tstag_free = Ta;
                end
            end
        end
        
        % Set number of steps

        g{1}.av.nstep = 100000; g{1}.av.dts = 0;
        %g{1}.av.nstep = 1; g{1}.av.dts = 0;
        %g{1}.av.nstep_save_start = g{1}.av.nstep-000;
        g{1}.av.nstep_save_start = g{1}.av.nstep-1000;
        g{1}.av.poisson_nstep = I.poisson_nstep;
        

        for bid = ts_get_bids(g)
            g{bid+1}.av.cfl = 0.4;
            g{bid+1}.av.cfl_mul = 0.5;
            g{bid+1}.av.turbvis_lim = 10000; %10000
        end

        % Debugging settings
%         g{1}.av.nstep = 18600;
%         g{1}.av.ilos = 1; des_name = [des_name '_ML'];
%         des_name = [des_name '_XL'];
        
        % Submit the job
        %job.name = [des_name '_PreNan_9750' '_O' num2str(n)];
        job.name = [des_name '_O' num2str(n)];
        Jobs{end+1,1} = ts_sbatch(job,g);

    end
    
    % Monitor all jobs until completed - allow reversed flow in transients
    ts_smonitor(Jobs,0,1);    

end

%% Restart all cases with upstream and downstream blocks set to flight velocity
if restart_case == 1
    
    % Loop over all operating points and run cases
    Jobs = cell(0,1);
    for n = 1:length(Vf)
        
        % Read previous solution
        g = ts_read_hdf5([dr.ts des_name '_O' num2str(n) '_avg.hdf5']);
        
        % Reset far field flow
        ro = mean(g{end-1}.ro(:)); roe = mean(g{end-1}.roe(:));
        for bid = length(g)-2:length(g)-1
            g{bid+1}.ro(:) = ro;
            g{bid+1}.rovx(:) = ro * Vf(n);
            g{bid+1}.rovr(:) = 0;      
            g{bid+1}.rorvt(:) = 0;  
            g{bid+1}.roe(:) = roe;
        end
        
        % Increase number of steps
        g{1}.av.nstep = 4e5; g{1}.av.nstep_save_start = 3e5; 
        
        % Submit the job
        job.name = [des_name '_O' num2str(n) '_restart'];
        Jobs{end+1,1} = ts_sbatch(job,g);
        
    end
    
    % Monitor all jobs until completed - allow reversed flow in transients
    ts_smonitor(Jobs,0,1); 
    
end

