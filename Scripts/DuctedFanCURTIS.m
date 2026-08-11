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
job.cards = 1;
job.walltime = 4;

% Script operation
create_cowl = 0;
run_des = 1;
mesh_lift = 0;
mesh_cowl = 0;
restart_cowl = 0;
plot_pressure = 0;
plot_chics = 1;
plot_thrust = 1;
write_pv = 1;
write_ibl = 0;

% Input design name and parameters for run_des

des_name = ['Phi06DUCTEDFANFIXRPM']; 


s = importdata([dr.geom append(des_name,'-mission.json')],'r');
%s= importdata([dr.geom append('mission.json')],'r');
%mission = dlmread([dr.geom 'mission.txt'],',',1,0);

%myOrder = struct(s{1});
mymission = jsondecode(s{1});
phi = mymission.phi_m; sig = mymission.sigma; rpm = mymission.RPM*mymission.ScaleFactor;

% Flight operating points for mesh_cowl
%Vf = [1 10 20 30]; alt = zeros(size(Vf)); leg = {'1' '10' '20' '30'}; T_des = [50 40 30 30];
Vf = mymission.Vfs*mymission.ScaleFactor; leg = mymission.Vfs*mymission.ScaleFactor; alt = zeros(size(Vf)); T_des  = mymission.Thrust;
% Post-processing settings
plot_avg = 1;


%% Create cowl geometry and modify meridional lines of blades 
if create_cowl == 1
    
    % Load blade designs
     load([dr.geom des_name '_B1_spline.mat']); B{1} = b;
     load([dr.geom des_name '_B2_spline.mat']); B{2} = b;
    %load([dr.geom des_name '_prelim.mat']); 
    
    % Move stator downstream
    %B{2}.x_ref = B{2}.x_ref + 0.005;
    
    % First geometry - too thick for Darren
%     r_cas = B{1}.xr_cas(1,2); r_hub = B{1}.xr_hub(1,2);
%     r_lip = 0.0712; t_te = 1e-3; x_lip = -0.01; xchord = 0.17;
%     x_nose = 0.01; e_nose = 2; x_tail = 0.14; l_exit = 0.09;
%     l_cone = 0.086; thick_max = 0.028; drdx_tail = -0.4; rad_in = 0.48;

%     % Second geometry - with 160mm OD
%     r_cas = B{1}.xr_cas(1,2); r_hub = B{1}.xr_hub(1,2);
%     r_lip = 0.0688; t_te = 1e-3; x_lip = 0; xchord = 0.16;
%     x_nose = 0.01; e_nose = 2; x_tail = 0.14; l_exit = 0.09;
%     l_cone = 0.086; thick_max = 0.0237; drdx_tail = -0.4; rad_in = 0.48; 
    
    % Third geometry - with 160mm fan diameter
%     r_cas = B{1}.xr_cas(1,2); r_hub = B{1}.xr_hub(1,2);
%     r_lip = 0.0892; t_te = 1e-3; x_lip = -0.02; xchord = 0.19;
%     x_nose = -0.01; e_nose = 2; x_tail = 0.17; l_exit = 0.09;
%     l_cone = 0.12; thick_max = 0.023; drdx_tail = -0.4; rad_in = 0.48; 
    
    % Fourth geometry - with 210mm diameter cowl
    r_cas = B{1}.xr_cas(1,2); r_hub = B{1}.xr_hub(1,2);
    r_lip = 0.095; t_te = 1e-3; x_lip = -0.02; xchord = 0.2;
    x_nose = 0; e_nose = 1.5; x_tail = 0.17; l_exit = 0.1;
    l_cone = 0.12; thick_max = 0.025; drdx_tail = -0.4; rad_in = 0.48; 
    
    % Cowl camberline
    c_out.chi_le = -2; c_out.qcam = 0; c_out.dcam_le = -2; c_out.dcam_te = -0.5;
    
    % Outer thickness distribution
    c_out.rad_le = 5; c_out.thick_max = thick_max; c_out.s_thick_max = 0.25; c_out.thick_te = 0;
    c_out.wedge_te = 63; c_out.rad_thick_max = 0.12; c_out.tchord = xchord;
    
    % Inner thickness distribution
    c_in = c_out; st = 0.9; c_in.thick_max = st * c_in.thick_max; c_in.rad_le = c_in.rad_le / st;
    c_in.s_thick_max = 0.2; c_in.rad_thick_max = 0.12; c_in.wedge_te = 50;
    
    % Figure window for geometry comparisons
    h = figure('position',[1 41 1920 970]); hold on; grid on; box on; cols = lines(7); 
    axis equal; axis([-0.03 0.18 0 0.11]);
     
    % Plot current blade designs
%     for n = 1:2
%         mesh(B{n}.xrrt(:,:,1),B{n}.xrrt(:,:,2),zeros(601,47),'edgecolor',cols(1,:))
% %         plot3(B{n}.xrrt([1 301],:,1)',B{n}.xrrt([1 301],:,2)',zeros(2,43)','color',cols(1,:))
%     end

    % Calculate exit radius from area ratio
    A_fan = pi * (r_cas^2 - r_hub^2);
    A_exit = sig * A_fan; r_exit = (A_exit / pi)^0.5;
    
    % Bisect on exit angle to find correct trailing edge radius
    chi_min = -45; chi_max = 0; r_end = 1;
    while abs(r_exit + t_te - r_end) > 1e-9
        
        % Generate camberline
        c_out.chi_te = 0.5 * (chi_min + chi_max);
        m = bl_construct_camber(c_out,linspace(0,1,301)');
        
        % Scale and shift camberline
        xr_cam = m.xy_cam * xchord / m.xy_cam(end,1);
        xr_cam(:,1) = xr_cam(:,1) + x_lip;
        xr_cam(:,2) = xr_cam(:,2) + r_lip;
        plot(xr_cam(:,1),xr_cam(:,2))

        % Update limits on exit angle
        r_end = xr_cam(end,2);
        if r_end < r_exit + t_te
            chi_min = c_out.chi_te;
        else
            chi_max = c_out.chi_te;
        end
        
    end
    
    % Generate outer cowl profile
    [xr_out,xr_cam] = bl_construct_section(c_out,1,1); xr_out = xr_out(301:end,:); 
    
    % Generate inner cowl profile
    c_in.chi_te = c_out.chi_te;
    xr_in = bl_construct_section(c_in,1,1); xr_in = flip(xr_in(101:301,:),1); 
    
    % Assemble inner and outer profiles
    xr = [flip(xr_in(2:end,:),1) ; xr_out];
    
    % Dimensionalise thickness distribution
    xr = xr * xchord / xr(end,1);
    xr_cam = xr_cam * xchord / xr_out(end,1); 
    c.xr_cowl = [xr(:,1) + x_lip xr(:,2) + r_lip];
    xr_cam = [xr_cam(:,1) + x_lip xr_cam(:,2) + r_lip];

    % Split out intake curve
    [~,i_lip] = min(c.xr_cowl(:,1));
    c.xr_int = flip(c.xr_cowl(1:i_lip,:),1);
    
    % Create elliptical nose cone
    th = linspace(0,90,201)';
    c.xr_nose = [x_nose + e_nose * r_hub * (1-cosd(th)) r_hub * sind(th)];
    
    % Create conical nozzle
    x_exit = x_lip + xchord;
    c.xr_exit = [linspace(x_exit-l_exit,x_exit,201)' linspace(r_cas,r_exit,201)'];
    
    % Create quartic tail cone
    A = [0 0 0 0 1 ; 0 0 0 1 0 ; 0 0 2 0 0 ; l_cone^4 l_cone^3 l_cone^2 l_cone^1 1 ; ...
        4*l_cone^3 3*l_cone^2 2*l_cone 1 0];
    b = [r_hub ; 0 ; 0 ; 0 ; drdx_tail]; p = A \ b;
    c.xr_tail = [linspace(x_tail-l_cone,x_tail,201)' polyval(p,linspace(0,l_cone,201)')];
    
    % Plot all external mesh control curves
%     figure(h); varnames = fieldnames(c);
%     for v = 1:length(varnames)
%         plot(c.(varnames{v})(:,1),c.(varnames{v})(:,2),'-','color',cols(5,:));
%     end

    % Modify meridional lines for rotor to line up with cowl and nose bullet
    r_1 = interp1(c.xr_nose(:,1),c.xr_nose(:,2),B{1}.xrrt(1,2,1));
    r_2 = interp1(c.xr_int(:,1),c.xr_int(:,2),B{1}.xrrt(1,31,1));
    %r_2 = interp1(c.xr_int(:,1),c.xr_int(:,2),B{1}.xrrt(1,42,1));
    B{1}.A_le = pi * (r_2^2 - r_1^2); B{1}.r_le = 0.5 * (r_1 + r_2);
    r_1 = r_hub; r_2 = r_cas;
    B{1}.A_te = pi * (r_2^2 - r_1^2); B{1}.r_te = 0.5 * (r_1 + r_2); 
    drdx_cowl = grad_mg(c.xr_int(:,1),c.xr_int(:,2));
    drdx_nose = grad_mg(c.xr_nose(:,1),c.xr_nose(:,2));
    B{1}.drdx_cas_in = interp1(c.xr_int(:,1),drdx_cowl,B{1}.x_in);
    B{1}.drdx_hub_in = interp1(c.xr_nose(:,1),drdx_nose,B{1}.x_in);
    %B{1}.drdx_cas_le = interp1(c.xr_int(:,1),drdx_cowl,B{1}.xrrt(1,32,1));
    B{1}.drdx_cas_le = interp1(c.xr_int(:,1),drdx_cowl,B{1}.xrrt(1,31,1));
    B{1}.drdx_hub_le = interp1(c.xr_nose(:,1),drdx_nose,B{1}.xrrt(1,2,1));
    B{1}.drdx_cas_te = 0;
    B{1} = bl_construct_blade(B{1});
    
    % Modify meridional lines for stator to line up with nozzle and tail cone
    r_1 = interp1(c.xr_tail(:,1),c.xr_tail(:,2),B{2}.xrrt(1,2,1)); 
    if B{2}.xrrt(1,31,1) > c.xr_exit(1,1)
        r_2 = interp1(c.xr_exit(:,1),c.xr_exit(:,2),B{2}.xrrt(1,31,1));
        B{2}.drdx_cas_le = interp1(c.xr_exit(:,1),drdx_exit,B{2}.xrrt(1,31,1));
    else
        r_2 = r_cas; drdx_cas_le = 0;
    end
    B{2}.A_le = pi * (r_2^2 - r_1^2); B{2}.r_le = 0.5 * (r_1 + r_2);
    r_1 = interp1(c.xr_tail(:,1),c.xr_tail(:,2),B{2}.xrrt(214,2,1));
    r_2 = interp1(c.xr_exit(:,1),c.xr_exit(:,2),B{2}.xrrt(214,31,1));
    B{2}.A_te = pi * (r_2^2 - r_1^2); B{2}.r_te = 0.5 * (r_1 + r_2); 
    drdx_exit = grad_mg(c.xr_exit(:,1),c.xr_exit(:,2));
    drdx_tail = grad_mg(c.xr_tail(:,1),c.xr_tail(:,2));
    B{2}.drdx_cas_out = interp1(c.xr_exit(:,1),drdx_exit,B{2}.x_out);
    B{2}.drdx_hub_out = interp1(c.xr_tail(:,1),drdx_tail,B{2}.x_out);
    B{2}.drdx_hub_le = interp1(c.xr_tail(:,1),drdx_tail,B{2}.xrrt(1,2,1));
    B{2}.drdx_hub_te = interp1(c.xr_tail(:,1),drdx_tail,B{2}.xrrt(214,2,1));
    B{2}.drdx_cas_te = interp1(c.xr_exit(:,1),drdx_exit,B{2}.xrrt(214,31,1));
    B{2} = bl_construct_blade(B{2});
    
    % Create all meridional lines for blade rows
    b = bl_construct_mer(B,5e-4); xr_hub = b.xr_hub; xr_cas = b.xr_cas;
    
    % Plot hub and casing lines
    figure(h);
    plot(xr_hub(:,1),xr_hub(:,2),'-','color',cols(3,:));
    plot(xr_cas(:,1),xr_cas(:,2),'-','color',cols(3,:));
        
    % Plot new blade designs
    for n = 1:2
        plot3(B{n}.xrrt([1 301],:,1)',B{n}.xrrt([1 301],:,2)',zeros(2,43)','color',cols(1,:))
    end
    
    % Plot all external mesh control curves
%     figure(h); varnames = fieldnames(c);
%     for v = 1:length(varnames)
%         plot(c.(varnames{v})(:,1),c.(varnames{v})(:,2),'-','color',cols(5,:));
%     end
    
    % Trim curves and morph to match blade passage coordinates    
    xr = [xr_hub([1 end],:) ; xr_cas([1 end],:)];
    varnames = {'xr_nose' 'xr_tail' 'xr_int' 'xr_exit'};
    for v = 1:length(varnames)
        
        % Find out of range points
        if v == 2 || v == 4
            q = c.(varnames{v})(:,1) < xr(v,1);
        else
            q = c.(varnames{v})(:,1) > xr(v,1);
        end
        
        % Delete out of range points
        c.(varnames{v})(q,:) = [];
        
        % Apply shifts to ensure continuity
        if v == 2 || v == 4
            dxr = xr(v,:) - c.(varnames{v})(1,:);
            c.(varnames{v}) = c.(varnames{v}) + linspace(1,0,size(c.(varnames{v}),1))' * dxr;
        else
            dxr = xr(v,:) - c.(varnames{v})(end,:);
            c.(varnames{v}) = c.(varnames{v}) + linspace(0,1,size(c.(varnames{v}),1))' * dxr;
        end        
        
    end
    
    % Join intake curve back to the cowl
    c.xr_cowl = [flip(c.xr_int,1) ; c.xr_cowl(i_lip+1:end,:)];
    
    % Create far field boundary curve
    th = linspace(0,90,201)';
    c.xr_in = [-rad_in * cosd(th) rad_in * sind(th) ; 1.5 * rad_in rad_in];
    
    % Outlet boundary condtion curve
    c.xr_out = [c.xr_in(end,:) ; c.xr_in(end,1) 0];
    
    % Plot all external mesh control curves
    figure(h);
    plot(xr_cam(:,1),xr_cam(:,2),'color',cols(2,:))
    varnames = fieldnames(c);
    for v = 1:length(varnames)
        plot(c.(varnames{v})(:,1),c.(varnames{v})(:,2),'-','color',cols(2,:));
    end
    
    % Save all curves and blades
    save([dr.geom des_name 'curves.mat'],'c')
    save([dr.geom des_name 'splines.mat'],'B','I')
    
end


%% Create mesh for the design and run open throttle point
if run_des == 1

    % Load blade geometries directly
    xrrt = dlmread([dr.geom des_name '.csv'],',',1,0);

    %mymap = containers.Map(fieldnames(kvp), struct2cell(kvp));
    
   
   
    % Reshape and split into rotor and stator
    B{1}.xrrt = reshape(xrrt(:,[1 3 2]),428,[],3);
    B{2}.xrrt = reshape(xrrt(:,[7 6 5]),428,[],3);

    % Reorder to begin at the leading edge
    %for n = 1:length(B)
    %    i_le = round(size(B{n}.xrrt,1)/2);
    %    B{n}.xrrt = B{n}.xrrt([i_le:end 2:i_le],:,:);
    %end

    %figure; hold on; grid on; box on;
    %plot3(B{1}.xrrt(:,1,1),B{1}.xrrt(:,1,3),size(B{1}.xrrt,1));
    %kjgnegnergnerg

    % Plot polar coordinates
    figure(); hold on; grid on; box on; axis equal; cols = lines(length(B));
    for n = 1:length(B)
        plot3(B{n}.xrrt(:,:,1),B{n}.xrrt(:,:,2),B{n}.xrrt(:,:,3),'-','color',cols(n,:))
    end

    % Scale domain size based on chord of first row
    dx_duct = 0.4 * max((max(B{1}.xrrt(:,:,1),[],1) - min(B{1}.xrrt(:,:,1),[],1))); 
    dx_row = 0.2 * dx_duct;

    %added by alex-----------
    xchord_max = (max(B{1}.xrrt(:,1,1)) - min(B{1}.xrrt(:,1,1)));
    dx_duct = 1 * xchord_max; M.dx_row = 0.2 * xchord_max; 
    %----------------------

    % Dummy domain size
    B{1}.x_in = min(min(B{1}.xrrt(:,:,1))) - dx_duct; 
    B{end}.x_out = max(max(B{end}.xrrt(:,:,1))) + dx_duct;

    % Dummy merdional curves
    %TEST = dlmread([dr.geom des_name 'Mission.csv'],',',1,0);
    x = linspace(B{1}.x_in,B{end}.x_out,101)'; 
    %TEST.

    %true_hub = mymission.Rh; 
    %true_cas = mymission.Rc;
    %xr_hub = [x ones(size(x)) * (true_hub)];
    %xr_cas = [x ones(size(x)) * (true_cas)];

    xr_hub = [x ones(size(x)) * (B{2}.xrrt(1,1,2) + 0.1e-3)];
    xr_cas = [x ones(size(x)) * (B{2}.xrrt(1,end,2) - 0.1e-3)];


    
    for n = 1:length(B); B{n}.xr_hub = xr_hub; B{n}.xr_cas = xr_cas; end;

    %r_1 = mymission.Rh; r_2 = mymission.Rc;
    r_1 = xr_hub(1,2); r_2 = xr_cas(1,2);
    B{1}.A_le = pi * (r_2^2 - r_1^2); B{1}.r_le = 0.5 * (r_1 + r_2);
    B{1}.A_te = pi * (r_2^2 - r_1^2); B{1}.r_te = 0.5 * (r_1 + r_2); 
    %drdx_cowl = grad_mg(c.xr_int(:,1),c.xr_int(:,2));
    B{2}.A_le = pi * (r_2^2 - r_1^2); B{2}.r_le = 0.5 * (r_1 + r_2);
    B{2}.A_te = pi * (r_2^2 - r_1^2); B{2}.r_te = 0.5 * (r_1 + r_2); 

    B{1}.drdx_hub_le = 0; 
    B{2}.drdx_hub_le = 0; 
    B{1}.drdx_hub_te = 0; 
    B{2}.drdx_hub_te = 0; 
    B{1}.drdx_cas_te = 0; 
    B{2}.drdx_cas_te = 0; 
    B{1}.drdx_cas_le = 0; 
    B{2}.drdx_cas_le = 0; 



    % Dummy blade counts
    B{1}.N = mymission.NB_r; B{2}.N = mymission.NB_s;

    % Mesh parameters for all rows
    M.stin = [0 0]; M.stout = [-1 0]; M.n_aid = [61 17]; M.n_aou = [61 17];
    M.n_aod = [17 17]; M.n_aiu = [17 17];
    M.n_bl = [81 81]; M.n_bu = [81 81]; M.n_so = [17 17]; M.n_zi = [17 17]; 
    M.n_zo = [17 17]; M.ui = 0; M.uo = 0; M.skew = 10; M.dx_row = dx_row;

    % Mesh parameters for all rows
    %M.stin = [0 0]; M.stout = [-1 0]; M.n_aid = [33 25]; M.n_aou = [33 25];
    %M.n_aod = [17 25]; M.n_aiu = [17 25];
    %M.n_bl = [53 65]; M.n_bu = [101 65]; M.n_so = [17 9]; M.n_zi = [9 9]; 
    %M.n_zo = [9 25]; M.ui = 0; M.uo = 0; M.skew = 10; M.dx_row = 2e-3;

    %M.zo is outlet


    % Turn off fillets and gaps
    for n = 1:2; F{n}.flag = [0 0]; D{n}.flag = [0 0]; end;

    % Target mass flow
    I.pout = mymission.Atm_Pressure; I.k = 10; mass_fac = 1.04; I.rpm = [mymission.RPM*mymission.ScaleFactor 0];
    r_mean = ((B{1}.xr_cas(1,2)^2 + B{1}.xr_hub(1,2)^2) / 2)^0.5;
    Vx_mean = phi * r_mean * pi * I.rpm(1) / 30; 
    %A_fan = pi * (xr_cas(1,2)^2 - xr_hub(1,2)^2);
    A_fan = mymission.Ax;
    I.m = 1.2 * Vx_mean * A_fan * mass_fac; 
    %I.m = mymission.mdot;

    % Set number of steps and initial guess
    I.nstep = 200000; I.poisson_nstep = 10000; 

     I.guess = [1.2 1.2 * Vx_mean 0 0 2e5 0.0002];
%     I.guess = [dr.ts 'IPM5_guess.dat'];
 %   I.guess = [dr.ts 'Geom2_avg.hdf5'];
    %I.guess = [dr.ts 'EPS-0.3_avg.hdf5'];
    
    % Mesh and run the design point in TURBOSTREAM
    job.name = des_name; job.run_sa = 1;
    job.walltime = 1;
    %if new job name, run autogrid
    if isempty(dir([dr.ag des_name '.*'])) == 1
        job = bl_create_machine(B,dr,[job.name '_'],job,I,M,F,D,1,1,1,0);
    end

    %job.name = [des_name 'changemdot''_'];
    %after manual meshing read in .g from autogrid
     if isempty(dir([dr.ts des_name  '_avg.*'])) == 1
        job = bl_create_machine(B,dr,[job.name  '_'],job,I,M,F,D,1,1,0,1);
    end
    
    % Define cut planes for evaluation and run the chic
    %coords = [2 7 ; 7 19];
    ts_run_chic(job,[],1,[],0,0.005)

    %run like the propeller cases
    %ts_sbatch(job);
    %Jobs{end+1,1} = ts_sbatch(job,g);

    end


%% Create around the cowl curves for liftfan configuration and run in TURBOSTREAM
if mesh_lift == 1
    
    % Mesh settings
    M.n_bl = 25; M.n_int = 45; M.ni_ext = 65; M.nj_ext = 41; M.n_noz = 41; 
    M.n_te = 25; M.n_k = 101; M.r_min = 0.002; M.d_fac = 1.5; M.d_lip = 0.8;
    M.dx = 2; M.dy = 2; M.alpha = 20; M.r_fac = 2; M.i_corn = 9;
    M.mesh_zero = 0;
    
    % Boundary conditions
    M.Vf = 80; M.h = 0;
    
    % Load in existing mesh and curves
    g = ts_read_hdf5([dr.ts des_name '.hdf5']);
    load([dr.geom des_name '_curves.mat']);
    c = rmfield(c,{'xr_in' 'xr_out'});
    
    % Update rpms to improve convergence
    rpm_ids = {'rpm' , 'rpmi1' , 'rpmi2' , 'rpmj1' , 'rpmj2' , 'rpmk1' , 'rpmk2'};
    for bid = ts_get_bids(g)
        for v = 1:length(rpm_ids)
            g{bid+1}.bv.(rpm_ids{v}) = 2 * g{bid+1}.bv.(rpm_ids{v});
        end
    end
        
    % Create the mesh
    g = mesh_liftfan(g,c,M,1);

    % Loop over all operating points and run cases
    Jobs = cell(0,1);

    % Write yplus values
    g{1}.av.write_yplus = 1;

    % Frozen rotor settings
    g{1}.av.nstep = 100000; g{1}.av.dts = 0; g{1}.av.nstep_save_start = 90000;

    % Submit the job
    job.cards = 2; job.walltime = 12;
    job.name = [des_name '_lift_froz'];
    job = ts_sbatch(job,g);

    % Monitor job until completed
    ts_smonitor(job,0,1);    
    
    % Read the job and export for PARAVIEW in cartesians
    g = ts_read_hdf5([dr.ts job.name '_avg.hdf5']);
    ts_export_paraview(g,[dr.pv job.name '.hdf5'],{'Po' 'mwall' 'To' 'P'},[],1);

end


%% Create external mesh around the cowl curves and run in TURBOSTREAM
if mesh_cowl == 1

    %-------------------------ADDED BY ALEX -----------------------------
    % Define overall nacelle parameters
    c.r_lip = 0.085; c.t_te = 1e-3; c.x_lip = -0.09;
    c.thick_max = 0.02; c.chi_le = 5; 

    %r_lip - radius of outer cowl
    %t_te - trailing edge thickness of nacelle
    %X_lip - axial starting point
    %chi_le = pitch angle at leading edge
    
    % Define detailed nacelle section parameters
    c.qcam = 0; c.dcam_le = -2; c.dcam_te = -0.5; c.qcam = 0;
    c.wedge_te = 55; c.rad_le = 5; c.rad_thick_max = 0.12; c.s_thick_max = 0.25; 
    c.rad_le = 10;

    % Define nozzle parameters
    c.l_noz = 0.05; c.sigma = mymission.sigma;

    %l_nozzle - length in meters
    
    % Define centrebody parameters
    c.e_nose = 1.5; c.l_tail = 0.05; c.drdx_tail = -0.4;
    
    % Modify the gradients to test the blends
    B{1}.drdx_hub_in = 0.1; B{1}.drdx_cas_in = 0; 
    B{end}.drdx_hub_out = -0.1; B{end}.drdx_cas_out = 0; 
    
    B{1}.drdx_hub_le = 0; 
    B{2}.drdx_hub_le = 0; 
    B{1}.drdx_hub_te = 0; 
    B{2}.drdx_hub_te = 0; 
    B{1}.drdx_cas_te = 0; 
    B{2}.drdx_cas_te = 0; 
    B{1}.drdx_cas_le = 0; 
    B{2}.drdx_cas_le = 0; 

    %r_1 = interp1(c.xr_nose(:,1),c.xr_nose(:,2),B{1}.xrrt(1,2,1));
    %r_2 = interp1(c.xr_int(:,1),c.xr_int(:,2),B{1}.xrrt(1,31,1));
    %r_2 = interp1(c.xr_int(:,1),c.xr_int(:,2),B{1}.xrrt(1,42,1));
    B{1}.A_le = pi * (r_2^2 - r_1^2); B{1}.r_le = 0.5 * (r_1 + r_2);
    %r_1 = true_hub; r_2 = true_cas;
    B{1}.A_te = pi * (r_2^2 - r_1^2); B{1}.r_te = 0.5 * (r_1 + r_2); 

    B{2}.A_le = pi * (r_2^2 - r_1^2); B{2}.r_le = 0.5 * (r_1 + r_2);
    r_1 = interp1(c.xr_tail(:,1),c.xr_tail(:,2),B{2}.xrrt(214,2,1));
    r_2 = interp1(c.xr_exit(:,1),c.xr_exit(:,2),B{2}.xrrt(214,31,1));
    B{2}.A_te = pi * (r_2^2 - r_1^2); B{2}.r_te = 0.5 * (r_1 + r_2); 

    

    c = bl_construct_nacelle(c,B,1);

    % Save all curves and blades
    save([dr.geom des_name '_curves.mat'],'c')
    save([dr.geom des_name '_splines.mat'],'B','I')

    % Remove unncessary 0D numbers from c
    varnames = fieldnames(c); 
    for v = 1:length(varnames)
        if numel(c.(varnames{v})) > 1
            duct.(varnames{v}) = c.(varnames{v});
        end
    end

    % Join the intake to the outer cowl surface
    duct.xr_cowl = [flip(duct.xr_int(2:end,:),1) ; duct.xr_out];

    % Create far field boundary curve
    th = linspace(0,90,201)'; rad_in = 0.48; 
    duct.xr_in = [-rad_in * cosd(th) rad_in * sind(th) ; 1.5 * rad_in rad_in];
    
    % Outlet boundary condtion curve
    duct.xr_out = [duct.xr_in(end,:) ; duct.xr_in(end,1) 0];

    % Rename the nozzle
    duct.xr_exit = duct.xr_noz;

    %-------------------------------------------
    % Mesh settings
    M.n_bl = 25; M.n_int = 45; M.n_ext = 101; M.n_noz = 41; M.n_out = 73; M.sj_jet = 0.4; M.n_in = 41; 
    M.n_te = 25; M.dt = 60; M.n_k = 25; M.r_min = 0.002; M.d_fac = 1.5; M.d_lip = 0.8; M.d_igv = 0.8;
    
    % Load in existing mesh and curves
    g = ts_read_hdf5([dr.ts des_name '.hdf5']);
    %load([dr.geom des_name '_curves.mat']);
    %load([dr.geom des_name '_spline.mat']);
        
    % Create the mesh
    g = mesh_external(g,duct,M,1);

    % Loop over all operating points and run cases
    Jobs = cell(0,1);
    for n = 1:length(Vf)
        
        % Read the initial guess from a previous successful case
%         g_guess = ts_read_hdf5([dr.ts 'IPM5_V1BE_O' num2str(n) '.hdf5']);
    
        % Air conditions at altitude
        [Ta,Pa,ro,~] = Altitude(alt(n));
        
        % Calculate stagnation pressure from flight speed
        Po = Pa + 0.5 * ro * Vf(n)^2;

        % Apply boundary conditions for different flight speeds
        g = ts_apply_inlet(g,[Po Ta 0 0],0);
        for bid = ts_get_bids(g)
            for pid = ts_get_pids(g,bid)
                if g{bid+1}.patch{pid+1}.attribute.kind == 1
                    g{bid+1}.patch{pid+1}.throttle_type = 0;
                    g{bid+1}.patch{pid+1}.ipout = 0;
%                     if bid == length(g)-1
%                         g{bid+1}.patch{pid+1}.ipout = 0;
%                     else
%                         g{bid+1}.patch{pid+1}.ipout = -3;
%                     end
                    g{bid+1}.patch{pid+1}.pout = Pa;
                elseif g{bid+1}.patch{pid+1}.attribute.kind == 13
                    g{bid+1}.patch{pid+1}.pstag_free = Po;
                    g{bid+1}.patch{pid+1}.pstat_free = Pa;
                    g{bid+1}.patch{pid+1}.sf_free = 0.5000;
                    g{bid+1}.patch{pid+1}.tstag_free = Ta;
                end
            end
        end
        
        % Apply rotational speeds for different operating points
%         rpm_ids = {'rpm' 'rpmi1' 'rpmi2' 'rpmj1' 'rpmj2' 'rpmk1' 'rpmk2'};
%         for bid = ts_get_bids(g)
%             for v = 1:length(rpm_ids)
%                 if g{bid+1}.bv.(rpm_ids{v}) ~= 0
%                     g{bid+1}.bv.(rpm_ids{v}) = omega(n) * 30 / pi;
%                 end
%             end
%         end

        % Set initial guess in far field from previous solution
%         varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe' 'trans_dyn_vis'};
%         for bid = length(g)-5:length(g)-1
%             for v = 1:length(varnames)
%                 g{bid+1}.(varnames{v})(:) = median(g_guess{bid+1}.(varnames{v})(:));
%             end
%         end
        
        % Reduce CFL number
        g{1}.av.cfl = 0.3;
        
        % Reduce problems at mixing planes
        g{1}.av.rfmix = 0.01;

        % Write yplus values
        g{1}.av.write_yplus = 1;
        
        % Increase number of steps
        g{1}.av.nstep = 100000; g{1}.av.dts = 0;
        g{1}.av.nstep_save_start = 90000;
        g{1}.av.nstep = 2440;
        
        % Submit the job
%         job.name = [des_name 'E_guess_O' num2str(n)];
        job.name = [des_name 'E_guess_O' num2str(n) '_froz'];
        Jobs{end+1,1} = ts_sbatch(job,g);

    end
    
    % Monitor all jobs until completed - allow reversed flow in transients
    ts_smonitor(Jobs,0,1);    

end


%% Restart all cowl cases with freestream block reset to flight velocity
if restart_cowl == 1
    
    % Number of hours to run for
    job.walltime = 20;
    
    % Loop over all operating points and run cases
    Jobs = cell(0,1);
    for n = 1:length(Vf)
        
        % Read previous solution
        g = ts_read_hdf5([dr.ts des_name 'E_guess_O' num2str(n) '_avg.hdf5']);
        
        % Reset far field flow
        ro = mean(g{end}.ro(:));
        g{end}.ro(:) = ro;
        g{end}.rovx(:) = ro * Vf(n);
        g{end}.rovr(:) = 0;      
        g{end}.rorvt(:) = 0;  
        g{end}.roe(:) = mean(g{end}.roe(:));

        % Increase number of steps
        g{1}.av.nstep = 1.5e6; g{1}.av.nstep_save_start = 5e5; 
        
        % Submit the job
        job.name = [des_name 'E_O' num2str(n)];
        Jobs{end+1,1} = ts_sbatch(job,g);
        
    end
    
    % Monitor all jobs until completed - allow reversed flow in transients
    ts_smonitor(Jobs,0,1); 
    
end


%% Calculate actual area ratios and plot pressure through the fan
if plot_pressure == 1

    % Figure window for pressure distributions
    figure(); hold on; grid on; box on; xlabel('Axial Distance / m'); ylabel('Pressure Coefficient');
    cols = lines(length(leg));
    
    % Block ids
    bid_noz = 16; bid_bl = 15; bid_in = 2; bid_sta = 7;
    
    % Read in all cases and plot nacelle pressure distributions
    for n = 1:length(leg)
    
        % Read the case
        if plot_avg == 1
            g = ts_read_hdf5([dr.ts des_name  '_avg.hdf5']);
        else
            g = ts_read_hdf5([dr.ts des_name '.hdf5']);
        end
        
        % Reference pressures based on flight speed
%         Po = g{end}.patch{6}.pstag_free; P = g{end}.patch{6}.pstat_free;
        
        % Reference pressures local to IPM lip & nose bullet
        inlet = ts_structured_cut(g,bid_in,1,1,1,'en',1,'en');
        Po = ts_mass_average(inlet,'Po',3); P = ts_area_average(inlet,'P',3);
        
        % Nacelle BL cut
        C{1} = ts_structured_cut(g,bid_bl,1,'en',1,1,1,'en');

        % Intake casing cut
        C{2} = ts_structured_cut(g,bid_in,1,'en','en','en',1,'en');
        
        % Nozzle casing cut
%         C{3} = ts_structured_cut(g,bid_noz,1,'en','en','en',1,'en');
        
        % Tail cone cuts
%         C{4} = ts_structured_cut(g,bid_noz,1,'en',1,1,1,'en');
%         C{5} = ts_structured_cut(g,bid_sta,1,'en',1,1,1,'en');

        % Plot pressure coefficients
        for m = 1:length(C)
            c = C{m};
            l(n) = plot(c.x(:,1),(P - ts_area_average(c,'P',2)) / (Po - P),'-','color',cols(n,:));
        end
        
        % Calculate area ratio from axial velocity ratio in the nozzle
        Vx_3 = ts_area_average(ts_structured_cut(g,bid_noz,1,1,1,'en',1,'en'),'Vx',3);
        Vx_4 = ts_area_average(ts_structured_cut(g,bid_noz,'en','en',1,'en',1,'en'),'Vx',3);
        sig(n) = Vx_3 / Vx_4;
        
    end
    
    % Add legend
    legend(l,leg)
    
end

    
%% Plot characteristics and calculate different fan operating points
if plot_chics == 1
    
        % Figure window for geometry comparisons
    h = figure('position',[1 41 1920 970]); hold on; grid on; box on; cols = lines(7); 
    axis equal; axis([-0.03 0.18 0 0.11]);

    % Define cut planes for evaluation
    coords = [2 7 ; 7 19];

    % Nine colours
    cols = [lines(7) ; 0.6 0.6 0.6 ; 0 0 0]; mar = {'x' '+' 'v' 'o' 's' 'p'};

    % Figure windows
    h.chic = figure(); hold on; grid on; box on; 
    xlabel('Flow Coefficient'); ylabel('Total-Total Pressure Rise Coefficient');
    h.eta = figure(); hold on; grid on; box on; 
    xlabel('Flow Coefficient'); ylabel('Total-Total Fan Efficiency');
    h.split = figure(); hold on; grid on; box on;
    xlabel('Flow Coefficient'); ylabel('Front Power / Rear Power Split');
    h.whirl = figure(); hold on; grid on; box on;
    xlabel('Flow Coefficient'); ylabel('Absolute Whirl Angle / ^\circ');    
    
    % Process all designs in test matrix
    phi_des = phi; nc = zeros(size(phi_des)); o = 1;
        
    % Process the chic
    [inlets,outlets,c,h_temp] = ts_plot_chic(dr.ts,[des_name 'E_M'],[],1,coords); 
    close(h_temp);

    % Plot the chic
    figure(h.chic); plot(c.phi,c.Cptt,'.-','color',cols(o,:));
    figure(h.eta); plot(c.phi,c.eta_poly,'.-','color',cols(o,:));

    % Plot exit contours for selected designs
    if numel(phi_des) < 4
        figure; is = unique(round(linspace(1,length(outlets),4)));
        for i = 1:length(is)
            a = subplot_tight(1,4,i); hold on; axis equal; set(a,'visible','off');
            cut = outlets{is(i)}; v = linspace(0.9,1,20);
            contourf(cut.z,cut.y,cut.sfunc,v); caxis(v([1 end]));
        end
    end

    % Plot power splits and absolute whirl angle changes
    RWx_chic = ones(size(inlets)); Alpha_2_chic = zeros(size(inlets)); 

    % Check the length of the chic
    nc(o) = length(inlets);
    if nc(o) > 1

        % Calculate tip Mach number for all points
        M_chic = zeros(size(inlets)); T_chic = zeros(size(inlets)); 
        omega_chic = pi * inlets{1}.bv.rpm / 30;
        for n = 1:length(inlets)
            M_chic(n) = max(ts_mass_average(inlets{n},'M_rel',2));
            T_chic(n) = ts_mass_average(inlets{n},'T',3);
        end
        M_chic = M_chic([1 1:end end]); T_chic = T_chic([1 1:end end]); 

        % Extend the chic to the axes
        phi_chic = c.phi; Cptt_chic = c.Cptt; 
        phi_chic = [interp1(Cptt_chic(1:2),phi_chic(1:2),0,'linear','extrap') ; phi_chic ; 0];
        eta_chic = [interp1(c.phi(1:2),c.eta_poly(1:2),phi_chic(1),'linear','extrap') ; ...
            c.eta_poly ; 0];
        Cptt_chic = [0 ; Cptt_chic ; 0]; RWx_chic = RWx_chic([1 1:end end]);

        % Plot the extended chic
%             figure(h.chic); plot(phi_chic,Cptt_chic,'--','color',cols(o,:))
%             figure(h.eta); plot(phi_chic,eta_chic,'--','color',cols(o,:))

        % Fan disc area
        A = pi * (inlets{1}.r(end,1)^2 - inlets{1}.r(1,1)^2);

        % Calculate fan operating point at all vehicle operating points
        eta_f = zeros(size(Vf)); eta_p = zeros(size(Vf)); M_tip = zeros(size(Vf)); 
        phi = zeros(size(Vf)); Cptt = zeros(size(Vf)); omega = zeros(size(Vf)); T = zeros(size(Vf));
        RWx = zeros(size(Vf));
        for n = 1:length(Vf)

            % Air conditions at altitude
            [Ta,Pa,ro,~] = Altitude(alt(n));

            % Calculate axial velocity from thrust quadratic formula 
            Vx = (ro * A * Vf(n) + (ro^2 * A^2 * Vf(n)^2 + 4 * ro * A * T_des(n) / sig(o))^0.5) / ...
                (2 * ro * A / sig(o));
            Vj = Vx / sig(o);

            % Bisect to find non-dimensional operating point
            phi_min = phi_chic(end); phi_max = phi_chic(1);
            for m = 1:15

                % Current flow coefficient and blade speed
                phi(n) = 0.5 * (phi_min + phi_max);
                U_1 = Vx / phi(n);

                % Pressure rise coefficient and speed
                Cptt(n) = interp1(phi_chic,Cptt_chic,phi(n),'linear','extrap');
                U_2 = ((Vj^2 - Vf(n)^2) / Cptt(n))^0.5;

                % Update boundaries on flow coefficient
                if U_1 > U_2
                    phi_min = phi(n);
                else
                    phi_max = phi(n);
                end

            end

            % Calculate fan and propulsive efficiencies
            eta_f(n) = interp1(phi_chic,eta_chic,phi(n),'linear','extrap');
            if Vf(n) > 0
                eta_p(n) = 2 * Vf(n) / (Vf(n) + Vj);
            end

            % Compressor inlet and outlet conditions
            Po_1 = Pa + 0.5 * ro * Vf(n)^2; Po_2 = Pa + 0.5 * ro * Vj^2;
            To_1 = Ta; To_2 = To_1 * (Po_2 / Po_1)^((ga-1)/(ga * eta_f(n)));

            % Calculate power requirements, shaft speed and torque
            Wx = ro * A * Vx * cp * (To_2 - To_1);
            r = ts_area_average(inlets{1},'r',3); 
            omega(n) = U_1 / r; T(n) = Wx / (2 * omega(n));

            % Calculate figure of merit
            if Vf(n) == 0
                eta_p(n) = F(n)^1.5 / (Wx * (2 * ro * A)^0.5); 
            end

            % Calculate motor power split
            RWx(n) = interp1(phi_chic,RWx_chic,phi(n),'linear','extrap');

            % Plot the operating point on the chic
            figure(h.chic); plot(phi(n),Cptt(n),'.','color',cols(o,:),'marker',mar{n})
            figure(h.eta); plot(phi(n),eta_f(n),'.','color',cols(o,:),'marker',mar{n})
            figure(h.split); plot(phi(n),RWx(n),'.','color',cols(o,:),'marker',mar{n})

            % Print motor requirements
            fprintf('Design %d\tPoint %d\tTotal Power = %.1f kW\tTotal Speed = %.1f rad/s\tSplit = %.1f\n',...
                [o n Wx / 1000 omega(n) RWx(n)]); 

            % Scale Mach number to actual value at the correct conditions
            M_tip(n) = interp1(phi_chic,M_chic,phi(n),'linear','extrap') * omega(n) * ...
                interp1(phi_chic,T_chic,phi(n),'linear','extrap')^0.5 / (omega_chic * Ta^0.5);

            % Print overall efficiency for cruise and climb cases
            fprintf('Design %d\tPoint %d\tFan Eta = %.1f\tProp. Eta = %.1f\tOv. Eta = %.1f\n',...
                [o n eta_f(n)*100 eta_p(n)*100 eta_f(n)*eta_p(n)*100]);   


        end

        % Save the chic
        save([dr.ts des_name '_chic.mat'],'c','phi','Cptt','eta_f')

        % Calculate stall margin at static conditions
        sm = (phi(1) - phi_chic(end-1)) / (phi(end) - phi_chic(end-1));  

        % Record all variables
        clear c; c.eta_f = eta_f(end); c.eta_p = eta_p(end); c.eta_ov = c.eta_f * c.eta_p;
        c.eta_av = mean(eta_f); c.M_tip = M_tip(1); c.sm = sm; c.eta_cl = eta_f(2) * eta_p(2);
        c.omega_cr = omega(end); c.omega_cl = omega(2); c.omega_st = omega(1); c.Mf = eta_p(1);
        c.T_cr = T(end); c.T_st = T(1); c.T_cl = T(2);

        % Record all data in arrays, preallocate on first case in loop
        scalnames = fieldnames(c);
        for v = 1:length(scalnames)
            if o == 1; p.(scalnames{v}) = nan(size(phi_des)); end;
            p.(scalnames{v})(o) = c.(scalnames{v}); 
        end

        % Print blade number
        disp([inlets{1}.bv.nblade outlets{1}.bv.nblade]);

    end        
    
    % Calculate fan pressure rise from external calc
    coords = [2 7 ; 7 19];
    [inlets,outlets,c,h_temp] = ts_plot_chic(dr.ts,[des_name 'E_O'],[],plot_avg,coords,[],1); 
    
    % Calculate fan flow coefficient at stator inlet
    coords = [9 2 ; 7 19];
    [~,~,c_phi,~] = ts_plot_chic(dr.ts,[des_name 'E_O'],[],plot_avg,coords,[],1); 
    
    % Plot operating points
    cols = lines(7); mar = {'x' '+' 'v' 'o'};
    for n = length(inlets):-1:1
        figure(h.chic); plot(c_phi.phi(n),c.Cptt(n),'.','color',cols(2,:),'marker',mar{end-n+1});
        figure(h.eta); plot(c_phi.phi(n),c.eta_poly(n),'.','color',cols(2,:),'marker',mar{end-n+1});
    end
    
    % Add legends
    figure(h.chic); legend(['2 Rows' leg leg]);
    figure(h.eta); legend(['2 Rows' leg leg]);
    
    % Plot exit swirl
    figure(); hold on; grid on; box on; xlabel('Exit Swirl / ^\circ'); ylabel('Span');
    for n = 1:length(inlets)
        plot(ts_mass_average(outlets{n},'Alpha',2),ts_mass_average(outlets{n},'r_nondim',2))
    end
    legend(flip(leg,2));
    
end


%% Plot thrust from nacelle cases and compare against target values
if plot_thrust == 1
    
    % Block ids and indices to assemble control volume
    bid_far = 18; bid_jet = 17; 
%     i_in = 7; i_out = 7;
    i_in = 41; i_out = 72;
    
    % Loop over all operating points
    for o = 1:length(leg)
    
        % Read the case
        g = ts_read_hdf5([dr.ts des_name '_O' num2str(o) '_avg.hdf5']);
%         g = ts_read_hdf5([dr.ts des_name 'E_O' num2str(o) '.hdf5']);
    
        % Take cuts around propulsor
        j = g{bid_far+1}.attribute.nj-i_out;
        i = g{bid_jet+1}.attribute.ni-i_out;
        C{1} = ts_structured_cut(g,bid_far,i_in,i_in,1,j,1,'en');
        C{2} = ts_structured_cut(g,bid_far,i_in+1,g{bid_far+1}.attribute.ni-1,j,j,1,'en');
        C{3} = ts_structured_cut(g,bid_jet,i,i,1,'en',1,'en');

        % Join cuts together
        varnames = fieldnames(C{1}); flips = [3 3 1];
        for v = 1:length(varnames)
            if isnumeric(C{1}.(varnames{v})) == 1
                c.(varnames{v}) = [];
                for n = 1:length(C)
                    c.(varnames{v}) = [c.(varnames{v}) ; flip(C{n}.(varnames{v}),flips(n))];
                end
            end
        end
        
        % Record other variables
        c.av = C{1}.av; c.bv = C{1}.bv; c.attribute = C{1}.attribute;
        
        % Calculate diagonal lengths
        c = ts_secondary(c); c.t = c.rt ./ c.r;
        t_av = 0.25 * (c.t(2:end,2:end) + c.t(1:end-1,1:end-1) + ...
            c.t(1:end-1,2:end) + c.t(2:end,1:end-1));
        x1 = c.x(2:end,2:end) - c.x(1:end-1,1:end-1); x2 = c.x(1:end-1,2:end) - c.x(2:end,1:end-1);
        r1 = c.r(2:end,2:end) - c.r(1:end-1,1:end-1); r2 = c.r(1:end-1,2:end) - c.r(2:end,1:end-1);
        rt1 = c.r(2:end,2:end) .* (c.t(2:end,2:end) - t_av)...
            - c.r(1:end-1,1:end-1) .* (c.t(1:end-1,1:end-1) - t_av);
        rt2 = c.r(1:end-1,2:end) .* (c.t(1:end-1,2:end) - t_av)...
            - c.r(2:end,1:end-1) .* (c.t(2:end,1:end-1) - t_av);

        % Correct velocities on walls
        q = c.mwall == 0; c.rovx(q) = 0; c.rovr(q) = 0;

        % Area components
        Ax = 0.5 * (r1.*rt2 - r2.*rt1);
        Ar = 0.5 * (x2.*rt1 - x1.*rt2);

        % Calculate pressure force components
        P = 0.25 * (c.P(1:end-1,1:end-1) + c.P(2:end,1:end-1) + c.P(1:end-1,2:end) + c.P(2:end,2:end));
        F = P .* Ax;

        % Calculate mass flux
        rovx_av = 0.25 * (c.rovx(1:end-1,1:end-1) + c.rovx(2:end,1:end-1) + ...
            c.rovx(1:end-1,2:end) + c.rovx(2:end,2:end));
        rovr_av = 0.25 * (c.rovr(1:end-1,1:end-1) + c.rovr(2:end,1:end-1) + ...
            c.rovr(1:end-1,2:end) + c.rovr(2:end,2:end));
        Vx_av = 0.25 * (c.Vx(1:end-1,1:end-1) + c.Vx(2:end,1:end-1) + ...
            c.Vx(1:end-1,2:end) + c.Vx(2:end,2:end));
        m = Ax.*rovx_av + Ar.*rovr_av;

        % Calculate axial momentum flux
        p = m .* Vx_av;

        % Calculate stagnation enthalpy flux
        To_av = 0.25 * (c.To(1:end-1,1:end-1) + c.To(2:end,1:end-1) + ...
            c.To(1:end-1,2:end) + c.To(2:end,2:end));    
        Ho = c.av.cp * m .* To_av;

        % Sum forces, momentum and enthalpy
        F_total = sum(F(:)); p_total = sum(p(:)); Ho_total = sum(Ho(:));

        % Calculate thrust and power
        T(o) = double(c.bv.nblade) * (F_total + p_total); 
        Wx(o) = double(c.bv.nblade) * Ho_total;
        
        % Calculate fan power from rotor row alone
        inlet = ts_structured_cut(g,2,1,1,1,'en',1,'en');
        outlet = ts_structured_cut(g,7,'en','en',1,'en',1,'en');
        [ho_in,m_in] = ts_mass_average(inlet,'ho',3);
        ho_out = ts_mass_average(outlet,'ho',3);
        Wx_fan(o) = m_in * double(inlet.bv.nblade) * (ho_out - ho_in);
        rpm(o) = inlet.bv.rpm;
        
    end

    % Plot the assembled cut
    figure(); hold on; grid on; box on; axis equal;
    plot(c.x(:,1),c.r(:,1),'.-');

    % Save thrust data
    save([dr.ts des_name '_thrust.mat'],'T','Wx','Wx_fan')
        
    % Plot thrust at all operating points
    h.thrust = figure(); hold on; grid on; box on;
    xlabel('Operating Point'); ylabel('Thrust / N');
    plot(-T,'.-');
    
    % Plot power at all operating points
    h.power = figure(); hold on; grid on; box on;
    xlabel('Operating Point'); ylabel('Power / kW');    
    plot(Wx_fan/1000,'.-'); 
    plot(Wx/1000,'.-'); 
    
    % Set operating point names
    varnames = fieldnames(h);
    for v = 1:length(varnames)
        figure(h.(varnames{v})); 
        set(gca,'xtick',1:length(leg),'xticklabel',leg,'xticklabelrotation',0)
        a = axis; axis([a(1:2) 0 a(4)]);
    end
    
    % Calculate non-dimensionals for Dave Filingham
    [Ta,Pa,ro,~] = Altitude(alt);
    n = rpm / 60; D_tip = max(inlet.r(:)) + max(outlet.r(:));
    J = Vf ./ (n * D_tip);
    Ct = -T ./ (ro .* n.^2 * D_tip^4);
    Cp = Wx_fan ./ (ro .* n.^3 * D_tip^5);
    
    % Print non-dimensionals
    for n = 1:length(J)
        fprintf('%12s\t%4.3f\t%4.3f\t%4.3f\t%4.1f\t%4.0f\t%4.2f\n',leg{n},...
            [J(n) Ct(n) Cp(n) Wx_fan(n) rpm(n) -T(n)])
    end
    
end


%% Write PARAVIEW files to share with BBSR
if write_pv == 1
    
    % Loop over all operating points
    for n = 1:length(leg)
        
        % Read the case
        if plot_avg == 1
            %g = ts_read_hdf5([dr.ts des_name 'E_O' num2str(n) '_avg.hdf5']);
            g = ts_read_hdf5([dr.ts des_name '_M9200_avg.hdf5']);
        else
            g = ts_read_hdf5([dr.ts des_name 'E_O' num2str(n) '.hdf5']);
        end
 
        % Write secondary flowfield in polars
        ts_export_paraview(g,[dr.pv des_name '_M9200_rotors' num2str(n) '.hdf5'],'HighSpeed',[0:6],1);
        ts_export_paraview(g,[dr.pv des_name '_M9200_stators' num2str(n) '.hdf5'],'HighSpeed',[7:13],1);
        fprintf('exported for Paraview')
    end
    
end


%% Create IBL files for drawing in Creo
if write_ibl == 1
    
    % Load the design
    load([dr.geom des_name 'E_splines.mat'])
    load([dr.geom des_name 'E_curves.mat'])
    load([dr.geom des_name 'E_xr_hub.mat'])
    load([dr.geom des_name 'E_xr_cas.mat'])
    
    % Create IBL files for all rows
    for n = 1:length(B)
        B{n}.xr_hub = xr_hub; B{n}.xr_cas = xr_cas; 
        B{n}.xrrt_F0 = B{n}.xrrt;
        bl_write_ibl([dr.cad des_name '_' num2str(n) '.ibl'],B{n},0,0);
        B{n}.N
    end
    
    % Create IBL files for nacelle geometry
    varnames = fieldnames(c);
    fid = fopen([dr.cad des_name '_curves.ibl'],'w');
    fprintf(fid,'%s\n','Closed Index arclength');
    for v = 1:length(varnames)
        fprintf(fid,'%s\n','begin section');
        fprintf(fid,'%s\n','begin curve');
        fprintf(fid,'%10.11f %10.11f %10.11f\n',...
            [c.(varnames{v}) zeros(size(c.(varnames{v}),1),1)]' * 1000);
    end
    fclose(fid);
    
end
