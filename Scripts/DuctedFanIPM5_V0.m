%   DuctedFanIPM5_V0
%
%   Generate blades, meshes and CFD solutions for 3 row IPM5 unit
%
%   Written by James Taylor & Sam Grimshaw               
%   May 2020

% Clear workspace and load eVTOL environment
clear; close all; clc; eVTOL;

% Script operation
create_geom = 0;
run_des = 0;
plot_chics = 0;
create_cowl = 0;
mesh_cowl = 0;
run_unsteady = 0;
run_rpms = 0;
read_probes = 0;
ensemble_probes = 1;
nc_bling = 0;
write_ibl = 0;

% Current job variables
job.priority = 'JTAYLOR-SL2-GPU';
job.cards = 1;
job.walltime = 4;

% Choose case to design, run or machine
rpm_fac = 1.8;
des_name = 'IPM5_PH60D'; rpm = 8022 * 1.8;
% des_name = 'IPM5_PH75D'; rpm = 6464 * 1.8;
% des_name = 'IPM5_PH90D'; rpm = 5641 * 1.8;

% Flight operating points
Vf = [1 10 20 30]; alt = zeros(size(Vf)); leg = {'1' '10' '20' '30'}; T_des = [50 40 30 30];

% Choose row to machine
geom_name = 'Rotor';
% geom_name = 'OGV';


%% Make blades to specified angle distributions from DuctedFanDesign
if create_geom == 1

    % Non-dimensional gapping and thicknesses 
    dx_gap = 0.6; dx_duct = 0.5; x_in = 0; thick_max = 0.06; chord_min = 0.80; chord_fac = 1;
    
    % Settings for printing with thicker blades
    thick_max = 0.08; thick_te = 0.6e-3;
    
    % Load the baseline design
    datum_name = 'NSUB85S100DF43HT0.31IGV22-8'; 
    load([dr.geom datum_name '_splines.mat']); B_datum = B;
    load([dr.nc 'otor.mat']); B_datum{2} = b;
    load([dr.nc 'OGV.mat']); B_datum{3} = b;
    
    for n = 1:3
        B_datum{n} = bl_construct_blade(B_datum{n});
    end
    
    % Load blade design from velocity triangles
    load([dr.geom des_name '.mat']);

    % Donor blade parameters
    b.dcam_le = -1.1 * ones(1,4); b.dcam_te = -0.4 * ones(1,4); b.qcam = 2 * ones(1,4); 
    b.thick_te = 0.25 * ones(1,4); b.s_thick_max = 0.3 * ones(1,4); 
    b.rad_thick_max = 0.2 * ones(1,4); b.rad_le = 7 * ones(1,4); b.wedge_te = 60 * ones(1,4);

    % Lean in rotor blades to improve frequency response
    lean = [-30 -10 10 30 ; 10 -25 -15 30] * 1e-4; 
    sweep = [0 30 45 20 ; 0 15 15 0] * 1e-4;

    % Meridional parameters
    A = pi * (d.g.rcas^2 - d.g.rhub^2); r = 0.5 * (d.g.rcas + d.g.rhub);
    b.A_le = A; b.A_te = A; b.r_le = r; b.r_te = r; b.x_ref = 0;
    b.drdx_hub_le = 0; b.drdx_hub_te = 0; b.drdx_cas_le = 0; b.drdx_cas_te = 0;

    % Average chords
    chord_av = mean(d.g.c,1);

    % Loop over both blades of fan design
    d.g.chi(:,:,1) = - d.g.chi(:,:,1);
    r_nondim = (d.g.r - d.g.r(1)) / (d.g.r(end) - d.g.r(1));
    for n = 1:2

        % Set minimum chord
        d.g.c(:,n) = max(d.g.c(:,n),chord_min * chord_av(n));

        % Scale thickness for machining purposes
        b.thick_max = [1.5 1 0.9 0.9 1 1.5] * thick_max * chord_av(n);
        b.thick_te =  [1.05 1 1.0 1.0 1 1.05] * 0.32 ;
        if n == 2
            b.thick_max = b.thick_max.*[1.25 1.15 1 1 1 1];
            b.thick_te = b.thick_te.*[1.25 1.15 1 1 1 1];
        elseif n == 1
            b.thick_te(1) = 0.4; b.thick_te(2) = 0.34;
            b.thick_max(1) = (b.thick_max(1) + 1e-3)*1.15; 
            b.thick_max(2) = (b.thick_max(2) + 0.5e-3).*1.10;
        end
        
        % Increase minimum trailing edge thickness for printing
        thick_abs = b.thick_max .* b.thick_te; b.thick_te = b.thick_te .* max(thick_te ./ thick_abs,1);

        % Set stacking
        b.lean = lean(n,:); b.sweep = sweep(n,:);

        % Fit with blade design from DuctedFanDesign
        b.tchord = spap2(1,4,r_nondim,d.g.c(:,n) * chord_fac);
        b.chi_le = spap2(1,4,r_nondim,d.g.chi(:,1,n));
        b.chi_te = spap2(1,4,r_nondim,d.g.chi(:,2,n));
        b.N = d.g.N(n);
        
        % Reduce rotor hub and stator casing chord
%         if n == 1
%             b.tchord = spap2(1,4,r_nondim,...
%                 fnval(b.tchord,r_nondim) .* interp1([0 0.5 1],[0.9 1 1],r_nondim,'pchip'));
%         end
        if n == 2
            b.tchord = spap2(1,4,r_nondim,...
                fnval(b.tchord,r_nondim) .* interp1([0 0.5 1],[1 0.93 0.85],r_nondim,'pchip'));
        end

        % Fit splines to numeric blade parameters
        B{n} = bl_spline_fit(b);

    end

    % Copy section geometry for IGVs in normally rotating machine
    B = B([2 1:2]); n = 1; B{n}.N = 4; 
    B{n}.lean = [0 0] * 1e-3; B{n}.chi_le = zeros(1,4); 
    B{n}.chi_te = spap2(1,4,r_nondim,d.a.alpha(:,1,1) - 1);
    B{n}.tchord = spap2(1,3,[0 0.5 1],[16 14 18] * 1e-3);
    B{n}.thick_max = [1.5 1 0.9 0.9 1 1.5] * thick_max * mean(chord_av(2));

    % Fit the splines to spanwise numbers
    B{n} = bl_spline_fit(B{n});
    N = [B{1}.N B{2}.N B{3}.N];

    % Record blade speed
    I.rpm = [0 d.d.rpm1 0];

    % Get leading and trailing edge coordinates
    chord = zeros(size(B));
    for n = 1:length(B)
        B{n} = bl_construct_blade(B{n},0,[],[],0,0);
        i_te = round(size(B{n}.xrrt,1) / 2);
        B{n}.i_te = i_te;
        B{n}.i_le = 1;
        x_le(:,n) = B{n}.xrrt(1,:,1); x_te(:,n) = B{n}.xrrt(i_te,:,1);
        c_x(:,n) = x_te(:,n)-x_le(:,n);
        c_x_max(n) = max(c_x(:,n));
        chord(n) = mean(B{n}.tchord.coefs);
    end

    % Set blade axial positions to be the same as the original datum design
    B{1}.x_ref = B_datum{1}.x_ref;
    B{2}.x_ref = B_datum{2}.x_ref - 1e-3;
    B{3}.x_ref = B_datum{3}.x_ref - 0.5e-3;
    B{1}.x_in = B_datum{1}.x_in; B{3}.x_out = B_datum{3}.x_out+10e-3;
    
    % Define fillet parameters
    s_cl = [0 1]; rad = [1.7 1.7] * 1e-3;
    for n = 1:3
        F{n}.hub.s1 = s_cl; F{n}.hub.s2 = s_cl; F{n}.hub.rad1 = rad; F{n}.hub.rad2 = rad; 
        F{n}.cas.s1 = s_cl; F{n}.cas.s2 = s_cl; F{n}.cas.rad1 = rad; F{n}.cas.rad2 = rad; 
        F{n}.flag = [1 1]; F{n}.tol = 1e-6;
    end
    
    % Create machine endwalls
    m = bl_construct_mer(B); xr_hub = m.xr_hub; xr_cas = m.xr_cas;

    % Recreate geometry
    for n = 1:length(B)
        
        % Make blade
        B{n} = bl_construct_blade(B{n},0,[],[],0,0);
        
        % Flip to correct orientation
        j_mid = 22; i_te = 301;
        if I.rpm(n) ~= 0 && B{n}.xrrt(1,j_mid,3) < B{n}.xrrt(i_te,j_mid,3)
            % Flip rotor direction
            B{n}.xrrt(:,:,3) = - B{n}.xrrt(:,:,3);
            disp('Flipping rotor');
        elseif I.rpm(n) == 0 && B{n}.xrrt(1,j_mid,3) > B{n}.xrrt(i_te,j_mid,3)
            if n == 1
                disp('Not flipping IGV'); 
            else
                % Flip stator direction
                B{n}.xrrt(:,:,3) = - B{n}.xrrt(:,:,3);
                disp('Flipping stator');
            end
        end
        
    end
    
    % Save all blades
    save([dr.geom des_name '_splines.mat'],'B','I');
            
    % Edit individual blades in preparation for machining
    for n = 1:length(B)
                
        % Add the fillet
        B{n}.xrrt = bl_construct_fillet(B{n}.xrrt,F{n},xr_hub,xr_cas,0);
        
        % Set meridional extents
        if n == 3
            B{n}.x_hub_le = 62.5/1000; B{n}.r_hub_le = 44/2/1000;
            B{n}.x_hub_te = 100.5/1000; B{n}.r_hub_te = 44/2/1000;
            B{n}.x_cas_le = 72.5/1000; B{n}.r_cas_le = 71.5/1000;
            B{n}.x_cas_te = 100.5/1000; B{n}.r_cas_te = 71/1000;
        elseif n == 2
            B{n}.x_hub_le = 32.5/1000; B{n}.r_hub_le = 44/2/1000;
            B{n}.x_hub_te = 60/1000; B{n}.r_hub_te = 44/2/1000;
            B{n}.x_cas_le = 33.2/1000; B{n}.r_cas_le = 143.6/2/1000;
            B{n}.x_cas_te = 56.2/1000; B{n}.r_cas_te = 143.6/2/1000;
        end
        B{n} = bl_construct_mer(B{n});
        
        % Save for machining
        if n == 1
            blade = '_IGV.mat';
        elseif n == 2
            blade = '_Rotor.mat';
        elseif n == 3
            blade = '_OGV.mat';
        end
        b = B{n}; save([dr.nc des_name blade],'b')
    end

    % Plot all splines
    h = figure(); cols = lines(length(B));
    for n = 1:length(B)
        bl_spline_plot(B_datum{n},h,0.5 * cols(n,:) + 0.5); 
        bl_spline_plot(B{n},h,cols(n,:)); 
    end

    % Plot all blades
    figure(); hold on; grid on; box on; axis equal; cols = lines(length(B)); view(0,0);
    for n = 1:length(B)
        if n > 1
            
            plot3(B{n}.xr_hub(:,1),zeros(size(B{n}.xr_hub,1),1),B{n}.xr_hub(:,2),'color',cols(n,:));
            plot3(B{n}.xr_cas(:,1),zeros(size(B{n}.xr_cas,1),1),B{n}.xr_cas(:,2),'color',cols(n,:));
        end
        plot3(B{n}.xrrt(:,2:2:end,1),-B{n}.xrrt(:,2:2:end,3),B{n}.xrrt(:,2:2:end,2),'k-');
        surf(B{n}.xrrt(:,:,1),-B{n}.xrrt(:,:,3),B{n}.xrrt(:,:,2),...
            'linestyle','none','facecolor',cols(n,:));
%         surf(B_datum{n}.xrrt(:,:,1),-B_datum{n}.xrrt(:,:,3),B_datum{n}.xrrt(:,:,2),...
%             ,'linestyle','none','facealpha',0.7,'facecolor',cols(n,:));
    end
end


%% Create mesh for the design and run design points
if run_des == 1
    
    % Load the design
    load([dr.geom des_name '.mat']);
    
    % Flight condition
    T = d.d.atm.T; P = d.d.atm.p; ro = d.d.atm.ro; V = 26.4; M1 = V / (ga * R * T)^0.5;

    % Calculate stagnation quantities
    Comp_Tables; Po = P / interp1(Mach,P_Po,M1); To = T / interp1(Mach,T_To,M1);
    save([dr.ts 'inlet' '.mat'],'T','P','ro','V','M','Po','To')

    % Generate inlet boundary condition with 2% high boundary layers
    r = [0 0.02 0.98 1]; Po = [P Po Po P]; ni = 401; r_new = linspace(0,1,ni)';
    A = [r_new To * ones(ni,1) interp1(r,Po,r_new,'pchip') zeros(ni,3)];
    dlmwrite([dr.ts 'inlet_ipm5_2percent.dat'],A,'\t',1,0);

    % Set job name
    job.name = des_name; job.run_sa = 1;

    % Check if the job has already been run
    if isempty(dir([dr.ts job.name '_avg.*'])) == 1

        % Load design
        load([dr.geom des_name '_splines.mat']); 
        I.inlet = [dr.ts 'inlet_ipm5_2percent.dat'];

        % Mesh parameters for all rows
        M.stin = [0 1 0]; M.stout = [0 -1 0]; M.n_aid = [25 33 25]; M.n_aou = [25 33 25];
        M.n_bl = [65 89 65]; M.n_bu = [65 89 65]; M.n_so = [9 17 9]; M.n_zi = [17 9 9]; 
        M.n_zo = [9 9 25]; M.ui = 0; M.uo = 0; M.skew = 5;
        
        % Turn off fillets and gaps
        for n = 1:3; F{n}.flag = [0 0]; F{n}.tol = 1e-6; D{n}.flag = [0 0]; end; 

        % Target mass flow 20% higher than design
        I.pout = 95e3; I.k = 50; mass_fac = 1.1;
        r_mean = ((B{1}.xr_cas(1,2)^2 + B{1}.xr_hub(1,2)^2) / 2)^0.5;
        Vx_mean = d.d.phi * r_mean * pi * (I.rpm(2) - I.rpm(3)) / 30;
        I.m = d.d.atm.ro * Vx_mean * B{1}.A_le * mass_fac; 
        
        % Remove IGV row
        B = B(2:3); B{1}.x_in = 0.02; I.rpm = I.rpm(2:3); F = F(2:3); D = D(2:3);
        varnames = fieldnames(M); 
        for v = 1:length(varnames)
            if length(M.(varnames{v})) == 3
                M.(varnames{v}) = M.(varnames{v})(2:3);
            end
        end

        % Set number of steps and initial guess
        I.nstep = 150000; I.poisson_nstep = 10000; I.guess = [];
        I.guess = [dr.ts 'IPM5_PH85B_avg.hdf5'];
        M.dx_row = 2e-3;

        % Mesh and run the design point in TURBOSTREAM
        job = bl_create_machine(B,dr,[job.name '_'],job,I,M,F,D,1,1);
        
    end
    
    % Run the chic
    ts_run_chic(job,[],1,[],0,0.002)

end


%% Plot characteristics and calculate different fan operating points
if plot_chics == 1
    
    % Operating speeds and thrusts - V19
    Vf = [0 26.4]; F = [50 6.8]; 
    alt = [0 0] / 3.281; mar = {'o' 's' 'v' '+' 'x' 'p'};
    
    % Run number
    run_no = 1;

    % Names to plot for designs with IGVs
%     new_name = {'IPM5_PH70' 'IPM5_PH85' 'IPM5_PH100'}; test_name = 'IPM5_PH'; coords = [2 7 ; 14 19];
%     phi = [0.7 0.85 1]; sig = [1 1 1];
    
    % Names to plot for designs with no IGVs
%     new_name = {'IPM5_PH70B' 'IPM5_PH85B' 'IPM5_PH100B'}; test_name = 'IPM5_PH'; coords = [2 7 ; 7 19];
%     phi = [0.7 0.85 1]; sig = [1 1 1];

    % Names to plot for designs with no IGVs and increased loading
%     new_name = {'IPM5_PH70C' 'IPM5_PH85C' 'IPM5_PH100C'}; test_name = 'IPM5_PH'; coords = [2 7 ; 7 19];
%     phi = [0.7 0.85 1]; sig = [1 1 1];
    
    % Names to plot for designs with reduced flow coefficient range
    new_name = {'IPM5_PH60D' 'IPM5_PH75D' 'IPM5_PH90D'}; test_name = 'IPM5_PH'; coords = [2 7 ; 7 19];
    phi = [0.6 0.75 0.9]; sig = [1 1 1];
    
    % Nine colours
    cols = [lines(7) ; 0.6 0.6 0.6 ; 0 0 0];

    % Figure windows
    h.chic = figure(); hold on; grid on; box on; 
    xlabel('Flow Coefficient'); ylabel('Total-Total Pressure Rise Coefficient');
    h.eta = figure(); hold on; grid on; box on; 
    xlabel('Flow Coefficient'); ylabel('Total-Total Fan Efficiency');
    h.whirl = figure(); hold on; grid on; box on;
    xlabel('Flow Coefficient'); ylabel('Absolute Whirl Angle / ^\circ');    
    
    % Process all designs in test matrix
    phi_des = phi; nc = zeros(size(phi_des)); leg = cell(numel(phi_des),1);
    for o = 1:numel(phi_des)
        
        % Process the chic
        phi = phi_des; ducted_desname;
        [inlets,outlets,c,h_temp] = ts_plot_chic(dr.ts,[des_name '_M'],[],1,coords); 
        close(h_temp);
        
        % Skip failed designs
        if isempty(c) == 1; continue; end;
        
        % Legend
        leg{o} = strrep(des_name,'_',' ');;
        
        % Smooth the chics
        c.phi = smooth(c.phi,3); c.Cptt = smooth(c.Cptt,5); c.eta_poly = smooth(c.eta_poly,5);
        
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
            Cptt_chic = [0 ; Cptt_chic ; 0];

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
                Vx = (ro * A * Vf(n) + (ro^2 * A^2 * Vf(n)^2 + 4 * ro * A * F(n) / sig(o))^0.5) / ...
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
                
                % Plot the operating point on the chic
                figure(h.chic); plot(phi(n),Cptt(n),'color',cols(o,:),'marker',mar{n},...
                    'MarkerFaceColor',cols(o,:))
                figure(h.eta); plot(phi(n),eta_f(n),'color',cols(o,:),'marker',mar{n},...
                    'MarkerFaceColor',cols(o,:))
                
                % Print motor requirements
                fprintf('Design %d\tPoint %d\tTotal Power = %.1f kW\tTotal Speed = %.1f rad/s\n',...
                    [o n Wx / 1000 omega(n)]); 
                
                % Scale Mach number to actual value at the correct conditions
                M_tip(n) = interp1(phi_chic,M_chic,phi(n),'linear','extrap') * omega(n) * ...
                    interp1(phi_chic,T_chic,phi(n),'linear','extrap')^0.5 / (omega_chic * Ta^0.5);
                
                % Print overall efficiency for cruise and climb cases
                fprintf('Design %d\tPoint %d\tFan Eta = %.1f\tProp. Eta = %.1f\tOv. Eta = %.1f\n',...
                    [o n eta_f(n)*100 eta_p(n)*100 eta_f(n)*eta_p(n)*100]);   
                
                % Print non-dimensionals for Dave's chart
%                 disp([F(n) / (A * ro * Vf(n)^2) F(n) * Vf(n) / Wx]);

                % Print jet speed and gas conditions for Ricardo
%                 disp([Vj ro To_2])

                % Print motor powers and speeds for DSD
%                 disp([Wx / (1000 * (1 + 1/RWx(n))) 0.5 * (1+rw) * omega(n)])
%                 disp([Wx / (1000 * (1 + RWx(n))) 0.5 * (1-rw) * omega(n)])
%                 disp([Wx / 2000 0.5 * omega(n)])              % CRF Designs
                disp([Wx / 1000 omega(n)])                      % NRF Designs
            
                
            end
            
            % Save the chic
            save([dr.ts des_name '_' num2str(run_no) '_chic.mat'],'c','phi','Cptt','eta_f')
            
            % Calculate stall margin at static conditions
            sm = (phi(1) - phi_chic(end-1)) / (phi(end) - phi_chic(end-1));  
            
            % Record all variables
            clear c; c.eta_f = eta_f(end); c.eta_p = eta_p(end); c.eta_ov = c.eta_f * c.eta_p;
            c.eta_av = mean(eta_f); c.M_tip = M_tip(1); c.sm = sm; c.eta_cl = eta_f(2) * eta_p(2);
            c.omega_cr = omega(end); c.omega_cl = omega(2); c.omega_st = omega(1); c.Mf = eta_p(1);
            c.T_cr = T(end); c.T_st = T(1); c.T_cl = T(2);
            
            % Note torque is 1/2 the real value for CRF fans
            
            % Record all data in arrays, preallocate on first case in loop
            scalnames = fieldnames(c);
            for v = 1:length(scalnames)
                if o == 1; p.(scalnames{v}) = nan(size(phi_des)); end;
                p.(scalnames{v})(o) = c.(scalnames{v}); 
            end
            
            % Print blade number
            disp([inlets{1}.bv.nblade outlets{1}.bv.nblade]);
            
        end        
    end  
    
    % Add the legends
    varnames = fieldnames(h); loc = {'northeast','southeast','northwest','northeast'};
    for v = 1:length(varnames)
        figure(h.(varnames{v}));
        l = flip(get(gca,'children'),1); 
        for a = length(l):-1:1; if strcmp(l(a).Marker,'.') == 0; l(a) = []; end; end;
        legend(l,leg(cellfun(@isempty,leg) == 0),'location',loc{v});
    end
    
    % Save operating point data to file
    save([dr.ts 'chic_' test_name '_' num2str(run_no) '.mat'],'p')
    
end


%% Create cowl geometry and modify meridional lines of blades 
if create_cowl == 1
    
    % Load blade designs
    load([dr.geom des_name '_B1_spline.mat']); B{1} = b;
    load([dr.geom des_name '_B2_spline.mat']); B{2} = b;
    
    % V0 cowl geometry traced from CAD
    r_cas = B{1}.xr_cas(1,2); r_hub = B{1}.xr_hub(1,2);
    r_lip = 0.085; t_te = 1e-3; x_lip = -0.026; xchord = 0.182;
    x_nose = -0.035; e_nose = 2; dx_tail = 0.025; x_tail = 0.18; l_exit = 0.055;
    l_cone = 0.08; thick_max = 0.031; drdx_tail = -0.4; rad_in = 0.58; 
    sig = 1.042; st = 0.75;
    
    % Smooth the transitions at casing end
    l_exit = 0.042; st = 0.8;
    
    % Cowl camberline
    c_out.chi_le = 0; c_out.qcam = 0; c_out.dcam_le = -2; c_out.dcam_te = -0.5;
    
    % Outer thickness distribution
    c_out.rad_le = 2; c_out.thick_max = thick_max; c_out.s_thick_max = 0.43; c_out.thick_te = 0;
    c_out.wedge_te = 70; c_out.rad_thick_max = 0.15; c_out.tchord = xchord;
    
    % Inner thickness distribution
    c_in = c_out; c_in.thick_max = st * c_in.thick_max; c_in.rad_le = c_in.rad_le / st;
    c_in.s_thick_max = 0.2; c_in.rad_thick_max = 0.12; c_in.wedge_te = 50;
    
    % Figure window for geometry comparisons
    h = figure('position',[1 41 1920 970]); hold on; grid on; box on; cols = lines(7); 
    axis equal; axis([-0.04 0.2 0 0.11]); grid minor;
     
    % Plot axis and casing lines
    plot(xlim,[0 0],'k-'); plot(xlim,[1 1] * B{2}.xr_cas(1,2),'k-'); 
    plot(xlim,[1 1] * B{2}.xr_hub(1,2),'k-'); 
    
    % Add screenshot of CAD
    I = imread([dr.fig 'IPM5_V0_Cowl.png']); 
    dy = 0.117; dx = dy * size(I,2) / size(I,1);
    a = image(-0.0495 + [0 dx],-0.013 + [0 dy],I); 
    uistack(a,'bottom'); set(a,'AlphaData',0.5);

    % Plot current blade designs
    for n = 1:2
        mesh(B{n}.xrrt(:,:,1),B{n}.xrrt(:,:,2),zeros(601,43),'edgecolor',cols(1,:))
%         plot3(B{n}.xrrt([1 301],:,1)',B{n}.xrrt([1 301],:,2)',zeros(2,43)','color',cols(1,:))
    end

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
%         plot(xr_cam(:,1),xr_cam(:,2))

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
    xr_in = bl_construct_section(c_in,1,1); xr_in = flip(xr_in(121:301,:),1); 
    
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
    c.xr_tail(:,1) = c.xr_tail(:,1) + dx_tail;
    
    % Plot all external mesh control curves
    figure(h); varnames = fieldnames(c);
    for v = 1:length(varnames)
        plot(c.(varnames{v})(:,1),c.(varnames{v})(:,2),'-','color',cols(5,:));
    end
    
    % Create all meridional lines for blade rows
    b = bl_construct_mer(B,5e-4); xr_hub = b.xr_hub; xr_cas = b.xr_cas;
    
    % Add passage coordinates to all curves
    c.xr_nose = [c.xr_nose ; xr_hub(1,:)];
    c.xr_tail = [xr_hub(end,:) ; c.xr_tail];
    c.xr_int = [c.xr_int ; xr_cas(1,:)];
    c.xr_exit = [xr_cas(end,:) ; c.xr_exit];
    
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
    save([dr.geom des_name 'E_curves.mat'],'c')
    save([dr.geom des_name 'E_splines.mat'],'B','I')
    
    
end


%% Create external mesh around the cowl curves and run in TURBOSTREAM
if mesh_cowl == 1
    
    % Mesh settings
    M.n_bl = 25; M.n_int = 45; M.n_ext = 101; M.n_noz = 41; M.n_out = 73; M.sj_jet = 0.4; M.n_in = 41; 
    M.n_te = 25; M.dt = 10; M.n_k = 9; M.r_min = 0.002; M.d_fac = 1.5; M.d_lip = 0.1; M.d_igv = 0.5;
    
    % Load in existing mesh and curves
    g = ts_read_hdf5([dr.ts des_name 'E_avg.hdf5']);
    load([dr.geom des_name 'E_curves.mat']);
        
    % Create the mesh
    g = mesh_external(g,c,M,1);

    % Loop over all operating points and run cases
    Jobs = cell(0,1);
    for n = 1%:length(Vf)
        
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
        job.name = [des_name 'E_guess_O' num2str(n)];
%         job.name = [des_name 'E_guess_O' num2str(n) '_froz'];
        Jobs{end+1,1} = ts_sbatch(job,g);

    end
    
    % Monitor all jobs until completed - allow reversed flow in transients
    ts_smonitor(Jobs,0,1);    
    
    % Number of hours to run for
    job.walltime = 20;
    
    % Loop over all operating points and run cases
    Jobs = cell(0,1);
    for n = 1:length(Vf)
        
        % Read previous solution
        g = ts_read_hdf5([dr.ts des_name 'E_guess_O' num2str(n) '.hdf5']);
        
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


%% Run full annulus unsteady calculations
if run_unsteady == 1
    
    % Read in previous case
    filename = [des_name 'E_O1'];
    g = ts_read_hdf5([dr.ts filename '.hdf5']);
    
    % Set RPMs
    varnames = {'rpm' , 'rpmi1' , 'rpmi2' , 'rpmj1' , 'rpmj2' , 'rpmk1' , 'rpmk2'};
    for bid = ts_get_bids(g)
        for v = 1:length(varnames)
            if g{bid+1}.bv.(varnames{v}) ~= 0
                g{bid+1}.bv.(varnames{v}) = rpm;
            end
        end
    end
    
    % Define probes
    bid_far = 18; bid_jet = 17; 
    i_in = [7 41]; i_out = [7 65];
    clear F; F.N = 6;
    for n = 1:2
        j = g{bid_far+1}.attribute.nj-i_out(n);
        i = g{bid_jet+1}.attribute.ni-i_out(n);
        F.bid(n*3-2:n*3) = [bid_far bid_far bid_jet]; 
        F.ist(n*3-2:n*3) = [i_in(n) i_in(n) i];
        F.ien(n*3-2:n*3) = [i_in(n) g{bid_far+1}.attribute.ni i];
        F.jst(n*3-2:n*3) = [1 j 1];
        F.jen(n*3-2:n*3) = [j j g{bid_jet+1}.attribute.nj];
        F.kst(n*3-2:n*3) = [1 1 1];
        F.ken(n*3-2:n*3) = [g{bid_far+1}.attribute.nk g{bid_far+1}.attribute.nk ...
            g{bid_jet+1}.attribute.nk];        
    end
    
    % Add the probes
    g = ts_face_add_patch(g,F,'probe');

    % Plot the probes
    ts_check_rpms(g)
    for bid = ts_get_bids(g)
        for pid = ts_get_pids(g,bid)
            if g{bid+1}.patch{pid+1}.attribute.kind == 8
                a = g{bid+1}.patch{pid+1}.attribute;
                i = a.ist+1:a.ien; j = a.jst+1:a.jen;
                plot3(g{bid+1}.x(i,j,1),g{bid+1}.r(i,j,1),100*ones(length(i),length(j)),'b.-')
            end
        end
    end
   
    % Create unsteady mesh
    b_row = ts_get_row_bids(g(1:14),[min(g{1}.r(:)),max(g{1}.r(:))]); b_row{1} = 14:18;
    g = ts_steady_to_unsteady(g,1,1,b_row);    
    
    % Set number of cycles and probe output
    g{1}.av.nstep_cycle = 100; g{1}.av.ncycle = 9*50; g{1}.av.nstep_save_probe = 1;
    job.cards = 4; job.walltime = 36;
    
    % Run frozen rotor first
    g{1}.av.dts = 0; g{1}.av.nstep = 40000;
    job.name = [filename '_froz']; 
    job = ts_sbatch(job,g); ts_smonitor(job);    
    g = ts_read_hdf5([dr.ts job.name '.hdf5']); 
    
    % Set the frequency properly
    g{1}.av.dts = 1; rpm_max = 0;
    for bid = ts_get_bids(g)
        if g{bid+1}.bv.rpm > rpm_max
            rpm_max = g{bid+1}.bv.rpm;
            nblade = double(g{bid+1}.bv.nblade);
        end
    end
    bpf_max = nblade * rpm_max / 60;
    g{1}.av.nstep_cycle = 200; g{1}.av.frequency = bpf_max; g{1}.av.ncycle = nblade*12; 
    g{1}.av.nstep_save_start_probe = g{1}.av.ncycle * g{1}.av.nstep_cycle + 1;
    g{1}.av.nstep_save_start = g{1}.av.ncycle * g{1}.av.nstep_cycle + 1;
    
    % Run the unsteady job
    job.name = [filename '_unst']; 
    job = ts_sbatch(job,g); ts_smonitor(job);    
    
    % Read the case and restart with production settings
    g = ts_read_hdf5([dr.ts job.name '.hdf5']);
    
    % Probe settings for final revolution
    g{1}.av.nstep_save_probe = 1; 
    g{1}.av.nstep_save_start_probe = (g{1}.av.ncycle - nblade) * g{1}.av.nstep_cycle;

    % Job name and submit
    job.name = [filename '_prod']; 
    job = ts_sbatch(job,g); ts_smonitor(job);
    
    % Run again with more revolutions and larger timestep
    g = ts_read_hdf5([dr.ts job.name '.hdf5']);
    g{1}.av.nstep_cycle = 100; g{1}.av.frequency = bpf_max; g{1}.av.ncycle = nblade*20; 
    g{1}.av.nstep_save_start_probe = (g{1}.av.ncycle - nblade) * g{1}.av.nstep_cycle;
    job.name = [filename '_nstep']; 
    job = ts_sbatch(job,g); ts_smonitor(job);           
    
end


%% Read and process unsteady probes
if read_probes == 1
    
    % Set filenames
%     file_prefix = 'IPM5_PH60DE_O1_prod'; nsteps = 13800:13999; 
%     file_prefix = 'IPM5_PH75DE_O1_prod'; nsteps = 15800:15999; 
%     file_prefix = 'IPM5_PH90DE_O1_prod'; nsteps = 21800:21999; 
%     file_prefix = 'IPM5_PH60DE_O1_nstep3'; nsteps = 16100:16799; 
%     file_prefix = 'IPM5_PH75DE_O1_nstep3'; nsteps = 18400:19199; 
    file_prefix = 'IPM5_PH90DE_O1_nstep'; nsteps = 20900:21999; 
    
    % Read in converged file
    g = ts_read_hdf5([dr.ts file_prefix '.hdf5']);
    
    % Export to PARAVIEW
    bp = {'P','ro','Vx','Vr','Vt'};
    ts_export_paraview(g,[dr.pv file_prefix '.hdf5'],bp,[],1,0);
    
    % Process the probes
    nstep_cycle = double(g{1}.av.nstep_cycle);
    ts_export_probes([dr.ts 'Probes/'],file_prefix,bp,g{1}.av.frequency,nstep_cycle,...
        1,[],nsteps,2*pi/9,0,0)
    disp(g{1}.av.frequency);

    % Plot the probes
    figure(); hold on; grid on; box on; cols = hsv(length(nsteps)); o = 1;
    for nstep = nsteps
        load([dr.ts 'Probes/' file_prefix '_probe_' num2str(nstep) '.mat'])
        g = ts_secondary(g); hold on; 
        for n = 50:2:70
            c = g{n};
%             surf(c.z,c.y,c.rovx,'edgecolor','none'); 
            plot(c.t(50,:),c.P(50,:),'color',cols(o,:));
        end
%         drawnow; 
        rovx_sample(o) = c.P(50,4);
        o = o + 1;
%         hold off; drawnow;
    end
    
    % Plot pressure at a single location
    figure(); hold on; grid on; box on;
    plot(nsteps,rovx_sample);
    
end


%% Ensemble average a single blade passing from a revolution of probe data
if ensemble_probes == 1
    
    % Set filenames and step numbers
    file_prefix = 'IPM5_PH60DE_O1_nstep3'; nsteps = 16100:16799; 
%     file_prefix = 'IPM5_PH75DE_O1_nstep3'; nsteps = 18400:19199; 
%     file_prefix = 'IPM5_PH90DE_O1_nstep'; nsteps = 20900:21999; 

    % Read in converged file
    g = ts_read_hdf5([dr.ts file_prefix '.hdf5']);
    nstep_cycle = double(g{1}.av.nstep_cycle); frequency = g{1}.av.frequency;
    
    % Read the first and last probes and store the offset
    varnames = {'ro' 'rovx' 'rovr' 'rorvt' 'roe'}; s = linspace(0,1,length(nsteps));
    load([dr.ts 'Probes/' file_prefix '_probe_' num2str(nsteps(1)) '.mat']); g_first = g;
    load([dr.ts 'Probes/' file_prefix '_probe_' num2str(nsteps(end)) '.mat']); g_last = g;
    for o = 1:length(g)
        for v = 1:length(varnames)
            d{o}.(varnames{v}) = g_last{o}.(varnames{v}) - g_first{o}.(varnames{v});
        end
    end

    % Loop over each position
    N = round((nsteps(end)-nsteps(1)) / nstep_cycle); 
    for n = 1:nstep_cycle
        
        % Loop over all rotors
        for m = 1:N
            
            % Read the probe
            nstep = nsteps(n) + (m-1) * nstep_cycle;
            load([dr.ts 'Probes/' file_prefix '_probe_' num2str(nstep) '.mat'])
            
            % Initialise empty arrays for first position
            if m == 1
                g_av = g;
                for o = 1:length(g)
                    for v = 1:length(varnames)
                        g_av{o}.(varnames{v}) = zeros(size(g_av{o}.(varnames{v}))); 
                    end
                end
            end
            
            % Correct the offset from the current value
            for o = 1:length(g)
                for v = 1:length(varnames)
                    g{o}.(varnames{v}) = g{o}.(varnames{v}) - s(n) * d{o}.(varnames{v}); 
                end
            end
            
            % Add contribution to average
            for o = 1:length(g)
                for v = 1:length(varnames)
                    g_av{o}.(varnames{v}) = g_av{o}.(varnames{v}) + g{o}.(varnames{v}) / N; 
                end
            end            
        end
        
        % Save the probe to file
        g = g_av; save([dr.ts 'Probes/' file_prefix '_av_probe_' num2str(n) '.mat'],'g')
        
    end
    
    % Save hdf5 files
    bp = {'P','ro','Vx','Vr','Vt'};
    ts_export_probes([dr.ts 'Probes/'],[file_prefix '_av'],bp,frequency,nstep_cycle,...
        1,[],1:nstep_cycle,2*pi/N,0,0)
        
    % Read the probes to record single point values
    for n = 1:nstep_cycle
        load([dr.ts 'Probes/' file_prefix '_probe_' num2str(nsteps(n)) '.mat'])
        c = ts_secondary(g{50}); rovx_sample(n) = c.rovx(50,4);       
        load([dr.ts 'Probes/' file_prefix '_av_probe_' num2str(n) '.mat'])
        c = ts_secondary(g{50}); rovx_av(n) = c.rovx(50,4);
    end
    
    % Plot flow property at a single location
    figure(); hold on; grid on; box on;
    plot(rovx_sample);
    plot(rovx_av);
    
end


%% Machine a bladed ring for the rotor or OGV rows
if nc_bling == 1
    
    % Machine settings and compensation for DMG
%     mach.type = 'DMG'; mach.dA = 0.070021;
%     mach.max = [inf inf inf inf inf]; mach.min = [-inf -inf -inf -inf -inf];
% 
%     % Machine origin in blade coordinate system
%     mach = {mach ; mach ; mach}; dx_flip = 0e-3; z_fix = 186.693e-3;
%     mach{1}.zero = [b.x_hub_te + z_fix + dx_flip 0 0]; mach{1}.abc = [0 1 0 ; 0 0 -1 ; -1 0 0];
%     mach{2}.zero = [b.x_hub_le - z_fix 0 0]; mach{2}.abc = [0 1 0 ; 0 0 1 ; 1 0 0];
%     mach{3}.zero = [b.x_hub_te + z_fix 0 0]; mach{3}.abc = [0 1 0 ; 0 0 -1 ; -1 0 0];
    
    % Machine settings and compensation for XYZ
    mach.type = 'XYZ'; mach.dA = 0; mach.long_part = 0;
    mach.max = [inf inf inf inf inf]; mach.min = [-inf -inf -inf -inf -inf];

    % Tool and holder geometries
    tool.BN3P0 = struct('name','BN3P0','r',1.5e-3,'l',3e-3,'nr',[65 6 ; 75 6 ; 75 12.5] * 1e-3);
    tool.EM3P0 = struct('name','EM3P0','r',1.5e-3,'l',3e-3,'nr',[65 6 ; 75 6 ; 75 12.5] * 1e-3);
    tool.EM8P0 = struct('name','QG8','r',4e-3,'l',43e-3,'nr',[0 13.5] * 1e-3);

    % Load the geometry
    load([dr.nc des_name '_' geom_name])
    
    % Choose which geometry to machine
    if strcmp(geom_name,'Rotor') == 1
        
        % Machine origin in blade coordinate system
        mach = {mach ; mach ; mach}; dx_flip = 0e-3;
        mach{1}.zero = [158.138 0 0] * 1e-3; mach{1}.abc = [0 0 -1 ; 1 0 0 ; 0 -1 0]';
        mach{2}.zero = [-65.6380 0 0] * 1e-3; mach{2}.abc = [0 0 1 ; 1 0 0 ; 0 1 0]';
        mach{3}.zero = [158.138 0 0] * 1e-3; mach{3}.abc = [0 0 -1 ; 1 0 0 ; 0 -1 0]';
        
        % Flip the rotor
        b.xrrt(:,:,3) = - b.xrrt(:,:,3);
        
    elseif strcmp(geom_name,'OGV') == 1
        
        % Machine origin in blade coordinate system
        mach = {mach ; mach ; mach}; dx_flip = 0e-3;
        mach{1}.zero = [198.638 0 0] * 1e-3; mach{1}.abc = [0 0 -1 ; 1 0 0 ; 0 -1 0]';
        mach{2}.zero = [-35.6448 0 0] * 1e-3; mach{2}.abc = [0 0 1 ; 1 0 0 ; 0 1 0]';
        mach{3}.zero = [198.638 0 0] * 1e-3; mach{3}.abc = [0 0 -1 ; 1 0 0 ; 0 -1 0]';
        
        % Flip the OGV
        b.xrrt(:,:,3) = - b.xrrt(:,:,3); b.xrrt = flip(b.xrrt,1);
        
    end
    
    % Create endwall surfaces
    b.x_in_hub = b.x_hub_le; b.x_out_hub = b.x_hub_te; 
    b.x_in_cas = b.x_cas_le; b.x_out_cas = b.x_cas_te; 
    b = bl_construct_platform(b,'hub'); b.xrrt_hub = b.xrrt_plat;
    b = bl_construct_platform(b,'cas'); b.xrrt_cas = b.xrrt_plat;
    b.dx_spigot = 0.003;
    
    % Calculate polar coordinates of surfaces
    b.xyz = m_pol2cart(b.xrrt); 
    b.xyz_hub = m_pol2cart(b.xrrt_hub); 
    b.xyz_cas = m_pol2cart(b.xrrt_cas);
    
    % Roughing operations
    edge = tool.EM8P0; edge.S = 7000; edge.F = 1000; edge.ds = 2e-3; 
    pass = tool.EM3P0; pass.S = 8000; pass.F = 500; pass.ds = 1e-3; pass.dz = 4e-3;
    rough = struct('pass',pass,'edge',edge); 
    rough.stock = 0.2e-3; rough.dx_flip = dx_flip; rough.type = 'break'; cycle.rough = rough;

    % Blade and meridional operations
    rough = tool.BN3P0; rough.S = 8000; rough.F = 800; rough.ds = 0.8e-3;
    finish = tool.BN3P0; finish.S = 8000; finish.F = 1000; finish.ds = 0.16e-3;
    blade = struct('rough',rough,'finish',finish);
    blade.stock = 0.1e-3; blade.tilt = 8; cycle.blade = blade; 

    % Fillet operations
    rough = tool.BN3P0; rough.S = 8000; rough.F = 800; rough.ds = 0.8e-3;
    finish = tool.BN3P0; finish.S = 8000; finish.F = 1000; finish.ds = 0.16e-3;
    fillet = struct('rough',rough,'finish',finish);
    fillet.stock = 0.1e-3; fillet.d_bl = 1e-3; fillet.d_mer = 1e-3;
    cycle.fillet = fillet; 

    % Write code
    [xyzbc,time] = g_bling(b,[dr.nc upper([des_name '_' geom_name]) '.MPF'],cycle,mach,1);
    disp(['Estimated run time = ' num2str(time / 60) ' / hours']);

    % Process instructions
    %   Bolt down and run roughing program
    %   Flip over and run back program which finishes the trailing edge
    %   Flip again and run front program which finishes all the passages
    
end


%% Write IBL files for input into Creo
if write_ibl == 1
    
    % Load the geometry
    load([dr.geom des_name '_splines.mat'])
    
    % Write the blade count to file
    dlmwrite([dr.cad lower([des_name '_N.dat'])],[B{2}.N B{3}.N])
    
    % Save to file
    for n = 2:3
        bl_write_ibl([dr.cad lower([des_name '_' num2str(n-1) '.ibl'])],B{n},0,0,22);
    end
    
end





