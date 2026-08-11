%   make_prop_thrust
%
%   Compute propeller thrust, thrust coefficient and efficiency vs blade count
%   from the swept CFD solutions produced by PropCURTIS_bladesweep.m.
%
%   For each blade count it reads the averaged solution and computes, across the
%   blade row (LE face = upstream H-block bid 2 i=end; TE face = downstream
%   H-block bid 0 i=1) over the blade span (hub -> tip Rc):
%
%     Thrust  T = nblade * [ (rovx*Vx + P)*Ax ]_TE - [ ... ]_LE   (rotor axial force)
%     Power   P_shaft = nblade * m * ( ho_TE - ho_LE )            (enthalpy rise)
%     C_T = T / (rho n^2 D^4),  C_P = P_shaft / (rho n^3 D^5)
%     J   = Vf / (n D),         eta = T*Vf / P_shaft
%
%   with n = rpm/60 (rev/s), D = 2*Rc. Thrust definition follows
%   entrop_loss_prop.m / DuctedFanIPM5.m, restricted to the blade annulus.
%
%   Saves the results to <TURBOSTREAM>/<des_base>_thrust_sweep.mat for
%   plot_prop_thrust.m .

clear; close all; clc; eVTOL;

% --- inputs ---
des_base = '02062026PropPhi06HYBRID';
suffix   = 'restart2_O1_avg';          % averaged-solution suffix
counts   = 2:21;                       % rotor blade counts in the sweep
Vf       = 20;                         % flight speed the sweep was run at [m/s]

s = importdata([dr.geom des_base '-mission.json'],'r'); mm = jsondecode(s{1});
Rc = mm.Rc; Rh = mm.Rh; Ndes = mm.NB_r;   % tip radius, hub radius, design blade count

% Rotor solidity s = C/S at the meanline (chord / pitch). The blade sections are
% held fixed across the sweep, so only the pitch S = 2*pi*r_m/N changes and the
% solidity scales linearly with blade count. NOTE mm.sigma is the disk-to-jet
% velocity ratio (Ch.3), NOT solidity - the chord is read from the blade geometry.
r_mean = sqrt((Rc^2 + Rh^2)/2);                       % meanline radius
xrrt = dlmread([dr.geom des_base '.csv'],',',1,0);    % blade point cloud
xrrt = reshape(xrrt(:,[1 3 2]), 428, [], 3);          % [pt, span, (x r rt)]
rsec = squeeze(mean(xrrt(:,:,2),1));                  % mean radius of each section
[~,jm] = min(abs(rsec - r_mean));                     % meanline section
xs = double(xrrt(:,jm,1)); ts = double(xrrt(:,jm,3)); % (axial, r*theta) of section
dx = xs - xs.'; dt = ts - ts.';                       % max chord = max LE->TE distance
Cm = max(sqrt(dx(:).^2 + dt(:).^2));                  % true meanline chord [m]

Tvec = nan(size(counts)); Pvec = nan(size(counts)); Nb = nan(size(counts));
mdot = nan(size(counts)); rho = nan(size(counts)); rpm = NaN;
turn = nan(size(counts)); Yp = nan(size(counts));   % rotor turning [deg] and loss coeff
for k = 1:numel(counts)
    N = counts(k);
    % Accept both the zero-padded naming (B_02_) used by new runs and the legacy
    % un-padded one (B_2_) that existing solutions were written with.
    f = [dr.ts 'B_' num2str(N,'%02d') '_' des_base suffix '.hdf5'];
    if isempty(dir(strrep(f,'.hdf5','.*')))
        f = [dr.ts 'B_' num2str(N) '_' des_base suffix '.hdf5'];
    end
    if isempty(dir(strrep(f,'.hdf5','.*'))); warning('missing %s', f); continue; end
    g = ts_read_hdf5(f);
    nblade = double(g{1}.bv.nblade); Nb(k) = nblade;

    % Rotor speed (max magnitude across blocks)
    rr = []; for bid = ts_get_bids(g); rr(end+1) = g{bid+1}.bv.rpm; end
    rpm = max(abs(rr));

    % Tip j-index (closest grid line to Rc) on the LE and TE faces
    [~,jtL] = min(abs(g{3}.r(end,:,1) - Rc));    % bid 2, i=end (LE)
    [~,jtT] = min(abs(g{1}.r(1,:,1)   - Rc));    % bid 0, i=1   (TE)

    % Blade-span cuts
    cLE = ts_secondary(ts_structured_cut(g, 2, 'en','en', 1, jtL, 1,'en'));
    cTE = ts_secondary(ts_structured_cut(g, 0, 1,1,     1, jtT, 1,'en'));

    % Rotor axial force (thrust)
    Tvec(k) = nblade * (prop_axflux(cTE) - prop_axflux(cLE));

    % Shaft power from the stagnation-enthalpy rise across the rotor
    [ho_LE, m_LE] = ts_mass_average(cLE,'ho',3);
    ho_TE = ts_mass_average(cTE,'ho',3);
    Pvec(k) = m_LE * nblade * (ho_TE - ho_LE);

    % Rotor flow turning and stagnation-pressure loss (relative frame, mass-avg)
    a1 = ts_mass_average(cLE,'Alpha_rel',3);        % relative flow angle in  [deg]
    a2 = ts_mass_average(cTE,'Alpha_rel',3);        % relative flow angle out [deg]
    turn(k) = abs(a1 - a2);                         % rotor deflection Delta-beta
    po1 = ts_mass_average(cLE,'Po_rel',3);          % rel. stagnation pressure in
    po2 = ts_mass_average(cTE,'Po_rel',3);          % rel. stagnation pressure out
    p1  = ts_mass_average(cLE,'P',3);               % static pressure in
    Yp(k) = (po1 - po2) / (po1 - p1);               % Yp = dPo_rel / (Po_rel - P)_in

    % Rotor mass flow + freestream density (for the coefficients)
    [Ax,rovx] = prop_cellax(cTE); mdot(k) = nblade * abs(sum(Ax(:).*rovx(:)));
    up  = ts_structured_cut(g, numel(g)-2, 1,1, 1,'en', 1,'en');   % far-upstream plane
    rho(k) = mean(up.ro(:));

    fprintf('N=%2d  mdot=%.3f  T=%6.2f N  P=%7.1f W  turn=%5.2f deg  Yp=%.4f\n', ...
        N, mdot(k), Tvec(k), Pvec(k), turn(k), Yp(k));
end
if nanmedian(Tvec) < 0; Tvec = -Tvec; end     % thrust positive
Pvec = abs(Pvec);

% --- non-dimensional coefficients ---
n = rpm/60; D = 2*Rc;                         % rev/s and diameter
solidity = Cm * counts / (2*pi*r_mean);       % meanline solidity s = C/S (~ N)
CT  = Tvec ./ (rho .* n^2 .* D^4);            % thrust coefficient
CP  = Pvec ./ (rho .* n^3 .* D^5);            % power coefficient
J   = Vf ./ (n * D);                          % advance ratio (constant across sweep)
eta = Tvec .* Vf ./ Pvec;                     % propulsive efficiency

fprintf('\n  N    T[N]    C_T      C_P      eta     (J=%.3f, rpm=%.0f)\n', J, rpm);
for k = 1:numel(counts)
    fprintf('  %2d  %6.2f  %.4f  %.4f  %.4f\n', counts(k), Tvec(k), CT(k), CP(k), eta(k));
end

save([dr.ts des_base '_thrust_sweep.mat'], ...
     'counts','Tvec','Pvec','CT','CP','J','eta','mdot','rho','Nb','Ndes','rpm','Vf','Rc','D', ...
     'solidity','Cm','r_mean','turn','Yp');
fprintf('\nsaved %s%s_thrust_sweep.mat\n', dr.ts, des_base);


% ================= local functions =================

function val = prop_axflux(c)
    % axial momentum + pressure flux over an i=const cut:  sum( (rovx*Vx + P)*Ax )
    [Ax,rovx,Vx,P] = prop_cellax(c);
    val = sum(sum((rovx.*Vx + P) .* Ax));
end

function [Ax,rovx,Vx,P] = prop_cellax(c)
    % axial cell-area (Ax) and cell-averaged rovx / Vx / P for an i=const cut
    c.t = c.rt ./ c.r;
    t_av = 0.25*(c.t(2:end,2:end)+c.t(1:end-1,1:end-1)+c.t(1:end-1,2:end)+c.t(2:end,1:end-1));
    r1 = c.r(2:end,2:end)-c.r(1:end-1,1:end-1); r2 = c.r(1:end-1,2:end)-c.r(2:end,1:end-1);
    rt1 = c.r(2:end,2:end).*(c.t(2:end,2:end)-t_av) - c.r(1:end-1,1:end-1).*(c.t(1:end-1,1:end-1)-t_av);
    rt2 = c.r(1:end-1,2:end).*(c.t(1:end-1,2:end)-t_av) - c.r(2:end,1:end-1).*(c.t(2:end,1:end-1)-t_av);
    Ax = 0.5*(r1.*rt2 - r2.*rt1);
    if isfield(c,'mwall'); c.rovx(c.mwall==0) = 0; end
    rovx = 0.25*(c.rovx(1:end-1,1:end-1)+c.rovx(2:end,1:end-1)+c.rovx(1:end-1,2:end)+c.rovx(2:end,2:end));
    Vx   = 0.25*(c.Vx(1:end-1,1:end-1)  +c.Vx(2:end,1:end-1)  +c.Vx(1:end-1,2:end)  +c.Vx(2:end,2:end));
    P    = 0.25*(c.P(1:end-1,1:end-1)   +c.P(2:end,1:end-1)   +c.P(1:end-1,2:end)   +c.P(2:end,2:end));
end
