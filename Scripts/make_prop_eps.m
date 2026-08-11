%   make_prop_eps
%
%   Compute the flow-contraction coefficient eps = V_fs / V_j for each blade
%   count of the sweep and append it to the thrust-sweep .mat, so eps can be
%   plotted against solidity (plot_prop_eps.py).
%
%   eps is the Chapter-3 propulsor loading ratio: freestream velocity over
%   developed-jet velocity. It is measured from the CFD exactly as in the
%   chic/jet-sweep pipeline (make_prop_chic.m + make_prop_jet_sweep.m +
%   plot_prop_chic.py):
%
%     V_fs = area-averaged axial velocity on the domain inlet plane (freestream)
%     V_j  = PEAK capture-streamtube axial velocity over a sweep of constant-x
%            planes downstream of the rotor. Behind the rotor the slipstream
%            contracts and accelerates to the developed jet (~0.8 Rc), then mixes
%            out; the peak is the physically meaningful jet. The capture
%            streamtube at each plane is the inner region carrying the rotor's
%            captured mass flow (accumulated hub -> tip to mdot_rotor).
%
%   Run make_prop_thrust.m first (this appends to its .mat).

clear; close all; clc; eVTOL;

des_base = '02062026PropPhi06HYBRID';
suffix   = 'restart2_O1_avg';
counts   = 2:21;

s = importdata([dr.geom des_base '-mission.json'],'r'); mm = jsondecode(s{1});
Rc = mm.Rc;

bid_in  = 2;        % upstream H-block: i=end is the rotor LE face
bid_out = 0;        % downstream H-block: i=1 is the rotor TE face
r_near  = 2.0;      % a block must reach within r_near*Rc of the axis to qualify
i_step  = 2;        % subsample every i_step planes per block

eps  = nan(size(counts));
Vfs  = nan(size(counts));
Vj   = nan(size(counts));
for k = 1:numel(counts)
    N = counts(k);
    f = [dr.ts 'B_' num2str(N,'%02d') '_' des_base suffix '.hdf5'];
    if isempty(dir(strrep(f,'.hdf5','.*')))
        f = [dr.ts 'B_' num2str(N) '_' des_base suffix '.hdf5'];
    end
    if isempty(dir(strrep(f,'.hdf5','.*'))); warning('missing %s', f); continue; end
    g = ts_read_hdf5(f);

    % --- freestream: area-average Vx on the domain inlet plane (kind 0 or 13) ---
    bid_din = []; i_din = [];
    for bid = ts_get_bids(g)
        for pid = ts_get_pids(g,bid)
            kd = g{bid+1}.patch{pid+1}.attribute.kind;
            if (kd == 0 || kd == 13) && isempty(bid_din)
                bid_din = bid; i_din = g{bid+1}.patch{pid+1}.attribute.ist;
            end
        end
    end
    din = ts_structured_cut(g, bid_din, i_din+1, i_din+1, 1,'en', 1,'en');
    din.Vx = din.rovx ./ din.ro;
    Vfs(k) = area_avg(din, 'Vx');

    % --- captured rotor mass flow: blade-span rotor LE plane ---
    [~,jt_in] = min(abs(g{bid_in+1}.r(end,:,1) - Rc));
    ni_in = g{bid_in+1}.attribute.ni;
    inlet = ts_structured_cut(g, bid_in, ni_in, ni_in, 1, jt_in, 1,'en');
    mdot_rotor = axial_massflow(inlet);

    % --- developed jet: peak capture-streamtube Vx over downstream planes ---
    x_te = mean(mean(g{bid_out+1}.x(1,:,:)));
    vj = 0;
    for bid = ts_get_bids(g)
        xb = g{bid+1}.x; rb = g{bid+1}.r;
        if max(xb(:)) <= x_te + 1e-9, continue; end            % must be downstream
        if min(rb(:)) > r_near * Rc, continue; end             % must reach near the axis
        dxi = abs(xb(end,:,:) - xb(1,:,:)); dri = abs(rb(end,:,:) - rb(1,:,:));
        if mean(dxi(:)) < mean(dri(:)), continue; end          % skip if i is not streamwise
        ni = double(g{bid+1}.attribute.ni);
        for i = 1:i_step:ni
            xi = xb(i,:,:); if mean(xi(:)) < x_te - 1e-9, continue; end
            c = ts_structured_cut(g, bid, i, i, 1,'en', 1,'en');
            c.Vx = c.rovx ./ c.ro;
            v = streamtube_vx(c, mdot_rotor);
            if isfinite(v) && v > vj; vj = v; end
        end
    end
    Vj(k) = vj;
    eps(k) = Vfs(k) / Vj(k);
    fprintf('N=%2d  Vfs=%.3f  Vj=%.3f  eps=%.4f\n', N, Vfs(k), Vj(k), eps(k));
end

% --- Chapter-3 new thrust coefficient, freshly measured from the CFD thrust
% budget: C_T = Thrust / (mdot * V_j). This is computed directly from the
% integrated rotor thrust, the rotor mass flow and the measured jet velocity
% (NOT as 1 - eps); it equals 1 - eps only to the extent the momentum relation
% Th = mdot*(V_j - V_fs) holds in the CFD, so the small difference from 1 - eps
% is the residual of that relation. ---
mf = [dr.ts des_base '_thrust_sweep.mat'];
S = load(mf); S.eps = eps; S.Vfs_eps = Vfs; S.Vj_eps = Vj;
S.CT_new = S.Tvec ./ (S.mdot .* Vj);          % Th/(mdot*Vj),   measured
S.CP_new = S.Pvec ./ (S.mdot .* Vj.^2);       % Power/(mdot*Vj^2), measured
save(mf, '-struct', 'S');
fprintf('\nappended eps, CT_new and CP_new to %s\n', mf);


% ================= local functions =================

function a = area_avg(c, var)
    [Ax,Ar,At] = cell_areas(c.x, c.r, c.rt);
    mag = sqrt(Ax.^2 + Ar.^2 + At.^2);
    p = cell_corners(c.(var));
    a = sum(p(:).*mag(:)) / sum(mag(:));
end

function m = axial_massflow(c)
    [Ax,~,~] = cell_areas(c.x, c.r, c.rt);
    rovx = cell_corners(c.rovx);
    m = sum(rovx(:) .* Ax(:));
end

function v = streamtube_vx(c, mdot_target)
    % Area-averaged Vx over the capture streamtube: accumulate |axial mass flux|
    % radially hub->tip until it reaches |mdot_target|, then area-average Vx over
    % just that inner (contracted-jet) region (partial band crossing included).
    [Ax,Ar,At] = cell_areas(c.x, c.r, c.rt);
    mag = sqrt(Ax.^2 + Ar.^2 + At.^2);
    rovx = cell_corners(c.rovx); Vx = cell_corners(c.Vx);
    band_m   = abs(sum(rovx .* Ax, 2));                % |axial mass flux| per radial band
    band_A   = sum(mag, 2);                            % face area per band
    band_VxA = sum(Vx .* mag, 2);                      % Vx . area per band
    rc = mean(0.5*(c.r(1:end-1,:) + c.r(2:end,:)), 2); % radial centre of each band
    [~,ord] = sort(rc);
    band_m = band_m(ord); band_A = band_A(ord); band_VxA = band_VxA(ord);
    target = abs(mdot_target); cum = cumsum(band_m);
    if target >= cum(end); v = sum(band_VxA) / sum(band_A); return; end
    n = find(cum >= target, 1);
    prev = 0; if n > 1; prev = cum(n-1); end
    frac = 0; if band_m(n) > 0; frac = (target - prev) / band_m(n); end
    A_inc   = sum(band_A(1:n-1))   + frac * band_A(n);
    VxA_inc = sum(band_VxA(1:n-1)) + frac * band_VxA(n);
    v = VxA_inc / A_inc;
end

function [Ax,Ar,At] = cell_areas(x, r, rt)
    x = squeeze(x); r = squeeze(r); rt = squeeze(rt);
    t = rt ./ r;
    t_av = 0.25*(t(2:end,2:end)+t(1:end-1,1:end-1)+t(1:end-1,2:end)+t(2:end,1:end-1));
    x1 = x(2:end,2:end)-x(1:end-1,1:end-1); x2 = x(1:end-1,2:end)-x(2:end,1:end-1);
    r1 = r(2:end,2:end)-r(1:end-1,1:end-1); r2 = r(1:end-1,2:end)-r(2:end,1:end-1);
    rt1 = r(2:end,2:end).*(t(2:end,2:end)-t_av) - r(1:end-1,1:end-1).*(t(1:end-1,1:end-1)-t_av);
    rt2 = r(1:end-1,2:end).*(t(1:end-1,2:end)-t_av) - r(2:end,1:end-1).*(t(2:end,1:end-1)-t_av);
    Ax = 0.5*(r1.*rt2 - r2.*rt1);
    Ar = 0.5*(x2.*rt1 - x1.*rt2);
    At = 0.5*(x1.*r2  - x2.*r1);
end

function o = cell_corners(v)
    v = squeeze(v);
    o = 0.25*(v(1:end-1,1:end-1)+v(2:end,1:end-1)+v(1:end-1,2:end)+v(2:end,2:end));
end
