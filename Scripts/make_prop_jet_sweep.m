%   make_prop_jet_sweep
%
%   Sweep constant-x planes downstream of the propeller rotor so plot_jet_sweep.py
%   can show how the capture-streamtube jet velocity Vj evolves (contracts, then
%   mixes out) with distance, and a sensible jet station can be picked.
%
%   For each operating point (_O1/_O2/_O3, averaged solution) this sweeps every
%   near-axis, streamwise block downstream of the rotor, taking a constant-i
%   structured cut (full span, all pitch) at each i-station. A block qualifies if
%   it extends downstream of the rotor TE, reaches near the axis (so it contains
%   the contracted slipstream), and its i-index runs streamwise (axial). Index
%   cuts are robust on any mesh (no cross-block interpolation) and each is a clean
%   2D radial x pitch plane. Planes are pooled across blocks and sorted by x, then
%   saved into <file>_jetsweep.mat; the capture-streamtube average of Vx at each is
%   computed in Python (plot_jet_sweep.py) from the rotor mass flow in <file>_chic.mat.

clear; close all; clc; eVTOL;

designs = {'02062026PropPhi06HYBRIDrestart2', ...
           '02062026PropPhi06n05restart2', ...
           '11052026PropPhi06restart2'};
ops = {'_O1' '_O2' '_O3'};

bid_out = 0;        % downstream H-block: i=1 is the rotor TE face
r_near = 2.0;       % a block must reach within r_near*Rc of the axis to qualify
i_step = 2;         % subsample every i_step planes per block to keep files modest

for d = 1:length(designs)
    des = designs{d};
    base = regexprep(des,'restart.*$','');
    s = importdata([dr.geom base '-mission.json'],'r'); mm = jsondecode(s{1});
    Rc = mm.Rc;

    for o = 1:length(ops)
        file = [dr.ts des ops{o} '_avg.hdf5'];
        if isempty(dir(strrep(file,'.hdf5','.*'))); warning('missing %s',file); continue; end
        g = ts_read_hdf5(file);

        x_te = mean(mean(g{bid_out+1}.x(1,:,:)));

        P = {}; X = [];
        for bid = ts_get_bids(g)
            xb = g{bid+1}.x; rb = g{bid+1}.r;
            if max(xb(:)) <= x_te + 1e-9, continue; end              % must be downstream
            if min(rb(:)) > r_near * Rc, continue; end               % must be near the axis
            dxi = abs(xb(end,:,:) - xb(1,:,:)); dri = abs(rb(end,:,:) - rb(1,:,:));
            if mean(dxi(:)) < mean(dri(:)), continue; end            % skip if i is not streamwise
            ni = double(g{bid+1}.attribute.ni);
            for i = 1:i_step:ni
                xi = xb(i,:,:); xm = mean(xi(:));
                if xm < x_te - 1e-9, continue; end                   % only downstream planes
                P{end+1} = ts_secondary(ts_structured_cut(g, bid, i, i, 1, 'en', 1, 'en'));
                X(end+1) = xm;
            end
        end
        if isempty(X); warning('%s%s: no downstream planes found', des, ops{o}); continue; end

        [X, ord] = sort(X); planes = P(ord); xs = X;
        save(strrep(file,'.hdf5','_jetsweep.mat'), 'planes', 'xs', 'x_te', 'Rc');
        fprintf('%s%s  swept %d planes  x/Rc = %.2f .. %.2f\n', des, ops{o}, numel(xs), ...
            (xs(1)-x_te)/Rc, (xs(end)-x_te)/Rc);
    end
end
