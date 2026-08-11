%   make_prop_chic
%
%   Generate rotor inlet/outlet cut files for propeller cases so that
%   plot_prop_chic.py can plot the flow-coefficient / pressure-rise
%   characteristic without MATLAB.
%
%   For each operating point (_O1/_O2/_O3, averaged solution) this takes a cut
%   at the rotor leading-edge face (upstream H-block, bid 2 i=end) and the
%   trailing-edge face (downstream H-block, bid 0 i=1), restricted to the blade
%   span (hub -> tip at Rc). It also cuts the full-domain inlet and outlet planes
%   (freestream and jet). The {inlet,outlet,dom_in,dom_out,pout,rpm} struct is
%   saved into  <file>_chic.mat .  phi / Cptt / eta_poly (across the blade row),
%   and the velocity ratios J / Chi / eps (from dom_in/dom_out and the blade-row
%   planes) are then computed in Python exactly as in TS/post/ts_calc_chic.m.

clear; close all; clc; eVTOL;

% Designs to process and the operating-point suffixes
designs = {'02062026PropPhi06HYBRIDrestart2', ...
           '02062026PropPhi06n05restart2', ...
           '11052026PropPhi06restart2'};
ops = {'_O1' '_O2' '_O3'};

% Rotor inlet/outlet blocks (mesh_prop O4H topology)
bid_in = 2;    % upstream H-block - i=end is the LE face
bid_out = 0;   % downstream H-block - i=1 is the TE face

for d = 1:length(designs)
    des = designs{d};

    % Blade tip radius from the mission file (strip any 'restart...' suffix)
    base = regexprep(des,'restart.*$','');
    s = importdata([dr.geom base '-mission.json'],'r'); mm = jsondecode(s{1});
    Rc = mm.Rc;

    for o = 1:length(ops)
        file = [dr.ts des ops{o} '_avg.hdf5'];
        if isempty(dir(strrep(file,'.hdf5','.*'))); warning('missing %s',file); continue; end
        g = ts_read_hdf5(file);

        % Tip j-index on each face (closest grid line to Rc)
        [~,jt_in]  = min(abs(g{bid_in+1}.r(end,:,1) - Rc));
        [~,jt_out] = min(abs(g{bid_out+1}.r(1,:,1)   - Rc));

        % Machine speed = max-min rpm across all blocks (as ts_cut_chic.m)
        rpm_v = {'rpm','rpmi1','rpmi2','rpmj1','rpmj2','rpmk1','rpmk2'};
        rr = [];
        for bid = ts_get_bids(g)
            for v = 1:length(rpm_v); rr(end+1) = g{bid+1}.bv.(rpm_v{v}); end
        end
        rpm = max(rr) - min(rr);

        % Blade-span cuts at the LE and TE faces
        inlet  = ts_secondary(ts_structured_cut(g,bid_in, g{bid_in+1}.attribute.ni, ...
                     g{bid_in+1}.attribute.ni, 1, jt_in, 1, 'en'));
        outlet = ts_secondary(ts_structured_cut(g,bid_out, 1, 1, 1, jt_out, 1, 'en'));
        pout = 0;   % not used downstream, kept for ts_cut_chic compatibility

        % Full-domain inlet/outlet planes for the freestream and jet velocities
        % (cf. ts_cut_chic.m): kind==0 (or kind==13 freestream) is the domain
        % inlet, kind==1 is the domain outlet. Cut the whole plane (hub->far field).
        bid_din = []; bid_dout = [];
        for bid = ts_get_bids(g)
            for pid = ts_get_pids(g,bid)
                k = g{bid+1}.patch{pid+1}.attribute.kind;
                if (k == 0 || k == 13) && isempty(bid_din)
                    bid_din = bid; i_din = g{bid+1}.patch{pid+1}.attribute.ist;
                end
                if k == 1 && isempty(bid_dout)
                    bid_dout = bid; i_dout = g{bid+1}.patch{pid+1}.attribute.ist;
                end
            end
        end
        if isempty(bid_din) || isempty(bid_dout)
            warning('%s%s: could not find domain inlet/outlet patches', des, ops{o});
        end
        dom_in  = ts_secondary(ts_structured_cut(g,bid_din,  i_din+1,  i_din+1,  1, 'en', 1, 'en'));
        dom_out = ts_secondary(ts_structured_cut(g,bid_dout, i_dout+1, i_dout+1, 1, 'en', 1, 'en'));

        save(strrep(file,'.hdf5','_chic.mat'),'inlet','outlet','dom_in','dom_out','pout','rpm');
        c = ts_calc_chic(inlet,outlet,rpm);
        fprintf('%s%s  rpm=%5.0f  phi=%.4f  Cptt=%.4f  eta=%.4f\n', ...
            des, ops{o}, rpm, c.phi, c.Cptt, c.eta_poly);
    end
end
