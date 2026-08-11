%   make_df_chic
%
%   Generate rotor- and stator-row inlet/outlet cut files for the ducted-fan
%   design points so that plot_df_loss.py can plot the spanwise loss without
%   MATLAB.
%
%   The ducted fan (bl_create_machine O4H topology) has two rows:
%     rotor  (rpm ~= 0): upstream H-block bid 2 (i=end = LE face),
%                        downstream H-block bid 0 (i=1 = TE face).
%     stator (rpm  = 0): upstream H-block bid 9 (i=end = LE face),
%                        downstream H-block bid 7 (i=1 = TE face).
%   All blocks span the full duct (hub -> casing), so no tip clipping is needed.
%
%   Each row's cut is saved with the row rpm (rotor speed for the rotor, 0 for the
%   stator) so that the loss is evaluated in the correct frame: relative
%   stagnation pressure for the rotor, absolute for the stator (Po_rel = Po when
%   rpm = 0).  Saved as <design>_rotor_chic.mat / <design>_stator_chic.mat with
%   the same {inlet,outlet,pout,rpm} struct as ts_cut_chic.m.

clear; close all; clc; eVTOL;

designs = {'Phi06DUCTEDFAN','Phi075DUCTEDFAN','Phi09DUCTEDFAN2'};

% row -> [upstream H-block (inlet), downstream H-block (outlet)]
rows = struct('name',{'rotor','stator'}, 'bid_in',{2,9}, 'bid_out',{0,7});

for d = 1:length(designs)
    des = designs{d};
    file = [dr.ts des '_avg.hdf5'];
    if isempty(dir(strrep(file,'.hdf5','.*'))); warning('missing %s',file); continue; end
    g = ts_read_hdf5(file);

    for k = 1:length(rows)
        bi = rows(k).bid_in; bo = rows(k).bid_out;

        % Row rpm (the moving blocks of this row); 0 for the stator
        rpm = g{bi+1}.bv.rpm;

        % Inlet at the LE face (upstream H-block i=end), outlet at the TE face
        % (downstream H-block i=1), full span (hub -> casing)
        inlet  = ts_secondary(ts_structured_cut(g,bi, g{bi+1}.attribute.ni, ...
                     g{bi+1}.attribute.ni, 1, 'en', 1, 'en'));
        outlet = ts_secondary(ts_structured_cut(g,bo, 1, 1, 1, 'en', 1, 'en'));
        pout = 0;

        save([dr.ts des '_' rows(k).name '_chic.mat'],'inlet','outlet','pout','rpm');

        % Quick overall loss check (relative for rotor, absolute for stator)
        Po_in  = ts_mass_average(inlet, 'Po_rel', 3);
        Po_out = ts_mass_average(outlet,'Po_rel', 3);
        P_in   = ts_area_average(inlet, 'P', 3);
        Yp = (Po_in - Po_out) / (Po_in - P_in);
        fprintf('%-16s %-6s  rpm=%5.0f  Yp_overall=%.4f\n', des, rows(k).name, rpm, Yp);
    end
end
