close all; clc; clear;

eVTOL;

filename = '0.2epsProp_O1_restart_avg.hdf5';
g = ts_read_hdf5([dr.ts filename]);

% Open figure
figure(); hold on; grid on; box on;

% Plot static pressure distribution radiall
bids = [9 0 10]; cols = lines(length(bids)); o = 1;
i_plot = [5 33 ; 1 13 ; 1 47];
for bid = bids
    for i = 1:2
        
        % Take structured cut
        c = ts_structured_cut(g,bid,i_plot(o,i),i_plot(o,i),1,'en',1,'en');

        % Area average pressure in pitchwise direction
        P = ts_area_average(c,'P',2);
        
        % Extract radius
        r = g{bid+1}.r(i_plot(o,i),:,1);

        % Plot
        plot(P,r,'color',cols(o,:));

    end
    o = o+1;
end

% Add the hbu and tip radius
axis([100400 101500 0 1.2]);
plot([0 200000],[0.325 0.325],'k-')
plot([0 200000],[0.131 0.131],'k-')

% Add the atmospheric pressure
plot([101350 101350],[0 2],'k-')

% Export static pressure for paraview
ts_export_paraview(g,[dr.pv strrep(filename,'.hdf5','_pstat.hdf5')],...
    {'M_rel', 'sfunc', 'Alpha_rel', 'Po_rel', 'Po','P','mwall'},[],1)