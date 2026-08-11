

des_name = ['FixRE_E0.31_phi0.75_t30_Bl4_DH0.95_rleanh10']; 
%g = ts_read_hdf5([dr.ts 'ASME_IPM5PROP_E0.31_phi0.75_t30_Bl4_DH0.95_rleanh10_3_O1norescalerestart_avg' '.hdf5']);

s = importdata([dr.geom des_name '-mission.json'],'r');
 mymission = jsondecode(s{1});

% Load blade geometry directly
xrrt = dlmread([dr.geom des_name '.csv'],',',1,0);

fileID = fopen('blade.ibl','w');
fwrite(fileID,xrrt );
fclose(fileID);


% Reshape and split into rotor and stator
B{1}.xrrt = reshape(xrrt(:,[1 3 2]),428,[],3);


if size(xrrt,2) == 7
    nn = 2;
    B{2}.xrrt = reshape(xrrt(:,[7 6 5]),428,[],3); %was 399
else
    nn = 1; 
end

% Initialise array to collect all surfaces and define CAD resolution
C = cell(0,1); nu = 81; nv = 21;

% Plot all blades
figure('units','normalized','position',[0.1 0.1 0.6 0.8],'color',[1 1 1]); 
axes('visible','off'); hold on; axis equal; view(-80,7); cols = lines(length(B));
for n = 1 %was originally 3
    
    % Downsample arrays
    i = 1:2:size(B{n}.xrrt,1); j = 1:1:size(B{n}.xrrt,2);
    u = round(linspace(1,size(B{n}.xrrt,1),nu)); v = round(linspace(1,size(B{n}.xrrt,2),nv));
    
    % Loop around the annulus
    for k = 1:mymission.NB_r 
        
        % Plot the surface
        dt = k * 2 * pi / mymission.NB_r + n * 0.2; 
        xrrt = m_shift(B{n}.xrrt,dt); 
        xyz = m_pol2cart(xrrt);
        %surf(xyz(i,j,1),-xyz(i,j,3),xyz(i,j,2),'facecolor',cols(n,:),'linestyle','none');
        
        % Record nurbs surface
        C{end+1,1} = nrbloft(1000 * xyz(u,v,1),1000 * xyz(u,v,2),1000 * xyz(u,v,3),2);
    end
end

% Plot cowl and hub surfaces
%nk = 1081; varnames = fieldnames(m);
%for o = 1:length(varnames)
    
    % Plot the surface
%    ni = size(m.(varnames{o}),1); xrrt = zeros(ni,nk,3); 
%    xrrt(:,:,1:2) = repmat(reshape(m.(varnames{o}),[ni 1 2]),[1 nk 1]);
%    t = linspace(0,2*pi,nk); xrrt(:,:,3) = xrrt(:,:,2) .* t; xyz = m_pol2cart(xrrt);
%    surf(xyz(:,:,1),-xyz(:,:,3),xyz(:,:,2),'facecolor',[0.5 0.5 0.5],'linestyle','none');
%    
%    % Record nurbs surface
%    u = round(linspace(1,size(xyz,1),nu)); v = round(linspace(1,size(xyz,2),nv));
%    C{end+1,1} = nrbloft(1000 * xyz(u,v,1),1000 * xyz(u,v,2),1000 * xyz(u,v,3),2);
%end

% Add lighting and zoom
axis tight; camlight;

% Save geometry to CAD readable file
igesout(C,'aero_surfaces');