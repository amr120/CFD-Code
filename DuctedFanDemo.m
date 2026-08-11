%   DuctedFanDemo
%                               
%   Demo the nacelle creation script
%
%   Written by James Taylor            
%   June 2024

% Clear workspace and load eVTOL environment
clear; close all; clc; eVTOL;

% Define the design to load
des_name = 'IPM5_V1I';

% Load blades from output of inverse design tool
load([dr.geom des_name '_prelim.mat']); 

% Define overall nacelle parameters
c.r_lip = 0.085; c.t_te = 1e-3; c.x_lip = -0.02;
c.thick_max = 0.02; c.chi_le = 5; 

% Define detailed nacelle section parameters
c.qcam = 0; c.dcam_le = -2; c.dcam_te = -0.5; c.qcam = 0;
c.wedge_te = 55; c.rad_le = 5; c.rad_thick_max = 0.12; c.s_thick_max = 0.25; 

% Define nozzle parameters
c.l_noz = 0.1; c.sigma = 1.4;

% Define centrebody parameters
c.e_nose = 1.5; c.l_tail = 0.05; c.drdx_tail = -0.4;

% Modify the gradients to test the blends
B{1}.drdx_hub_in = 0.1; B{1}.drdx_cas_in = 0; 
B{end}.drdx_hub_out = -0.1; B{end}.drdx_cas_out = 0; 

% Create the nacelle
c = bl_construct_nacelle(c,B,1);
