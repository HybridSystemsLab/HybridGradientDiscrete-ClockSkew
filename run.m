%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all
close all
clc

%% Parameters

% plant parameters
epsilon = 0.1; theta = 1 + epsilon;

A = [0 1 + epsilon; 0 0]; 
B = [0; 0] ; 
H = [0 0; 0 -1];
M = [0; 1];

% estimator parameters
s = 0.05;               % sample period during flows
gammac = 20;
gammad = 2;

parameters.A = A;
parameters.B = B;
parameters.H = H;
parameters.M = M;
parameters.s = s;
parameters.gammac = gammac;
parameters.gammad = gammad;


%% Create hybrid system
sys = ClkSk(parameters);

tspan = [0, 10];
jspan = [0, 2*tspan(end)/s];
config = HybridSolverConfig('RelTol', 1e-4, 'AbsTol', 1e-4,  ...
                            'MaxStep', s/4, 'Refine', 4);

z0 = [0; 1];             % plant state
thetahat0 = zeros(1,1); % parameter estimate
tau0 = 0;                % timer for samples during flows
%x0 = tau0;
%x0 = [z0;tau0;tau0];
x0 = [z0;tau0;tau0;thetahat0;thetahat0;thetahat0];
sol = sys.solve(x0, tspan, jspan, config);

%%{

%% Plots
HybridPlotBuilder.defaults.set('flow line width', 3, ...
                               'jump line width', 2,...
                               'label size', 34,...
                               'tick label size', 20,...
                               't_label', '$t$ [s]')
%legendSize = 20;

figure; clf;
tlt = tiledlayout(1, 1);

nexttile(1)
hold on; grid on;
HybridPlotBuilder()....
    .flowColor('blue')...
    .jumpMarker('none')...
    .labels('$\tau$')...
    .plotFlows(sol.select(1))
xlim(tspan)
set(gca,'Box','on');

tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), 1.3/2*pos(4)])
movegui(gcf,'north')

%%
inD1 = sol.x(:,4) >= 0.5;    % get points in the jump set
isFlowing = ~inD1;                               % get points during flows
isJumping = inD1 + cat(1,0,inD1(1:end-1));       % get points either side of each jump

figure; clf;
hold on; grid on;
HybridPlotBuilder().... % plot only the flows
    .flowColor('blue')...
    .jumpColor('none')...
    .labels('$|\tilde\theta_s|$')...
    .plotFlows(sol.select(5),@(x) norm(theta(:) - x))
HybridPlotBuilder().... % plot the jumps due to samples during flows
    .flowColor('none')...
    .jumpColor('blue')...
    .jumpLineStyle('-')...
    .jumpMarker('none')...
    .labels('$|\tilde\theta_s|$')...
    .filter(isFlowing)...
    .plotFlows(sol.select(5),@(x) norm(theta(:) - x))
HybridPlotBuilder().... % plot the jumps due to the jump set
    .flowColor('none')...
    .jumpColor('red')...
    .jumpMarker('none')...
    .labels('$|\tilde\theta_s|$')...
    .filter(isJumping)...
    .plotFlows(sol.select(5),@(x) norm(theta(:) - x))
xlim(tspan)
ylim([0 1.8])
set(gca,'Box','on');

xlim(tspan)
ylim([0 1.5])
set(gca,'Box','on');

tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), pos(4)])
movegui(gcf,'north')
set(gca, 'LooseInset', get(gca,'TightInset'))

%}
