%%
%-------------------------------------------------------------------------%
% Plot slider/create animation
%-------------------------------------------------------------------------%
%
%{
This code creates an animation comparing the temperature distribution 
predicted by the dual extended Kalman filter algorithm to experimental
measurements. Note that only the surface and core temperatures are
measured and so a PDE solution to the heat equation was fitted to the 
experimental data in order to provide the full experimental temperature 
distribution for illustration purposes.

A published html document showing this code and its results is
located in the 'html' sub-folder within this package.

Usage:
Select whether to run the slider tool (which allows the user to manually
step back and forth through the simulation) or to create a .gif image by
selecting the appropriate value for the 'selectAnimation' property below.
Then run the simulation as standard.


I would ask that you cite this paper as Richardson, Robert R., and 
David A. Howey. "Sensorless battery internal temperature estimation using a
kalman filter with impedance measurement." Sustainable Energy, IEEE 
Transactions on 6.4 (2015): 1190-1199. if you want to use this code for 
your own research. For further details on the work of the Energy Power 
Group at Oxford, please see epg.eng.ox.ac.uk.

Copyright (c) 2016 by Robert Richardson, David Howey
and The Chancellor, Masters and Scholars of the University of Oxford.
See the licence file LICENCE.txt for more information.


%}


%%
%-------------------------------------------------------------------------%
% Clear environment and required data
%-------------------------------------------------------------------------%

% Clear environment
clear;
close all;
clc;

addpath(genpath('./Functions'));
addpath(genpath('./Data'));

% Set plotting preferences
set(0,'defaultlinelinewidth',1.5)
set(0,'DefaultFigureColor','White');
set(0,'DefaultAxesBox','On');

% Load results
load('TemperatureDistributionData.mat')


%%
%-------------------------------------------------------------------------%
% Choose to create animation or run slider tool
%-------------------------------------------------------------------------%
% Choose whether to create animation or run slider tool

% selectAnimation = 1   -> Run slider tool;
% selectAnimation = 2   -> Create .gif animation (note: heavy memory req.);
selectAnimation = 1;

%%
%-------------------------------------------------------------------------%
% Plots
%-------------------------------------------------------------------------%

% Position and size of figure
x_fig = 350;
y_fig = 300;
w_fig = 900;
h_fig = 500;
hFig = figure;
set(hFig, 'Position', [x_fig y_fig w_fig h_fig])
k = 3450;

% Plot core/surface temperature
haxis(1)    = subplot(6,10,[1:4 11:14]);
hold on;
xlim([0 3500])
ylim([5 30])
hylab = ylabel('T / ^\circC');
set(hylab,'Position',(get(hylab,'Position') + [+135 0 0]))
set(gca, 'xticklabel', [])

hplot(1)    = plot(RAW.TEMP.t,  RAW.TEMP.T_core ,       'k-');
hplot(2)    = plot(RAW.TEMP.t,  RAW.TEMP.T_surf ,       'k--');
hplot(3)    = plot(INTERP.t,   RESULTS.DEKF.T_core,    'b-');
hplot(4)    = plot(INTERP.t,   RESULTS.DEKF.T_surf,    'b--');
hplot(5)    = plot([k k],       [0 100],                'r-','linewidth',0.5);


% Plot convection coefficient
haxis(2)    = subplot(6,10,21:24);
hold on;
xlim([0 3500])
ylim([20 100])      
hylab = ylabel('h / Wm^-^2k^-^1');
set(hylab,'Position',(get(hylab,'Position') + [35 0 0]))
set(gca, 'xticklabel', [])

hplot(6)    = plot(INTERP.t,   P.h_true*ones(size(INTERP.t)), 'k:');
hplot(7)    = plot(INTERP.t,   [RESULTS.DEKF.h_mat],           'b-');
hplot(8)    = plot([k k],       [0 100],                        'r-','linewidth',0.5);


% Plot voltage
haxis(3)    = subplot(6,10,31:34);
hold on
xlim([0 3500])
hylab = ylabel('V (V)');
set(hylab,'Position',(get(hylab,'Position') + [-4 2.2 0]))
set(gca, 'xticklabel', [])

hplot(9)    = plot(INTERP.t,       INTERP.V,    'k','linewidth',0.5);
hplot(10)   = plot([k k],           [2 4],      'r-','linewidth',0.5);


% Plot current
haxis(4)    = subplot(6,10,41:44);
hold on
xlim([0 3500])
ylim([-50 50])
hylab = ylabel('I (A)');
set(hylab,'Position',(get(hylab,'Position') + [5 0 0]))

set(gca, 'xticklabel', [])
hplot(11)   = plot(INTERP.t,   INTERP.I,        'k','linewidth',0.5);
hplot(12)   = plot([k k],       [-50 50],       'r-','linewidth',0.5);


% Real Impedance
haxis(5)    = subplot(6,10,51:54);
hold on
xlim([0 3500])
ylim([2e-3 5e-3])
hylab = ylabel('Z'' (\Omega)');
xlabel('t / s');
set(hylab,'Position',(get(hylab,'Position') + [-120 0 0]))

hplot(13)   = plot(RAW.EIS.t,   RAW.EIS.Re_Z,   'k','linewidth',0.5);
hplot(14)   = plot([k k],       [2e-3 5e-3],    'r-','linewidth',0.5);


% Plot T(r,t)
haxis(6)    = subplot(6,10,[6:10 16:20 26:30 36:40]);
hold on
xlim([-P.r_o P.r_o])
ylim([8 28])
ylabel 'T (\circC)'
set(gca, 'xticklabel', [])

hplot(15)   = plot([0 0],               [8 28],                    'k-.');
hplot(16)   = plot(P.r_sym(:,k),        RESULTS.DEKF.T_rt_sym(:,k), '-');
hplot(17)   = plot(PDEPE.r_sym(:,k),    PDEPE.T_rt_sym(:,k),        'k-');


% Plot heat generation
haxis(7)    = subplot(6,10,[46:50 56:60]);
hold on
xlim([-P.r_o P.r_o]);
ylim([-0.5 14]);
xlabel('r (m)');
ylabel('Q (W)');

hplot(18)   = plot([-P.r_o P.r_o],      [0 0],                  'k:');
hplot(19)   = plot([0 0],               [0 14],                 'k-.');
hplot(20)   = plot(P.r_sym(:,k),        RESULTS.Q_rt_sym(:,k),     'b-');

% Link axes
linkaxes([haxis(1) haxis(2) haxis(3) haxis(4) haxis(5)],'x');
linkaxes([haxis(6) haxis(7)],'x');
% linkaxes([haxis(6) haxis(1)],'y');

%%
%-------------------------------------------------------------------------%
% Implement slider
%-------------------------------------------------------------------------%

if selectAnimation == 1;
P.sliderTime = 3501;

h1 =    uicontrol(  'style','slider',...
                'units','pixel',...
                'position',[(w_fig-200)*0.025 10 200 20],...
                'Min',1,'Max',size(INTERP.t(1:P.sliderTime),1),'Value',1,...
                'SliderStep',[1/size(INTERP.t(1:P.sliderTime),1) 25/size(INTERP.t(1:P.sliderTime),1)]);
            
h2 =    uicontrol('Style','text',...
        'Position',[(w_fig-150)*0.975 10 150 20],...
        'String',['Time step: ' num2str(INTERP.t(k+1)) ' s' ' / ' num2str(INTERP.t(P.sliderTime)) ' s']);
 
    
% addlistener(h1,'ActionEvent',@(hObject, event) updatePlot(hObject,hplot,h2,result,xx));
addlistener(h1,'ContinuousValueChange',@(hObject, event) updatePlot(hObject,hplot,h2,P,INTERP,RESULTS,PDEPE));
end

%%
%-------------------------------------------------------------------------%
% Implement animation
%-------------------------------------------------------------------------%
if selectAnimation == 2;
P.animationTime = 3501;
axes(haxis(6))
ht = text(P.r_o*0.98,28*0.98,...
	['', num2str(INTERP.t(3200)), ' s / ', num2str(INTERP.t(P.animationTime)), ' s'],...
	'HorizontalAlignment','right','VerticalAlignment','top',...
        'fontsize',9,...
	'BackgroundColor',[.85 .85 .85]);

% Run and save animation
animate(P,INTERP,RESULTS,PDEPE,hplot,ht)
end


