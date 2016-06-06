function updatePlot(hObject,hplot,h2,P,INTERP,RESULTS,PDEPE)
% This function upodates the data on the graph depending on the time step
% considered (given by the value 'k' of the slider)
% 
% Copyright (c) 2016 by Robert Richardson, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

    k = round(get(hObject,'Value'));    % Ensure 'k' is integer   
    set(h2,'String',['Time step: ' num2str(INTERP.t(k)) ' s' ' / ' num2str(INTERP.t(P.sliderTime)) ' s']);

    
    % Update core/surface temperature
    set(hplot(5) ,'xdata',[k k]);
    
    % Update convection coefficient
    set(hplot(8) ,'xdata',[k k]);
    
    % Update voltage
    set(hplot(10),'xdata',[k k]);
    
    % Update current
    set(hplot(12),'xdata',[k k]);
    
    % Update impedance
    set(hplot(14),'xdata',[k k]);
    
    % Update T(r,t)
    set(hplot(16),'xdata',P.r_sym(:,k));
    set(hplot(16),'ydata',RESULTS.DEKF.T_rt_sym(:,k));
    set(hplot(17),'xdata',PDEPE.r_sym(:,k));
    set(hplot(17),'ydata',PDEPE.T_rt_sym(:,k));
    
    % Update Q(r,t)
    set(hplot(20),'xdata',P.r_sym(:,k));
    set(hplot(20),'ydata',RESULTS.Q_rt_sym(:,k));
    
    drawnow;