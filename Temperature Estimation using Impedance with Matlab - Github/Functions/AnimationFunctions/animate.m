function animate(P,INTERP,RESULTS,PDEPE,hplot,ht)
% This function runs an animation of the battery temperature distirbution
% during the drive cycle, and saves it as a .gif file.
% 
% Copyright (c) 2016 by Robert Richardson, David Howey
% and The Chancellor, Masters and Scholars of the University of Oxford.
% See the licence file LICENCE.txt for more information.

% myVideo = VideoWriter('myfile.avi');
% myVideo.FrameRate=10;
% open(myVideo);

for k = 1:20:round(INTERP.t(P.animationTime))
    % Update core/surface temperature
    set(hplot(5) ,'xdata',[k k]);
    k
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
    
    set(ht, 'String', sprintf(['', num2str(INTERP.t(k)), ' s / ', num2str(INTERP.t(P.animationTime)), ' s']));
    drawnow
    % Get frame as an image
    f = getframe(gcf);
    
    im = frame2im(f);
    [A,map] = rgb2ind(im,256); 

	if k == 1;
        imwrite(A,map,'DEKF_animation.gif','gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,'DEKF_animation.gif','gif','WriteMode','append','DelayTime',0.1);
    end
    
end
end



