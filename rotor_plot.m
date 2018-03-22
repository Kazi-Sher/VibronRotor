% Copyright (C) 2018 Kazi Sher Ahmed <kazisherahmed@gmail.com> & SM Ahmad.
% This file is a part of VibronRotor - A finite-element code for rotordynamic analysis of flexible rotor-bearing systems.
% VibronRotor is released under the terms of GNU General Public License 3.0.

function rotor_plot (l_segments, dist_segments, d_segments, num_discs, l_discs, d_discs, dist_discs, num_bearings, dist_bearings)

    figure
    title('Rotor Schematic');
    for a=1:length(l_segments)
        rectangle('Position',[dist_segments(a)  0-d_segments(a)/2 l_segments(a) d_segments(a)],'edgecolor',[223/255 0/255 0/255],'LineWidth',0.8); % shaft        
        str = num2str(a);
        tline = text ( dist_segments(a)+ 0.01, 0-d_segments(a)/2 + d_segments(a) + 0.02 , str );
%        tline.Color = [223/255 0/255 0/255];
    end
    hold on;
    
    if num_discs~=0
        for a=1:length(l_discs)
            rectangle('Position',[ dist_discs(a)-l_discs(a)/2  0-d_discs(a)/2  l_discs(a)  d_discs(a) ],'edgecolor',[0 0 159/255],'LineWidth',0.8,'LineStyle','--'); % first disc
            str = ['d' num2str(a) ''];
            tline = text ( dist_discs(a)-l_discs(a)/2 + 0.01, 0-d_discs(a)/2 + d_discs(a) + 0.02 , str );
%            tline.Color = [0 0 159/255];
        end
    end
    hold on;
    
    modeshape_axislength = [-1.1*max(d_discs)/2 0];
    for a=1:num_bearings
        dist_bearings_plot = [ dist_bearings(a) dist_bearings(a)];  
        plot( dist_bearings_plot, modeshape_axislength, ('m-.'));
    end
    hold off;      
    
%     rectangle('Position',[0 0-0.02616/2 0.0116 0.02616],'edgecolor','m','LineWidth',0.8,'facecolor','w'); % coupling
    
    axis equal;
    set(gca,'box','off');    
    set(gcf,'color','w');
    set(gca,'fontsize',8)
    xlabel('Rotor Length (m)')
    ylabel('Rotor Radius (m)')
    
    grid minor;
%     set(gca,'MinorGridLineStyle','-.');
    grid on;
    set(gca,'GridLineStyle','-');
    axis ([ 0 dist_segments(end) -1.1*max(d_discs)/2 1.1*max(d_discs)/2 ]);
    set(gca,'YMinorTick','off');  
%     set(gca,'XColor','r');
%     axis([0 lbeam -0.06 0.06]);    
        
    disp('Press Enter to continue to the next selected functionality.'); pause
end