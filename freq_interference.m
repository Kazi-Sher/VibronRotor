% Copyright (C) 2018 Kazi Sher Ahmed <kazisherahmed@gmail.com> & SM Ahmad.
% This file is a part of VibronRotor - A finite-element code for rotordynamic analysis of flexible rotor-bearing systems.
% VibronRotor is released under the terms of GNU General Public License 3.0.

% Campbell Diagram
function freq_interference(num_modes, cd_range, increments, mbb, kbb, cbb, gbb)

    clear A evec1 evalu evalud evaludhz firstcolumns evalorder indexhz evalhzr;

    b=1;
    I = eye(length(mbb));      % identity matrix of (4n+4)*(4n+4) order or (4m X 4m)
    O = zeros(length(mbb));    % null matrix of (4n+4)*(4n+4) order 
 
    for arpm = 0:increments:cd_range
        arad = (arpm*pi)/30;   % rotational velocity of rotor in rad/s

        a(1,b) = arpm;

        A = [ O I ; -mbb\kbb -mbb\(cbb+(arad*gbb)) ];
        [evec1,evalu] = eig(A);    
        evalud = diag(evalu);
        evaludhz = evalud/(2*pi);    

        firstcolumns = 1:2:(length(A));
        evaludhz(firstcolumns) = []; % 4mX1
        [evalorder,indexhz] = sort(abs((evaludhz)));

        for cnt = 1:length(evaludhz)
            evalhzr(cnt,b) = (evaludhz(indexhz(cnt)));
        end
        damped_freq (:,b) = -imag (evalhzr(1:num_modes*2,b));  
        b = b+1;
    end
    figure
    for i = 1:2:num_modes*2
        c_modes_b = plot(a, damped_freq(i,:)*60, ('--'),'LineWidth',1); hold on;
        c_modes_f = plot(a, damped_freq(i+1,:)*60, ('-'),'LineWidth',1); hold on;
    end
     
    line_1x = plot (a,a,('k-.'),'LineWidth',0.6); hold off;
%    line_1x.Color = [0 0 0];

%     str = '1X';
%     tline = text ( 700,1000+400,str );
%     tline.Color = [0 0 0]; hold off;
    
    title(['Campbell Diagram']);
    xlabel('Rotor Speed (RPM)')
    ylabel('Damped Natural Frequency (CPM)')
%     axis([ 0 ana_range 0 14000 ])
    
    grid off
    set(gca,'box','off');    
    set(gcf,'color','w');
    set(gca,'fontsize',9)
    disp('Press Enter to continue to the next selected functionality.'); pause
end