%% With Cubic Spline: Instability Threshold ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function instab_threshold(th_modes, th_increments, th_range, support_coeff_speed, kxx_speed, kyy_speed,kxy_speed,kyx_speed, cxx_speed, cxy_speed, cyy_speed, cyx_speed, k, c_mat, num_bearings, node_bearings, mbb, gbb )

    % clear A evec1 evalu evalud evaludhz firstcolumns evalorder indexhz evalhzr;
    
    b = 1;
    I = eye(length(mbb));      % identity matrix of (4n+4)*(4n+4) order or (4m X 4m)
    O = zeros(length(mbb));    % null matrix of (4n+4)*(4n+4) order 

    for arpm = 0:th_increments:th_range
        arad = (arpm*pi)/30;   % rotational velocity of rotor in rad/s
        a(1,b) = arpm;

        %Cubic spline curve fitting
        kxx_fit = spline(support_coeff_speed,kxx_speed,arpm);
        kyy_fit = spline(support_coeff_speed,kyy_speed,arpm);
        kxy_fit = spline(support_coeff_speed,kxy_speed,arpm);
        kyx_fit = spline(support_coeff_speed,kyx_speed,arpm);

        cxx_fit = spline(support_coeff_speed,cxx_speed,arpm);
        cyy_fit = spline(support_coeff_speed,cyy_speed,arpm);
        cxy_fit = spline(support_coeff_speed,cxy_speed,arpm);
        cyx_fit = spline(support_coeff_speed,cyx_speed,arpm);

        clear k_temp c_temp;
        k_temp = k;
        c_temp = c_mat;

        for num_b = 1:num_bearings
            % Bearing nodes have already been located. Bearing stiffness addition in global stiffness matrix. Order below is xx,yy,xy,yx.
            k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) = k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) + kxx_fit;
            k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) = k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) + kyy_fit;
            k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) = k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) + kxy_fit;
            k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) = k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) + kyx_fit;

            c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) = c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) + cxx_fit;
            c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) = c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) + cyy_fit;
            c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) = c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) + cxy_fit;
            c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) = c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) + cyx_fit;
        end

        A = [ O I ; (-mbb\k_temp) (-mbb\(c_temp+arad*gbb)) ];
        [evec1,evalu] = eig(A);    
        evalud = diag(evalu);
        evaludhz = evalud;    

        firstcolumns = 1:2:(length(A));
        evaludhz(firstcolumns) = []; % 4mX1
        [evalorder,indexhz] = sort(abs((evaludhz)));

        for cnt = 1:length(evaludhz)
            evalhzr(cnt,b) = (evaludhz(indexhz(cnt)));
        end
        text
        damped_freq (:,b) = real (evalhzr(1:th_modes*2,b));  
        b = b+1;
    end
    z=1; 
    figure
    for i = 1:2:th_modes*2
        plot(a, damped_freq(i,:), ('--'),'LineWidth',1,'DisplayName', ['' num2str(z) 'b']); hold on;
        plot(a, damped_freq(i+1,:), ('-'),'LineWidth',1,'DisplayName', ['' num2str(z) 'f']); hold on;
        z=z+1;
    end 

    legend('show','Location','southwest');
    title('Instability Threshold Map');
    xlabel('Rotor Speed (RPM)')
    ylabel('Damping Exponent (rad/s)')
    grid on;
    set(gca,'box','off');    
    set(gcf,'color','w');
    set(gca,'fontsize',9);  
##    export_fig instab.png;
%     disp('Press Enter for the next selected functionality.'); pause
end