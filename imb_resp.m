function imb_resp(orb_wrpm, fx_node_out, fx_node_in, ini_phase, ini_imb, imb_interval, ana_range, node_discs, kbb, mbb, cbb, gbb)

    % Imbalance Response Code [Equations: Chen pg321]

    % Max. allowable residual unbalance. API 617 standard [Vance 160]
    % rotor_mass = (pi*rbeam^2*overall_length*density + md_1 + md_2) * 2.20462 ; % rotor mass in lbm 
    % max_speed = ana_range; % max operating speed in rpm
    % res_imb = (4*rotor_mass/max_speed) * 0.0007200778893234; % in kg m

    % ecc = 0.04123; % in m, for 10mm imb disc, it is 41.23mm for outermost holes
    % m_imb = 0.005; % in kg

    imb_amount = ini_imb .* (1/(10^3)) * (1/39.3701); % a row vect containing imbalance amount. Also converts from g-in to kg-m
    ini_phase  = ini_phase.*(pi/180); % conversion from deg to rad
    
    IR = zeros(2*length(mbb));
    Force = zeros((2*length(mbb)),1);
    Imb_out = zeros((2*length(mbb)),1);

    c=1; % column counter
    for a = 0:imb_interval:ana_range % in RPM
        w_imb = (a*pi)/30;   % rotational velocity of rotor in rad/s
        RPM(1,c) = w_imb; 

        f_imb = imb_amount .* (w_imb^2); % in N

            for z=1:length(ini_imb) % [Chen321]
            % real component of force both in x and y directions
            Force(4*fx_node_in(z)-3,1) = f_imb(z)*cos(ini_phase(z));
            Force(4*fx_node_in(z)-2,1) = f_imb(z)*sin(ini_phase(z));

            % imag component of force both in x and y directions
            Force((4*fx_node_in(z)-3)+(4*length(mbb)/4),1) = -f_imb(z)*sin(ini_phase(z));
            Force((4*fx_node_in(z)-2)+(4*length(mbb)/4),1) = f_imb(z)*cos(ini_phase(z));
            end

        % Bearing Reactions
        % Force(4*node_bearings(1)-2,1) =  -116.31*0.453592*9.81 *sin(-90*2*pi/360); 
        % Force((4*node_bearings(1)-2)+(4*mnu),1) =  -116.31*0.453592*9.81 *cos(-90*2*pi/360); 
        % 
        % Force(4*node_bearings(2)-2,1) =  -222.23*0.453592*9.81 *sin(-90*2*pi/360); 
        % Force((4*node_bearings(2)-2)+(4*mnu),1) =  -222.23*0.453592*9.81 *cos(-90*2*pi/360); 

        % Force(2,1) = (density*lbeam*pi*(rbeam^2)*9.81)/2;
        % Force(2+4*mnu,1) = (density*lbeam*pi*(rbeam^2)*9.81)/2;
        % Force(end-2-4*mnu,1) = (density*lbeam*pi*(rbeam^2)*9.81)/2;
        % Force(end-2,1) = (density*lbeam*pi*(rbeam^2)*9.81)/2;

        IR = [ kbb-(mbb*(w_imb^2)) (cbb+gbb*w_imb)*w_imb ; -(cbb+gbb*w_imb)*w_imb (kbb-(mbb*(w_imb^2))) ];
        %IR = [ kbb-(mbb*(w_imb^2)) (cbb+gbb)*w_imb ; -(cbb+gbb)*w_imb (kbb-(mbb*(w_imb^2))) ];
        %Imb_out(:,a) = IR\Force;
        Imb_out(:,c) = IR\Force;
        c=c+1;

    end

    RPM = RPM*30/pi; % from rad/s into RPM

    for a= 1: length(RPM)
        Imb_out_magx (:,a) = sqrt ( ((Imb_out(4*fx_node_out-3,a) )^2) + ((Imb_out(4*fx_node_out-3+(4*length(mbb)/4),a))^2));
        % Imb_out_phasx (:,a) = atan2d ( -(Imb_out(4*fx_node_out-3,a)) , (Imb_out((4*fx_node_out-3)+(4*mnu),a)) );

        temp_phasex = (atan2d ( (Imb_out((4*fx_node_out-3)+(4*length(mbb)/4),a)) , (Imb_out(4*fx_node_out-3,a)) ) );
            if (temp_phasex<0)
                temp_phasex = temp_phasex + 360;
            end    
        Imb_out_phasx (:,a) = temp_phasex;

        Imb_out_magy (:,a) = sqrt ( ((Imb_out(4*fx_node_out-2,a) )^2) + ((Imb_out((4*fx_node_out-2)+(4*length(mbb)/4),a))^2));

        temp_phasey = (atan2d ( (Imb_out((4*fx_node_out-2)+(4*length(mbb)/4),a)) , (Imb_out(4*fx_node_out-2,a)) ) );
            if (temp_phasey<0)
                temp_phasey = temp_phasey + 360;
            end    
        Imb_out_phasy (:,a) = temp_phasey;

    %     Imb_out_phasy (:,a) = atan2d ( -(Imb_out(4*fx_node_out-2,a)) , (Imb_out((4*fx_node_out-2)+(4*mnu),a)) );
    %     Imb_out_phasy (:,a) = ( atan2d ( (Imb_out((4*fx_node_out-2)+(4*mnu),a)) , (Imb_out(4*fx_node_out-2,a)) ) );
    end
    % Not putting the minus sign with phase will lead to inverted plot about 0 degree line.
    figure
    subplot(2,1,1)
    FE_phasex = plot ( RPM , Imb_out_phasx(1,:) , ('-'),'LineWidth',1); hold on;
    FE_phasey = plot ( RPM , Imb_out_phasy(1,:) , ('--'),'LineWidth',1); hold on;
    legend([FE_phasex FE_phasey], {['X-Phase'],['Y-Phase']}, 'FontSize',10,'Location','northwest' );
    title(['Response at Node: ', num2str(fx_node_out)]);
    xlabel('Rotor Speed (RPM)');
    ylabel('Phase (Degrees)');

    grid off;
    set(gca,'box','off');    
    set(gcf,'color','w');
    set(gca,'fontsize',9);   

    subplot(2,1,2) 
    FE_ampx = plot ( RPM , Imb_out_magx(1,:)*1000, ('-'),'LineWidth',1); hold on; % if need mils, *39370.1 
    FE_ampy = plot ( RPM , Imb_out_magy(1,:)*1000 , ('--'),'LineWidth',1); hold on;

    legend([FE_ampx FE_ampy], {['X-Amplitude'], ['Y-Amplitude']}, 'FontSize',10,'Location','northwest' );

    xlabel('Rotor Speed (RPM)');
    ylabel('Amplitude (mm)');
    grid off;
    set(gca,'box','off');    
    set(gcf,'color','w');
    set(gca,'fontsize',9);
##    export_fig imb_balance.png;
%     disp('Enter to continue to the orbit plot'); pause

    %% Orbit Plotting

        % the more away rom is from critical speed, the better resolution (lesser time intervals) is needed to accurate orbit
        orb_wrad = (orb_wrpm*pi)/30;
        t_col = 1;
        imb_out_target = orb_wrpm/imb_interval + 1;

        for t = 0:0.0001:0.02 % the less the interval, the better resolution in orbit

        x_time (1,t_col) =  ( (Imb_out(4*fx_node_out-3,imb_out_target)*cos(orb_wrad*t)) + (Imb_out((4*fx_node_out-3)+(4*length(mbb)/4),imb_out_target)*sin(orb_wrad*t)) );
        y_time (1,t_col) =  ( (Imb_out(4*fx_node_out-2,imb_out_target)*cos(orb_wrad*t)) + (Imb_out((4*fx_node_out-2)+(4*length(mbb)/4),imb_out_target)*sin(orb_wrad*t)) );
        % pg 61 Chen

        t_col = t_col + 1; % time columns in integers
        end
            figure
            subplot(1,1,1)
            labels = cellstr( num2str([1:t_col-1]') ); % numbering the points
            labels = labels';

            FE_code_plus_sign = plot ( x_time(1)*1000, y_time(1)*1000, ('w+') ); hold on;
            plot ( x_time(1)*1000, y_time(1)*1000, ('k+') ); hold on;
            x_time_plot = plot ( x_time*1000, y_time*1000, ('-'),'LineWidth',1 ); hold on;
    %         text ( x_time*1000, y_time*1000, labels, 'VerticalAlignment','bottom', 'HorizontalAlignment','right');
           
            legend([ FE_code_plus_sign ], {['Orbit starts at +']}, 'FontSize',10,'Location','northwest' );

            title(['Orbit Plot at Node: ', num2str(fx_node_out), '   RPM: ', num2str(orb_wrpm)]);

            xlabel('X axis (mm)');
            ylabel('Y axis (mm)');
%             axis([-0.01 0.01 -0.01 0.01]);
            grid on;
            axis equal;

        set(gca,'box','off');    
        set(gcf,'color','w');
        set(gca,'fontsize',9);  
##        export_fig orbit.png;
%         disp('Press Enter to continue to the next selected functionality.'); pause
        
    %%
    % % Major / Minor axis length calculations
    % xr = Imb_out(4*fx_node-3,orb_wrpm);
    % xi = Imb_out(4*fx_node-3+(4*mnu),orb_wrpm);
    % yr = Imb_out(4*fx_node-2,orb_wrpm);
    % yi = Imb_out(4*fx_node-2+(4*mnu),orb_wrpm);
    % 
    % minor_axis = sqrt ( (((xr^2)+(xi^2)+(yr^2)+(yi^2))/2) - sqrt ( ((((xr^2)-(xi^2)+(yr^2)-(yi^2))^2)/4) + (xr*xi + yr*yi)^2 ) );
    % major_axis = sqrt ( (((xr^2)+(xi^2)+(yr^2)+(yi^2))/2) + sqrt ( ((((xr^2)-(xi^2)+(yr^2)-(yi^2))^2)/4) + (xr*xi + yr*yi)^2 ) );
end