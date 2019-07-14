%% Eigenanalysis
function [evaludhz, evalhzr] = mode_shapes(num_plot, lvec, w, mbb, kbb, cbb, gbb, num_bearings, dist_bearings, l_segments, dist_segments, d_segments, num_discs, l_discs, dist_discs, d_discs, lbeam)

    I = eye(length(mbb));      % identity matrix of (4n+4)*(4n+4) order or (4m X 4m)
    O = zeros(length(mbb));    % null matrix of (4n+4)*(4n+4) order

    A = [ O I ; -mbb\kbb -mbb\(cbb+w*gbb) ] ;

    % MATLAB help: [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V 
    % whose columns are the corresponding right eigenvectors, so that A*V = V*D
    % Each column in evec1 corresponds to a mode of vibration. Each column has
    % relative magnitudes of all DOFs in the system. [Inman 4th ed pg. 350]
    % Order of elements in displacement matrix gives info about which 
    % eigenvector relates to which DOF.
    [evec1,evalu] = eig(A);

    % Ideally, the eigenvalue decomposition satisfies the relationship. 
    % Since eig performs the decomposition using floating-point computations, 
    % then A*V can, at best, approach V*D. In other words, A*V - V*D is close 
    % to, but not exactly, 0.
    % check = A*evec1 - evec1*evalu;

    evalud = diag(evalu); %square matrix containing eigenvalues converted to a 
    % vector of eigenvalues, order 8mX1
    evaludhz0 = evalud/(2*pi); % Conversion into Hz (order: 8mX1)

    evec2 = evec1( (1:((length(evec1))/2)) , 1:length(evec1) ); 
    % Lower part contains mode shapes multiplied by lambda.
    % [Inman 4th ed. pg. 403,402] Taking upper half. evec2 is 4mX8m.

    firstcolumns = 1:2:(length(A));
    evaludhz = evaludhz0;
    evaludhz(firstcolumns) = []; % 4mX1
    evec2(:,firstcolumns) = []; % 4mX4m

    % now reorder the eigenvalues and eigenvectors from low to high freq
    [evalorder,indexhz] = sort(abs((evaludhz))); % abs gives magnitude of 
    % complex number. Indexhz is vector of indexes pointing to unsorted
    % locations.

    for cnt = 1:length(evaludhz)
        evalhzr(cnt,1) = (evaludhz(indexhz(cnt)));
    end

    zeta = -(real(evalhzr(:,1))) ./ abs(evalhzr(:,1));
    c_factor = ((1-(zeta.^2)).*(1-(2*(zeta.^2)))).^0.5; % to account for difference b/w damped nat freq & critical speeds as per Vance 2010

    for cnt = 1:length(evec2) % loops runs for all columns
        evec(:,cnt) = evec2(:,indexhz(cnt)); % Unsorted columns of eigenvectors 
        % are linked to corresponding unsorted eigenvalue rows ( since all columns 
        % in original 'evalu' matrix have been merged into a vector. So we can't 
        % have a separate sorting for eigenvectors. evec is 4mX4m.
    end

     % Extraction of displacements for plotting

     % Raw matrices
     % evec_disp_mode_x_raw = zeros(length(evec2)/8,length(evec2));
     % evec_disp_mode_y_raw = zeros(length(evec2)/8,length(evec2));

     evec_disp_mode_absx = zeros(length(evec)/4,length(evec)); % new matrix to 
     % contain nodal translations for plotting, matrix of 
     % order mX4m with each column corresponding to a mode shape
     % Only one plane translations.      
     evec_disp_mode_absy = zeros(length(evec)/4,length(evec));
     evec_disp_mode_phasex = zeros(length(evec)/4,length(evec));
     evec_disp_mode_phasey = zeros(length(evec)/4,length(evec));

     % Extraction of translations from eigenvectors in evec. Only real part
     % of x and y translations are taken in evec_disp_mode_x & 
     % evec_disp_mode_x respectively, i.e. first, fifth, ninth elements and
     % so on for each column for x translations. evec_disp_mode is filled 
     % row-wise (i.e. first a complete row is filled in and then second 
     % and so on. No. of modes remain same but we took out eigenvectors of our interest only.
     for cnt = 1:length(evec)/(4)
         for cnt1 = 1:length(evec)

             evec_disp_mode_absx(cnt,cnt1) = abs ( evec( (4*cnt-3) , cnt1));
             evec_disp_mode_absy(cnt,cnt1) = abs ( evec( (4*cnt-2) , cnt1));

             evec_disp_mode_phasex(cnt,cnt1) = angle ( evec( (4*cnt-3) , cnt1));
             evec_disp_mode_phasey(cnt,cnt1) = angle ( evec( (4*cnt-2) , cnt1));

             evec_disp_mode_xr(cnt,cnt1) = evec_disp_mode_absx(cnt,cnt1) .* cos (evec_disp_mode_phasex(cnt,cnt1));
             evec_disp_mode_xi(cnt,cnt1) = evec_disp_mode_absx(cnt,cnt1) .* sin (evec_disp_mode_phasex(cnt,cnt1));

             evec_disp_mode_yr(cnt,cnt1) = evec_disp_mode_absy(cnt,cnt1) .* cos (evec_disp_mode_phasey(cnt,cnt1));
             evec_disp_mode_yi(cnt,cnt1) = evec_disp_mode_absy(cnt,cnt1) .* sin (evec_disp_mode_phasey(cnt,cnt1)); 

         end
     end

%% Separate Repeated Plotting
    modeshape_axislength = [-1.3 1.6];
    
    for mode_cnt = 1:num_plot

        % Backward Whirl ======================================================
        evec_cnt = 2*mode_cnt-1;

        max_xr = max(abs(evec_disp_mode_xr(:,evec_cnt)));
        max_xi = max(abs(evec_disp_mode_xi(:,evec_cnt)));
        max_xm = max(max_xr,max_xi);

        max_yr = max(abs(evec_disp_mode_yr(:,evec_cnt)));
        max_yi = max(abs(evec_disp_mode_yi(:,evec_cnt)));
        max_ym = max(max_yr,max_yi);

        figure
        subplot(1,2,1) %+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        for a=1:num_bearings
            dist_bearings_plot = [ dist_bearings(a) dist_bearings(a)];  
            plot( dist_bearings_plot, modeshape_axislength, ('m--'));
            hold on;
        end

        xr_axis = plot(lvec,(evec_disp_mode_xr(:,evec_cnt) / max_xm),('r-')); hold on;
        xi_axis = plot(lvec,(evec_disp_mode_xi(:,evec_cnt) / max_xm),('r--')); hold on;        
        yr_axis = plot(lvec,(evec_disp_mode_yr(:,evec_cnt) / max_ym ),('b-')); hold on;
        yi_axis = plot(lvec,(evec_disp_mode_yi(:,evec_cnt) / max_ym),('b--')); hold on;
        legend ([xr_axis xi_axis yr_axis yi_axis], {['X-axis real'], ['X-axis imag.'], ['Y-axis real'], ['Y-axis imag.']}, 'FontSize',7,'Location','northwest' );
        % lvec is 1Xm.  
        title(['Backward Whirl: Mode ', ...
        num2str(mode_cnt),': ',num2str(-imag(evalhzr(evec_cnt))),' Hz, ',num2str((-imag(evalhzr(evec_cnt)))*60),' RPM  Rotation: ',num2str(w),' RPM']);
        %num2str(mode_cnt),': ',num2str(-imag(evalhzr(evec_cnt))./c_factor(evec_cnt)),' Hz, ',num2str((-imag(evalhzr(evec_cnt))./c_factor(evec_cnt))*60),' RPM  Rotation: ',num2str(w1),' RPM']);    
        % Since eig values appear in conjugate pairs, we have pairs of nat freq.
        % Each mode has two nat freq. Seems like one for forward and one for
        % backward whirl. So first, third, fifth and so on rows of evalhzr and
        % first, third, fifth and so on columns of evec_disp_mode are being used.

        % Plotting the rotor for mode shapes
        for a=1:length(l_segments)
            rectangle('Position',[dist_segments(a)  0-d_segments(a)*4.347/2  l_segments(a)  d_segments(a)*4.347],'edgecolor','k','LineWidth',0.4); % shaft
        end

        if num_discs~=0;
            for a=1:length(l_discs)
                rectangle('Position',[ dist_discs(a)-l_discs(a)/2  0-d_discs(a)*4.347/2  l_discs(a)  d_discs(a)*4.347 ],'edgecolor','k','LineWidth',0.4); % first disc
            end
        end

        grid off;
        set(gca,'box','off');    
        set(gcf,'color','w');
        set(gca,'fontsize',10);
        axis([0 lbeam -1.3 1.6]);
        xlabel('Rotor Length [m]');
        ylabel('Normalized Displacement');


        % Forward Whirl =======================================================
        evec_cnt = 2*mode_cnt;

        max_xr = max(abs(evec_disp_mode_xr(:,evec_cnt)));
        max_xi = max(abs(evec_disp_mode_xi(:,evec_cnt)));
        max_xm = max(max_xr,max_xi);

        max_yr = max(abs(evec_disp_mode_yr(:,evec_cnt)));
        max_yi = max(abs(evec_disp_mode_yi(:,evec_cnt)));
        max_ym = max(max_yr,max_yi);

        subplot(1,2,2) %+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        for a=1:num_bearings
            dist_bearings_plot = [ dist_bearings(a) dist_bearings(a)];  
            plot( dist_bearings_plot, modeshape_axislength, ('m--'));
            hold on;
        end
        
        xr_axis = plot(lvec,(evec_disp_mode_xr(:,evec_cnt) / max_xm),('r-')); hold on;
        xi_axis = plot(lvec,(evec_disp_mode_xi(:,evec_cnt) / max_xm),('r--')); hold on;        
        yr_axis = plot(lvec,(evec_disp_mode_yr(:,evec_cnt) / max_ym ),('b-')); hold on;
        yi_axis = plot(lvec,(evec_disp_mode_yi(:,evec_cnt) / max_ym),('b--')); hold on;
        legend ([xr_axis xi_axis yr_axis yi_axis], {['X-axis real'], ['X-axis imag.'], ['Y-axis real'], ['Y-axis imag.']}, 'FontSize',7,'Location','northwest' );
        title(['Forward Whirl: Mode ', ...
        num2str(mode_cnt),': ',num2str(-imag(evalhzr(evec_cnt))),' Hz, ',num2str((-imag(evalhzr(evec_cnt)))*60),' RPM  Rotation: ',num2str(w),' RPM']);
        
        % Plotting the rotor for mode shapes
        for a=1:length(l_segments)
            rectangle('Position',[dist_segments(a)  0-d_segments(a)*4.347/2  l_segments(a)  d_segments(a)*4.347],'edgecolor','k','LineWidth',0.4); % shaft
        end

        if num_discs~=0;
            for a=1:length(l_discs)
                rectangle('Position',[ dist_discs(a)-l_discs(a)/2  0-d_discs(a)*4.347/2  l_discs(a)  d_discs(a)*4.347 ],'edgecolor','k','LineWidth',0.4); % first disc
            end
        end
        hold off;

        grid off;
        set(gca,'box','off');    
        set(gcf,'color','w');
        set(gca,'fontsize',10);
        axis([0 lbeam -1.3 1.6]);
        xlabel('Rotor Length [m]');
        ylabel('Normalized Displacement');
##    disp('Press Enter to display next mode'); pause
    end
##    disp('Press Enter to continue to the next selected functionality.'); pause
end