function [node_discs, seg_dia_repeated, node_bearings, c_mat, k, lvec, mbb, kbb, cbb, gbb] = global_assembly(num_discs, dist_discs, dist_segments, d_segments, density, d_discs, l_discs, dist_bearings, l_d_ratio, E, num_bearings, Kxx, Kyy, Kxy, Kyx, Damp_xx, Damp_yy, Damp_xy, Damp_yx, lumped_mass, consistent_mass, G, k_sh, beam_th,disc_density)
    
    % Disc and general calculations

        for a=1:num_discs % disc counter
            clear lowest;

            if dist_discs(a)==0
                lowest = 1;
            else
                lowest = find(dist_discs(a) > dist_segments(1:end-1) & dist_discs(a) <= dist_segments(2:end), 1);
            end

            % temp = sort(abs( dist_discs(a) - dist_segments));   
            % lowest = find(abs( dist_discs(a) - dist_segments ) == temp(1) );
            % sec_lowest = find(abs( dist_discs(a) - dist_segments ) == temp(2) ); % finds the range in dist_segments where location of disc lies.

            shaft_dia_used(a) = d_segments(lowest); 
            md(a) = disc_density(a) * pi * (((d_discs(a)/2)^2) - (d_segments(lowest)/2)^2) * l_discs(a);  % mass of discs [kg]
            Ip(a) = (md(a)/8) * ((d_discs(a)^2)+(d_segments(lowest)^2)); % polar mass moment of inertia [kg m2]
            Id(a) = md(a)* ( (((d_discs(a)^2)+(d_segments(lowest)^2)) /16) + ((l_discs(a)^2)/12) ); % diametral transverse mass moment of inertia in kg m2 
            % Id = [1.924567 0.119595 0.08221606];
        end
        clear lowest;

    %% Meshing Starts

    dist_all_segments = [dist_segments, dist_bearings, dist_discs ];
    dist_all_segments = sort(dist_all_segments);
    dist_all_segments(1) = [];

    diff_dist_all_elements = [ dist_all_segments(1),diff(dist_all_segments) ];
    lvec = [];
    seg_dia_repeated = [];
    % Figure out lengths of each mesh group. Each is a collection
    % of elements. Change in group occurs @ location of discs and bearings.
    % Mesh adjustment ensures each element ends exactly in middle of disc or bearing
    % centre. Ceil function ensures l/d ratio of each element stays below 1.

    for a=1:length(diff_dist_all_elements) % segment by segment loop

    if dist_all_segments(a)==0
       lowest = 1;
    else
       lowest = find(dist_all_segments(a) > dist_segments(1:end-1) & dist_all_segments(a) <= dist_segments(2:end), 1);
    end

    l_oneelement = l_d_ratio*d_segments(lowest); % initial element length without mesh adjustment
    d_all_segments(a) = d_segments(lowest);
    num_elements(a) = ceil (diff_dist_all_elements(a)/l_oneelement); % vector of number of elements for each segment
    f_l_oneelement(a) = diff_dist_all_elements(a)/num_elements(a); % vector of final lengths of one element for each segment
    
    for seg=1:num_elements(a)
        seg_dia_repeated = [ seg_dia_repeated, d_all_segments(a)]; 
    end
    
    % define length vector for plotting, left-to-right numbering. Contains locations of nodes in terms of lenghts
    if a==1
        lvec = [ lvec, 0 : f_l_oneelement(a) : f_l_oneelement(a)*num_elements(a) ] ;
    else
    lvec = [ lvec, lvec(end)+f_l_oneelement(a) : f_l_oneelement(a) : lvec(end)+f_l_oneelement(a)*num_elements(a) ] ;
    end

    end

    all_num_elements = sum(num_elements); % total number of elements
    sum_num_elements = cumsum(num_elements); % vector of cumulative sums of segment number of elements

    % Meshing Ends*************************************************************

    % Guyan Reduction is not used
    num_plot_max = 2*all_num_elements;
    num_plot_default = all_num_elements;

%     num_plot = input(['enter the number of modes to plot, max�,' ...  '
%     num2str(num_plot_max),', default ',num2str(num_plot_default),' ... ']);
%     if (isempty(num_plot))
%     num_plot = num_plot_default;
%     else
%     end

    %% Global matrix formulations

    % define the node numbers
    % n = 1 : (final_all_num_elements+1);

    % number the nodes for the elements
    node1 = 1 : all_num_elements; % this vector contains first node of each element
    node2 = 2 : all_num_elements+1; % this vector contains second node of each element

    % size the stiffness and mass matrices to have 4 times the number of nodes
    % to allow for translation and rotation dof's for each node
    mnu = all_num_elements+1; % total no. of nodes
    k = zeros(4*mnu);
    m = zeros(4*mnu);
    mr = zeros(4*mnu);
    g = zeros(4*mnu);
    c_mat = zeros(4*mnu);
    
    % Building up the global stiffness, mass and gyroscopic matrices, element
    % by element. First only lower triangular part is made, then symmetry/ 
    % skew-symmetry is used to create the other part.

    % I, area, and mass per unit length of beam
    z=1; % to refer to first segment

    for i = 1 : all_num_elements
        dof1 = 4*node1(i)-3;   % index number of node displacement for a specific element i 
        dof2 = 4*node1(i)-2;   
        dof3 = 4*node1(i)-1;   
        dof4 = 4*node1(i);     

        dof5 = 4*node2(i)-3;   
        dof6 = 4*node2(i)-2;   
        dof7 = 4*node2(i)-1;   
        dof8 = 4*node2(i);     

    if (i > sum_num_elements(z) )
        z=z+1; % segment counter
    end

    if (i <= sum_num_elements(z) ) % this does not need to run for all elements in a segment as is the case now. Optimize!
        l = f_l_oneelement(z);
        dia_for_lumped = d_all_segments(z); 
        I = pi*((d_all_segments(z)/2)^4)/4;  % cross-sectional area moment of inertia [m^4]
        area = pi*((d_all_segments(z)/2)^2); % [m^2]
        mpl = density*area;  % [kg/m]
        j = ((d_all_segments(z)/2)^2)/4;
        j1 = 2*j;
        sh_form = (12*E*I*k_sh)/(G*area*(l^2)); % transverse shear effect [Nelson 1980]
%         sh_form = 0;
    end
    
    if (beam_th==0)
        % Only symm matrix is being defined here
        k(dof1,dof1) = k(dof1,dof1)+(12*E*I/(l^3));       % matrix element definition for stifness matrix for element i (all rows below)
        k(dof4,dof1) = k(dof4,dof1)+(6*E*I/(l^2));
        k(dof5,dof1) = k(dof5,dof1)+(-12*E*I/(l^3));
        k(dof8,dof1) = k(dof8,dof1)+(6*E*I/(l^2));

        k(dof2,dof2) = k(dof2,dof2)+(12*E*I/(l^3));
        k(dof3,dof2) = k(dof3,dof2)+(-6*E*I/(l^2));
        k(dof6,dof2) = k(dof6,dof2)+(-12*E*I/(l^3));
        k(dof7,dof2) = k(dof7,dof2)+(-6*E*I/(l^2));

        k(dof3,dof3) = k(dof3,dof3)+(4*E*I/l);
        k(dof6,dof3) = k(dof6,dof3)+(6*E*I/(l^2));
        k(dof7,dof3) = k(dof7,dof3)+(2*E*I/l);

        k(dof4,dof4) = k(dof4,dof4)+(4*E*I/l);
        k(dof5,dof4) = k(dof5,dof4)+(-6*E*I/(l^2));
        k(dof8,dof4) = k(dof8,dof4)+(2*E*I/l);

        k(dof5,dof5) = k(dof5,dof5)+(12*E*I/(l^3));
        k(dof8,dof5) = k(dof8,dof5)+(-6*E*I/(l^2));

        k(dof6,dof6) = k(dof6,dof6)+(12*E*I/(l^3));
        k(dof7,dof6) = k(dof7,dof6)+(6*E*I/(l^2));

        k(dof7,dof7) = k(dof7,dof7)+(4*E*I/l);

        k(dof8,dof8) = k(dof8,dof8)+(4*E*I/l);

        if (lumped_mass == 1)
            m(dof1,dof1) = m(dof1,dof1) + mpl*l/2;
            m(dof2,dof2) = m(dof2,dof2) + mpl*l/2;
            m(dof5,dof5) = m(dof5,dof5) + mpl*l/2;
            m(dof6,dof6) = m(dof6,dof6) + mpl*l/2;

            m(dof3,dof3) = m(dof3,dof3) + (1/24)*mpl*l* ( 3*(dia_for_lumped^2)/4 + (l/2)^2 ) + (mpl*l/2)*((l/4)^2);
            m(dof4,dof4) = m(dof4,dof4) + (1/24)*mpl*l* ( 3*(dia_for_lumped^2)/4 + (l/2)^2 ) + (mpl*l/2)*((l/4)^2);
            m(dof7,dof7) = m(dof7,dof7) + (1/24)*mpl*l* ( 3*(dia_for_lumped^2)/4 + (l/2)^2 ) + (mpl*l/2)*((l/4)^2);
            m(dof8,dof8) = m(dof8,dof8) + (1/24)*mpl*l* ( 3*(dia_for_lumped^2)/4 + (l/2)^2 ) + (mpl*l/2)*((l/4)^2);
        end

        if (consistent_mass == 1)
            m(dof1,dof1) = m(dof1,dof1)+(mpl/420)*(156*l)+       (mpl/(30*l))*j*36;       % matrix element definition for mass matrix for element i (all rows below)
            m(dof4,dof1) = m(dof4,dof1)+(mpl/420)*(22*(l^2))+    (mpl/(30*l))*j*(3*l);    % sum of rotational and translational mass matrices
            m(dof5,dof1) = m(dof5,dof1)+(mpl/420)*(54*l)+        (mpl/(30*l))*j*-36;
            m(dof8,dof1) = m(dof8,dof1)+(mpl/420)*(-13*(l^2))+   (mpl/(30*l))*j*(3*l);

            m(dof2,dof2) = m(dof2,dof2)+(mpl/420)*(156*l)+       (mpl/(30*l))*j*36;
            m(dof3,dof2) = m(dof3,dof2)+(mpl/420)*(-22*(l^2))+   (mpl/(30*l))*j*(-3*l);
            m(dof6,dof2) = m(dof6,dof2)+(mpl/420)*(54*l)+        (mpl/(30*l))*j*(-36);
            m(dof7,dof2) = m(dof7,dof2)+(mpl/420)*(13*(l^2))+    (mpl/(30*l))*j*(-3*l);

            m(dof3,dof3) = m(dof3,dof3)+(mpl/420)*(4*(l^3))+     (mpl/(30*l))*j*(4*(l^2));
            m(dof6,dof3) = m(dof6,dof3)+(mpl/420)*(-13*(l^2))+   (mpl/(30*l))*j*(3*l);
            m(dof7,dof3) = m(dof7,dof3)+(mpl/420)*(-3*(l^3))+    (mpl/(30*l))*j*(-(l^2));

            m(dof4,dof4) = m(dof4,dof4)+(mpl/420)*(4*(l^3))+     (mpl/(30*l))*j*(4*(l^2));
            m(dof5,dof4) = m(dof5,dof4)+(mpl/420)*(13*(l^2))+    (mpl/(30*l))*j*(-3*l);
            m(dof8,dof4) = m(dof8,dof4)+(mpl/420)*(-3*(l^3))+    (mpl/(30*l))*j*(-(l^2));

            m(dof5,dof5) = m(dof5,dof5)+(mpl/420)*(156*l)+       (mpl/(30*l))*j*36;
            m(dof8,dof5) = m(dof8,dof5)+(mpl/420)*(-22*(l^2))+   (mpl/(30*l))*j*(-3*l);

            m(dof6,dof6) = m(dof6,dof6)+(mpl/420)*(156*l)+       (mpl/(30*l))*j*36;
            m(dof7,dof6) = m(dof7,dof6)+(mpl/420)*(22*(l^2))+    (mpl/(30*l))*j*(3*l);

            m(dof7,dof7) = m(dof7,dof7)+(mpl/420)*(4*(l^3))+     (mpl/(30*l))*j*(4*(l^2));

            m(dof8,dof8) = m(dof8,dof8)+(mpl/420)*(4*(l^3))+     (mpl/(30*l))*j*(4*(l^2));
        end

        g(dof2,dof1) = g(dof2,dof1)+(mpl/(30*l))*j1*-36;       %matrix element definition for gyroscopic matrix for element i (all rows below)
        g(dof3,dof1) = g(dof3,dof1)+(mpl/(30*l))*j1*(3*l);     % According to Ishida. Signs are inverted in Nelson 1976.
        g(dof6,dof1) = g(dof6,dof1)+(mpl/(30*l))*j1*36;
        g(dof7,dof1) = g(dof7,dof1)+(mpl/(30*l))*j1*(3*l);

        g(dof4,dof2) = g(dof4,dof2)+(mpl/(30*l))*j1*(3*l);
        g(dof5,dof2) = g(dof5,dof2)+(mpl/(30*l))*j1*-36;
        g(dof8,dof2) = g(dof8,dof2)+(mpl/(30*l))*j1*(3*l);

        g(dof4,dof3) = g(dof4,dof3)+(mpl/(30*l))*j1*(-4*(l^2));
        g(dof5,dof3) = g(dof5,dof3)+(mpl/(30*l))*j1*(3*l);
        g(dof8,dof3) = g(dof8,dof3)+(mpl/(30*l))*j1*(l^2);

        g(dof6,dof4) = g(dof6,dof4)+(mpl/(30*l))*j1*(3*l);
        g(dof7,dof4) = g(dof7,dof4)+(mpl/(30*l))*j1*(-(l^2));

        g(dof6,dof5) = g(dof6,dof5)+(mpl/(30*l))*j1*-36; 
        g(dof7,dof5) = g(dof7,dof5)+(mpl/(30*l))*j1*(-3*l);

        g(dof8,dof6) = g(dof8,dof6)+(mpl/(30*l))*j1*(-3*l);

        g(dof8,dof7) = g(dof8,dof7)+(mpl/(30*l))*j1*(-4*(l^2));
    end
    
    if (beam_th==1) % Timoshenko beam elements [Nelson 1980]. Only symm matrix is being defined here
        cn = E*I/((l^3)*(1+sh_form)); % constant to be multiplied with each element below. 
        % This changes with each element so is not being taken out of element wise loop.
        k(dof1,dof1) = k(dof1,dof1)+ cn*12;       % matrix element definition for stifness matrix for element i (all rows below)
        k(dof4,dof1) = k(dof4,dof1)+ cn*6*l;
        k(dof5,dof1) = k(dof5,dof1)+ cn*-12;
        k(dof8,dof1) = k(dof8,dof1)+ cn*6*l;

        k(dof2,dof2) = k(dof2,dof2)+ cn*12;
        k(dof3,dof2) = k(dof3,dof2)+ cn*-6*l;
        k(dof6,dof2) = k(dof6,dof2)+ cn*-12;
        k(dof7,dof2) = k(dof7,dof2)+ cn*-6*l;

        k(dof3,dof3) = k(dof3,dof3)+ cn* ((4*l^2)+sh_form*l^2);
        k(dof6,dof3) = k(dof6,dof3)+ cn*6*l;
        k(dof7,dof3) = k(dof7,dof3)+ cn* ((2*l^2)-sh_form*l^2);

        k(dof4,dof4) = k(dof4,dof4)+ cn* ((4*l^2)+sh_form*l^2);
        k(dof5,dof4) = k(dof5,dof4)+ cn*-6*l;
        k(dof8,dof4) = k(dof8,dof4)+ cn* ((2*l^2)-sh_form*l^2);

        k(dof5,dof5) = k(dof5,dof5)+ cn*12;
        k(dof8,dof5) = k(dof8,dof5)+ cn*-6*l;

        k(dof6,dof6) = k(dof6,dof6)+ cn*12;
        k(dof7,dof6) = k(dof7,dof6)+ cn*6*l;

        k(dof7,dof7) = k(dof7,dof7)+ cn* ((4*l^2)+sh_form*l^2);

        k(dof8,dof8) = k(dof8,dof8)+ cn* ((4*l^2)+sh_form*l^2);

        
        cnm = mpl*l/(420*(1+sh_form)^2); % translational mass matrix constant to be multiplied with each element below. 
        cnr = mpl*((d_all_segments(z)/2)^2)/(120*l*(1+sh_form)^2);
        m(dof1,dof1) = m(dof1,dof1)+ cnm*(156 + sh_form*294 + 140*sh_form^2)+           cnr*36;       % matrix element definition for mass matrix for element i (all rows below)
        m(dof4,dof1) = m(dof4,dof1)+ cnm*(22*l + sh_form*38.5*l + 17.5*l*sh_form^2)+    cnr*(3*l + sh_form*-15*l);    % sum of rotational and translational mass matrices
        m(dof5,dof1) = m(dof5,dof1)+ cnm*(54 + sh_form*126 + 70*sh_form^2)+             cnr*-36;
        m(dof8,dof1) = m(dof8,dof1)+ cnm*(-13*l - sh_form*31.5*l - 17.5*l*sh_form^2)+   cnr*(3*l + sh_form*-15*l);

        m(dof2,dof2) = m(dof2,dof2)+ cnm*(156 + sh_form*294 + 140*sh_form^2)+           cnr*36;
        m(dof3,dof2) = m(dof3,dof2)+ cnm*(-22*l - sh_form*38.5*l - 17.5*l*sh_form^2)+   cnr*(-3*l + sh_form*15*l);
        m(dof6,dof2) = m(dof6,dof2)+ cnm*(54 + sh_form*126 + 70*sh_form^2)+             cnr*-36;
        m(dof7,dof2) = m(dof7,dof2)+ cnm*(13*l + sh_form*31.5*l + 17.5*l*sh_form^2)+    cnr*(-3*l + sh_form*15*l);

        m(dof3,dof3) = m(dof3,dof3)+ cnm*(4*(l^2) + sh_form*7*(l^2) + 3.5*(l^2)*sh_form^2)+     cnr*(4*(l^2) + sh_form*5*(l^2) + 10*(l^2)*sh_form^2);
        m(dof6,dof3) = m(dof6,dof3)+ cnm*(-12*l - sh_form*31.5*l - 17.5*l*sh_form^2)+           cnr*(3*l + sh_form*-15*l);
        m(dof7,dof3) = m(dof7,dof3)+ cnm*(-3*(l^2) - sh_form*7*(l^2) - 3.5*(l^2)*sh_form^2)+    cnr*(-(l^2) + sh_form*-5*(l^2) + 5*(l^2)*sh_form^2);

        m(dof4,dof4) = m(dof4,dof4)+ cnm*(4*(l^2) + sh_form*7*(l^2) + 3.5*(l^2)*sh_form^2)+     cnr*(4*(l^2) + sh_form*5*(l^2) + 10*(l^2)*sh_form^2);
        m(dof5,dof4) = m(dof5,dof4)+ cnm*(13*l + sh_form*31.5*l + 15.7*l*sh_form^2)+            cnr*(-3*l + sh_form*15*l);
        m(dof8,dof4) = m(dof8,dof4)+ cnm*(-3*(l^2) + sh_form*-7*(l^2) - 3.5*(l^2)*sh_form^2)+   cnr*(-(l^2) + sh_form*-5*(l^2) + 5*(l^2)*sh_form^2);

        m(dof5,dof5) = m(dof5,dof5)+ cnm*(156+sh_form*294+140*sh_form^2)+                       cnr*36;
        m(dof8,dof5) = m(dof8,dof5)+ cnm*(-22*l + sh_form*-38.5*l - 17.5*l*sh_form^2)+          cnr*(-3*l + sh_form*15*l);

        m(dof6,dof6) = m(dof6,dof6)+ cnm*(156+sh_form*294+140*sh_form^2)+                       cnr*36;
        m(dof7,dof6) = m(dof7,dof6)+ cnm*(22*l + sh_form*38.5*l + 17.5*l*sh_form^2)+            cnr*(3*l + sh_form*-15*l);
        
        m(dof7,dof7) = m(dof7,dof7)+ cnm*(4*(l^2) + sh_form*7*(l^2) + 3.5*(l^2)*sh_form^2)+     cnr*(4*(l^2) + sh_form*5*(l^2) + (sh_form^2)*10*l^2);

        m(dof8,dof8) = m(dof8,dof8)+ cnm*(4*(l^2) + sh_form*7*(l^2) + 3.5*(l^2)*sh_form^2)+     cnr*(4*(l^2) + sh_form*5*(l^2) + (sh_form^2)*10*l^2);
        
%         % Separate rotational mass matrix definition
%         mr(dof1,dof1) = mr(dof1,dof1)+   cnr*36;       % matrix element definition for mass matrix for element i (all rows below)
%         mr(dof4,dof1) = mr(dof4,dof1)+   cnr*(3*l + sh_form*-15*l);    % sum of rotational and translational mass matrices
%         mr(dof5,dof1) = mr(dof5,dof1)+   cnr*-36;
%         mr(dof8,dof1) = mr(dof8,dof1)+   cnr*(3*l + sh_form*-15*l);
% 
%         mr(dof2,dof2) = mr(dof2,dof2)+   cnr*36;
%         mr(dof3,dof2) = mr(dof3,dof2)+   cnr*(-3*l + sh_form*15*l);
%         mr(dof6,dof2) = mr(dof6,dof2)+   cnr*-36;
%         mr(dof7,dof2) = mr(dof7,dof2)+   cnr*(-3*l + sh_form*15*l);
% 
%         mr(dof3,dof3) = mr(dof3,dof3)+   cnr*(4*(l^2) + sh_form*5*(l^2) + 10*(l^2)*sh_form^2);
%         mr(dof6,dof3) = mr(dof6,dof3)+   cnr*(3*l + sh_form*-15*l);
%         mr(dof7,dof3) = mr(dof7,dof3)+   cnr*(-(l^2) + sh_form*-5*(l^2) + 5*(l^2)*sh_form^2);
% 
%         mr(dof4,dof4) = mr(dof4,dof4)+  cnr*(4*(l^2) + sh_form*5*(l^2) + 10*(l^2)*sh_form^2);
%         mr(dof5,dof4) = mr(dof5,dof4)+  cnr*(-3*l + sh_form*15*l);
%         mr(dof8,dof4) = mr(dof8,dof4)+  cnr*(-(l^2) + sh_form*-5*(l^2) + 5*(l^2)*sh_form^2);
% 
%         mr(dof5,dof5) = mr(dof5,dof5)+  cnr*36;
%         mr(dof8,dof5) = mr(dof8,dof5)+  cnr*(-3*l + sh_form*15*l);
% 
%         mr(dof6,dof6) = mr(dof6,dof6)+  cnr*36;
%         mr(dof7,dof6) = mr(dof7,dof6)+  cnr*(3*l + sh_form*-15*l);
%         
%         mr(dof7,dof7) = mr(dof7,dof7)+  cnr*(4*(l^2) + sh_form*5*(l^2) + (sh_form^2)*10*l^2);
% 
%         mr(dof8,dof8) = mr(dof8,dof8)+  cnr*(4*(l^2) + sh_form*5*(l^2) + (sh_form^2)*10*l^2);
                    
        % Friswell 2010 matrices. Shear factor in Nelson 1980 is missing. Also
        % Nelson 1980 matrices have inverted signs since equation has negative sign with
        % gyroscopic term.
        cng= mpl*((d_all_segments(z)/2)^2)/(60*l*(1+sh_form)^2); % [Friswell 2010]
        g(dof2,dof1) = g(dof2,dof1)+    cng*-36;       % matrix element definition for gyroscopic matrix for element i (all rows below)
        g(dof3,dof1) = g(dof3,dof1)+    cng*(3*l - 15*l*sh_form);     
        g(dof6,dof1) = g(dof6,dof1)+    cng*36;
        g(dof7,dof1) = g(dof7,dof1)+    cng*(3*l - 15*l*sh_form);

        g(dof4,dof2) = g(dof4,dof2)+    cng*(3*l - 15*l*sh_form);
        g(dof5,dof2) = g(dof5,dof2)+    cng*-36;
        g(dof8,dof2) = g(dof8,dof2)+    cng*(3*l - 15*l*sh_form);

        g(dof4,dof3) = g(dof4,dof3)+    cng*(-4*(l^2) - 5*(l^2)*sh_form - 10*(l^2)*(sh_form^2));
        g(dof5,dof3) = g(dof5,dof3)+	cng*(3*l - 15*l*sh_form);
        g(dof8,dof3) = g(dof8,dof3)+	cng*((l^2) + 5*(l^2)*sh_form - 5*(l^2)*(sh_form^2));

        g(dof6,dof4) = g(dof6,dof4)+	cng*(3*l - 15*l*sh_form);
        g(dof7,dof4) = g(dof7,dof4)+	cng*(-(l^2) - 5*(l^2)*sh_form + 5*(l^2)*(sh_form^2));

        g(dof6,dof5) = g(dof6,dof5)+	cng*-36;
        g(dof7,dof5) = g(dof7,dof5)+	cng*(-3*l + 15*l*sh_form);

        g(dof8,dof6) = g(dof8,dof6)+	cng*(-3*l + 15*l*sh_form);

        g(dof8,dof7) = g(dof8,dof7)+	cng*(-4*(l^2) - 5*(l^2)*sh_form - 10*(l^2)*(sh_form^2));
    end
    
    end

    % Completing global stiffness, global mass and global gyroscopic matrices using symmetry and skew-symmetry
    g = g - g'; % diag elements are zero in skew-sym matrix
    % m+m' means diag elements added twice. diag(m) extracts diag elements of m into a vector. Another diag converts diag elements of vector
    % into a diag matrix for subtraction
    if (consistent_mass ==1)
        m = m + m' - diag(diag(m));
    end

    k = k + k' - diag(diag(k));

    %% Addition of discs and supports in global matrices

    % Addition of disc mass and gyro matrices & bearing and coupling stiffness
    % matrices in global mass, gyro, & stiffness matrices. For 8DOF beam 
    % element, elements at locations 55,66,77,88 correspond to node 2. 
    % Elements at 99,1010,1111,1212, correspond to node % 3. So for nth node, 
    % locations are (4n-3)(4n-3),(4n-2)(4n-2),(4n-1)(4n-1), (4n)(4n).

    for a=1:num_discs;
        node_discs(a) = sum_num_elements( find( dist_discs(a) == dist_all_segments,1)) + 1;

    % Disc placement in global mass matrix [Ishida Pg. 330]
        m(4*node_discs(a)-3,4*node_discs(a)-3) = m(4*node_discs(a)-3,4*node_discs(a)-3) + md(a);
        m(4*node_discs(a)-2,4*node_discs(a)-2) = m(4*node_discs(a)-2,4*node_discs(a)-2) + md(a);

        m(4*node_discs(a)-1,4*node_discs(a)-1) = m(4*node_discs(a)-1,4*node_discs(a)-1) + Id(a);
        m(4*node_discs(a),4*node_discs(a))     = m(4*node_discs(a),4*node_discs(a)) + Id(a);

    % Disc placement in global gyro matrix.% -Ip is added in same row but a previous column
        g(4*node_discs(a),4*node_discs(a)-1) = g(4*node_discs(a),4*node_discs(a)-1) - Ip(a); 
        g(4*node_discs(a)-1,4*node_discs(a)) = g(4*node_discs(a)-1,4*node_discs(a)) + Ip(a);
    end

    % Let m be number of nodes for comments only
    kbb = k;     % this retain the k matrix in original form without additions of bearing stiffness
    % matrix of (4n+4)*(4n+4) order or (4m X 4m) order
    cbb = c_mat; % Retains c_mat in original form
    mbb = m;                % matrix of (4n+4)*(4n+4) order
    gbb = g;                % matrix of (4n+4)*(4n+4) order

    for a=1:num_bearings;

        node_bearings(a) = sum_num_elements( find( dist_bearings(a) == dist_all_segments,1)) + 1;

    % Bearing stiffness addition in global stiffness matrix
        kbb(4*node_bearings(a)-3,4*node_bearings(a)-3) = kbb(4*node_bearings(a)-3,4*node_bearings(a)-3) + Kxx(a); % these are all const support coeff
        kbb(4*node_bearings(a)-2,4*node_bearings(a)-2) = kbb(4*node_bearings(a)-2,4*node_bearings(a)-2) + Kyy(a);
        kbb(4*node_bearings(a)-2,4*node_bearings(a)-3) = kbb(4*node_bearings(a)-2,4*node_bearings(a)-3) + Kxy(a);
        kbb(4*node_bearings(a)-3,4*node_bearings(a)-2) = kbb(4*node_bearings(a)-3,4*node_bearings(a)-2) + Kyx(a);

        cbb(4*node_bearings(a)-3,4*node_bearings(a)-3) = cbb(4*node_bearings(a)-3,4*node_bearings(a)-3) + Damp_xx(a);
        cbb(4*node_bearings(a)-2,4*node_bearings(a)-2) = cbb(4*node_bearings(a)-2,4*node_bearings(a)-2) + Damp_yy(a);
        cbb(4*node_bearings(a)-2,4*node_bearings(a)-3) = cbb(4*node_bearings(a)-2,4*node_bearings(a)-3) + Damp_xy(a);
        cbb(4*node_bearings(a)-3,4*node_bearings(a)-2) = cbb(4*node_bearings(a)-3,4*node_bearings(a)-2) + Damp_yx(a);
    end
%     cbb=-cbb;
end