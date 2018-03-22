% Copyright (C) 2018 Kazi Sher Ahmed <kazisherahmed@gmail.com> & SM Ahmad.
% This file is a part of VibronRotor - A finite-element code for rotordynamic analysis of flexible rotor-bearing systems.
% VibronRotor is released under the terms of GNU General Public License 3.0.

function [node_discs, seg_dia_repeated, node_bearings, c_mat, k, lvec, mbb, kbb, cbb, gbb] = global_assembly(num_discs, dist_discs, dist_segments, d_segments, density, d_discs, l_discs, dist_bearings, l_d_ratio, E, num_bearings, Kxx, Kyy, Kxy, Kyx, Damp_xx, Damp_yy, Damp_xy, Damp_yx, lumped_mass, consistent_mass)
    
    % Disc and general calculations

        for a=1:num_discs
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
            md(a) = density * pi * (((d_discs(a)/2)^2) - (d_segments(lowest)/2)^2) * l_discs(a) ; % mass of discs in kg
            Ip(a) = (md(a)/8) * ((d_discs(a)^2)+(d_segments(lowest)^2)); % polar mass moment of inertia in kg m2
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
    num_elements(a) = ceil (diff_dist_all_elements(a)/l_oneelement);
    f_l_oneelement(a) = diff_dist_all_elements(a)/num_elements(a); % matrix of final lengths of one element for each segment
    
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

    all_num_elements = sum(num_elements);
    sum_num_elements = cumsum(num_elements);

    % Meshing Ends*************************************************************

    % Guyan Reduction is not used
    num_plot_max = 2*all_num_elements;
    num_plot_default = all_num_elements;

%     num_plot = input(['enter the number of modes to plot, max’,' ...  '
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
    mnu = all_num_elements+1;
    k = zeros(4*mnu);
    m = zeros(4*mnu);
    g = zeros(4*mnu);
    c_mat = zeros(4*mnu);
    % Building up the global stiffness, mass and gyroscopic matrices, element
    % by element. First only lower triangular part is made, then symmetry/ 
    % skew-symmetry is used to create the other part.

    % I, area and mass per unit length of beam
    z=1;

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
        z=z+1;
    end

    if (i <= sum_num_elements(z) )
        l = f_l_oneelement(z);
        dia_for_lumped = d_all_segments(z); 
        I = pi*((d_all_segments(z)/2)^4)/4;  % cross-sectional area moment of inertia in m^4
        area = pi*((d_all_segments(z)/2)^2); % in m^2
        mpl = density*area;  % in kg/m
        j = ((d_all_segments(z)/2)^2)/4;
        j1 = 2*j;
    end

    % Only symm matrix is being defined here
    k(dof1,dof1) = k(dof1,dof1)+(12*E*I/(l^3));       %matrix element definition for stifness matrix for element i (all rows below)
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
        m(dof1,dof1) = m(dof1,dof1)+(mpl/420)*(156*l)+       (mpl/(30*l))*j*36;       %matrix element definition for mass matrix for element i (all rows below)
        m(dof4,dof1) = m(dof4,dof1)+(mpl/420)*(22*(l^2))+    (mpl/(30*l))*j*(3*l);   %sum of rotational and translational mass matrices
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
    g(dof3,dof1) = g(dof3,dof1)+(mpl/(30*l))*j1*(3*l);     % According to Ishida.
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

    % making final global stiffness, global mass and global gyroscopic matrices
    % using symmetry and skew-symmetry

    g = g - g'; % diag elements are zero in skew-sym matrix

    % m+m' means diag elements added twice. diag(m) extracts diag elements
    % of m into a vector. Another diag converts diag elements of vector
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
    kbb = k;     % this retain the k matrix in original form without additions of bearing sitffness
    % matrix of (4n+4)*(4n+4) order or (4m X 4m) order
    cbb = c_mat; % Retains c_mat in original form
    mbb = m;                % matrix of (4n+4)*(4n+4) order
    gbb = g;                % matrix of (4n+4)*(4n+4) order

    for a=1:num_bearings;

        node_bearings(a) = sum_num_elements( find( dist_bearings(a) == dist_all_segments,1)) + 1;

    % Bearing stiffness addition in global stiffness matrix
        kbb(4*node_bearings(a)-3,4*node_bearings(a)-3) = kbb(4*node_bearings(a)-3,4*node_bearings(a)-3) + Kxx(a);
        kbb(4*node_bearings(a)-2,4*node_bearings(a)-2) = kbb(4*node_bearings(a)-2,4*node_bearings(a)-2) + Kyy(a);
        kbb(4*node_bearings(a)-2,4*node_bearings(a)-3) = kbb(4*node_bearings(a)-2,4*node_bearings(a)-3) + Kxy(a);
        kbb(4*node_bearings(a)-3,4*node_bearings(a)-2) = kbb(4*node_bearings(a)-3,4*node_bearings(a)-2) + Kyx(a);

        cbb(4*node_bearings(a)-3,4*node_bearings(a)-3) = cbb(4*node_bearings(a)-3,4*node_bearings(a)-3) + Damp_xx(a);
        cbb(4*node_bearings(a)-2,4*node_bearings(a)-2) = cbb(4*node_bearings(a)-2,4*node_bearings(a)-2) + Damp_yy(a);
        cbb(4*node_bearings(a)-2,4*node_bearings(a)-3) = cbb(4*node_bearings(a)-2,4*node_bearings(a)-3) + Damp_xy(a);
        cbb(4*node_bearings(a)-3,4*node_bearings(a)-2) = cbb(4*node_bearings(a)-3,4*node_bearings(a)-2) + Damp_yx(a);
    end
end