% VibronRotor - A finite-element code for rotordynamic analysis of flexible rotor-bearing systems.
% Authors: Kazi Sher Ahmed, Prof. S. M. Ahmad.
% Release of VibronRotor source code is licensed under GNU GPL 3.0.

clear all;
format long;

%% User Inputs [SI unit version [metres, kilograms, seconds]

    % Beam theory selection [0 for Euler Bernoulli, 1 for Timoshenko]
        beam_th = 1;

    % Mass matrix formulation [0 is for No, 1 is for Yes]
        lumped_mass     = 0;
        consistent_mass = 1;

    % Bearing coefficients configuration [0 is for No, 1 is for Yes]
        const_supp_coeff   = 1;   % Constant [keep it Yes always]
        vari_supp_coeff    = 1;   % Speed-Dependent [only for instability threshold and critical speed map]

    % User-defined parameters
        w1        = 0;      % rpm of rotor
        l_d_ratio = 0.4;    % element length-to-diameter ratio

    % Rotor geometrical properties
        num_segments = 7;                                % number of segments
        l_segments   = 0.0254 * [ 4 6 10 2 4 2 4 ];      % length of segments [m]
        d_segments   = 0.0254 * [ 1 2 3 2 1.5 1.75 2 ];  % diameter of segments [m]

        num_discs    = 3;                              % number of load discs
        l_discs      = 0.0254 * [ 2 2 2 ];             % length of discs [m]
        d_discs      = 0.0254 * [ 22 11 10 ];          % diameter of discs [m]
        dist_discs   = 0.0254 * [ 15 17 31 ];          % distance of disc centres from left-most end of shaft [m]
        disc_density = 27679.9 * [0.283 0.283 0.283];  % densities of discs

        num_bearings  = 2;                  % number of bearing supports
        dist_bearings = 0.0254 * [ 3 25 ];  % distance of bearing centres from left-most end of shaft [m]

    % Rotor mechanical properties
       if (const_supp_coeff == 1)               % select if constant bearing coefficients 
            Kxx = 175.126835* [31.6e3 31.6e3];  % bearing stiffness in N/m
            Kyy = 175.126835* [31.6e3 31.6e3];
             
%             Kxx = [0.00001 0.00001];  % bearing stiffness in N/m
%             Kyy = [0.00001 0.00001];
            Kxy = 175.126835* [0 0];
            Kyx = 175.126835* [0 0];

            Damp_xx = 175.126835* [10 10];      % bearing damping in Ns/m
            Damp_yy = 175.126835* [10 10];
             
%             Damp_xx = [0.00001 0.00001];      % bearing damping in Ns/m
%             Damp_yy = [0.00001 0.00001];
            Damp_xy = 175.126835* [0 0];
            Damp_yx = 175.126835* [0 0];
        end

        if (vari_supp_coeff == 1)   % select if speed-dependent bearing coefficients
            kxx_speed = 175.126835* [ 9.00E+03 8.00E+03 7.00E+03 6.00E+03 5.00E+03 4.00E+03 ];
            kyy_speed = 175.126835* [ 9.00E+03 8.00E+03 7.00E+03 6.00E+03 5.00E+03 4.00E+03 ];
            kxy_speed = 175.126835* [ 0 0 0 0 0 0 ];
            kyx_speed = 175.126835* [ 0 0 0 0 0 0 ];

            cxx_speed = 175.126835* [ 10 30 50 70 90 110 ];
            cyy_speed = 175.126835* [ 10 30 50 70 90 110 ];
            cxy_speed = 175.126835* [ 0 0 0 0 0 0 ];
            cyx_speed = 175.126835* [ 0 0 0 0 0 0 ];

            support_coeff_speed = [ 1 1000 3000 5000 7000 10000 ];
        end

        E       = 6894.76 * 30e6;           % modulus of elasticity in N/m2
        density = 27679.9 * 0.283;          % density of shaft material [kg/m3]
        G       = 6894.76 * 12e6;           % Shear modulus [Pa]
        k_sh    = 1.128;                    % shear form factor for circular cross section

    % Functionalities [0 is for No, 1 is for Yes]
        plot_rotor         = 1;
        mesh_plotting      = 1;
            numbering      = 1;
        
        campbell_diag      = 1;
            cd_range       = 6000;  % analysis range in rpm
            increments     = 500;   % increments to calculate results after
            num_modes      = 3;     % number of modes to plot

        combined_modes     = 1;
            num_plot       = 3;  % number of modes to plot
        
        imb_resp_control   = 1;
            ana_range      = 6000;      % analysis range in rpm
            imb_interval   = 10;        % steps to calculate results
            ini_imb        = [ 10 2 ];  % initial imbalances in g-in
            ini_phase      = [ 0 -90 ]; % phase differences in imbalance masses in degrees
            fx_node_in     = [27 48];   % input nodes
            fx_node_out    = [27];      % output node
            orb_wrpm       = 2750;

        instability_threshold = 1;
            th_range          = 5000;   % analysis range [rpm]
            th_increments     = 1000;   % steps
            th_modes          = 3;      % number of modes to plot

        critical_speed_map = 1;
            crt_rpm_range  = 12000;
            crt_rpm_inc    = 1000;
            dyn_k          = 175.127 .* [ 10000 31623 100000 316228 1000000 3162278 10000000 2e7 3e7 4e7 ];    % in N/m
                       
%% Code Flow [ users do not need to edit below ]
    
    % Initial calculations
        dist_segments = cumsum(l_segments);
        dist_segments = [ 0,dist_segments ];

        overall_length = sum(l_segments);
        lbeam = overall_length;

        w = (w1*2*pi)/60;   % rotational velocity of rotor [rad/s]

    % Function calls
         if (plot_rotor==1)
             rotor_plot(numbering, l_segments, dist_segments, d_segments, num_discs, l_discs, d_discs, dist_discs, num_bearings, dist_bearings);
         end
         [node_discs, seg_dia_repeated, node_bearings, c_mat, k, lvec, mbb, kbb, cbb, gbb] = global_assembly(num_discs, dist_discs, dist_segments, d_segments, density, d_discs, l_discs, dist_bearings, l_d_ratio, E, num_bearings, Kxx, Kyy, Kxy, Kyx, Damp_xx, Damp_yy, Damp_xy, Damp_yx, lumped_mass, consistent_mass,G,k_sh,beam_th,disc_density);
         if (mesh_plotting==1)
             mesh_plot(numbering, lvec, seg_dia_repeated, l_segments, dist_segments, d_segments, num_discs, l_discs, d_discs, dist_discs, num_bearings, dist_bearings);
         end
         if (campbell_diag==1)
            freq_interference(num_modes, cd_range, increments, mbb, kbb, cbb, gbb);
         end
         if(combined_modes==1)
            [evaludhz, evalhzr] = mode_shapes(num_plot, lvec, w, mbb, kbb, cbb, gbb, num_bearings, dist_bearings, l_segments, dist_segments, d_segments, num_discs, l_discs, dist_discs, d_discs, lbeam);
         end
         if(instability_threshold==1)
            instab_threshold(th_modes, th_increments, th_range, support_coeff_speed, kxx_speed, kyy_speed,kxy_speed,kyx_speed, cxx_speed, cxy_speed, cyy_speed, cyx_speed, k, c_mat, num_bearings, node_bearings, mbb, gbb );
         end
         if(imb_resp_control==1)
            imb_resp(orb_wrpm, fx_node_out, fx_node_in, ini_phase, ini_imb, imb_interval, ana_range, node_discs, kbb, mbb, cbb, gbb);
         end
         if(critical_speed_map==1)
%             crt_speed(support_coeff_speed, kxx_speed, kyy_speed,kxy_speed,kyx_speed, cxx_speed, cxy_speed, cyy_speed, cyx_speed, k, c_mat, num_bearings, node_bearings, mbb, gbb);
              crt_speed(crt_rpm_range, crt_rpm_inc, dyn_k, k, num_bearings, node_bearings, mbb, gbb,cbb );

         end

