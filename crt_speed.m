% Copyright (C) 2018 Kazi Sher Ahmed <kazisherahmed@gmail.com> & SM Ahmad.
% This file is a part of VibronRotor - A finite-element code for rotordynamic analysis of flexible rotor-bearing systems.
% VibronRotor is released under the terms of GNU General Public License 3.0.

% function crt_speed(support_coeff_speed, kxx_speed, kyy_speed,kxy_speed,kyx_speed, cxx_speed, cxy_speed, cyy_speed, cyx_speed, k, c_mat, num_bearings, node_bearings, mbb, gbb )
% 
%     clear A evec1 evalu evalud evaludhz firstcolumns evalorder indexhz evalhzr;
% 
%     ana_range = 10000; % analysis range in rpm
%     b=1;
%     I = eye(length(mbb));      % identity matrix of (4n+4)*(4n+4) order or (4m X 4m)
%     O = zeros(length(mbb));    % null matrix of (4n+4)*(4n+4) order 
% 
% 
%     for arpm = 0:1000:6000
%         arad = (arpm*pi)/30;   % rotational velocity of rotor in rad/s
%         a(1,b) = arpm;
% 
%         %Cubic spline curve fitting
%         kxx_fit(1,b) = spline(support_coeff_speed,kxx_speed,arpm);
%         kyy_fit(1,b) = spline(support_coeff_speed,kyy_speed,arpm);
%         kxy_fit(1,b) = spline(support_coeff_speed,kxy_speed,arpm);
%         kyx_fit(1,b) = spline(support_coeff_speed,kyx_speed,arpm);
% 
%         cxx_fit(1,b) = spline(support_coeff_speed,cxx_speed,arpm);
%         cyy_fit(1,b) = spline(support_coeff_speed,cyy_speed,arpm);
%         cxy_fit(1,b) = spline(support_coeff_speed,cxy_speed,arpm);
%         cyx_fit(1,b) = spline(support_coeff_speed,cyx_speed,arpm);
% 
%         clear k_temp c_temp;
%         k_temp = k;
%         c_temp = c_mat;
% 
%         for num_b = 1:num_bearings
%             % Bearing nodes have already been located. Bearing stiffness addition in global stiffness matrix. Order below is xx,yy,xy,yx.
%             k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) = k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) + kxx_fit(1,b);
%             k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) = k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) + kyy_fit(1,b);
%             k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) = k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) + kxy_fit(1,b);
%             k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) = k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) + kyx_fit(1,b);
% 
%             c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) = c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) + cxx_fit(1,b);
%             c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) = c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) + cyy_fit(1,b);
%             c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) = c_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-3) + cxy_fit(1,b);
%             c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) = c_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-2) + cyx_fit(1,b);
%         end
% 
%         A = [ O I ; (-inv(mbb)*k_temp) (-inv(mbb)*(c_temp+arad*gbb)) ];
%         [evec1,evalu] = eig(A);    
%         evalud = diag(evalu);
% 
%     %     for_deletion = find(imag(evalud) == 0);
%     %     evalud(for_deletion,:) = [];
% 
%         evaludhz = evalud/(2*pi);    
% 
%         firstcolumns = 1:2:length(evaludhz);
%         evaludhz(firstcolumns) = []; % 4mX1
%         [evalorder,indexhz] = sort(abs((evaludhz)));
% 
%         for cnt = 1:length(evaludhz)
%             evalhzr(cnt,b) = (evaludhz(indexhz(cnt)));
%         end
%         text
%         damped_freq (:,b) = -imag (evalhzr(1:6,b));  
%         b = b+1;
%     end
% 
%         subplot(1,1,1)
%         c_modes_1b = plot(a, damped_freq(1,:)*60, ('-'),'LineWidth',1); hold on;
%         c_modes_2b = plot(a, damped_freq(3,:)*60, ('-'),'LineWidth',1); hold on;
%         c_modes_3b = plot(a, damped_freq(5,:)*60, ('-'),'LineWidth',1); hold on;
% 
%         c_modes_1f = plot(a, damped_freq(2,:)*60, ('-'),'LineWidth',1); hold on;
%         c_modes_2f = plot(a, damped_freq(4,:)*60, ('-'),'LineWidth',1); hold on;
%         c_modes_3f = plot(a, damped_freq(6,:)*60, ('-'),'LineWidth',1); hold on;
%         legend([c_modes_1b c_modes_2b c_modes_3b c_modes_1f c_modes_2f c_modes_3f], {['1b'], ['2b'], ['3b'], ['1f'], ['2f'], ['3f']}, 'FontSize',9,'Location','northeast' );
% 
%         title(['Critical Speed Map']);
%         xlabel('Bearing Stiffness (N/m)')
%         ylabel('Critical Speed (RPM)')
%     %     axis([ 0 ana_range 0 (ana_range)+8100 ])
% 
%         grid off
%         set(gca,'box','off');    
%         set(gcf,'color','w');
%         set(gca,'fontsize',9)
% 
%         disp('Press Enter to continue to the next selected functionality.'); pause

    
% Critical Speed Map - Without spline
function crt_speed(crt_rpm_range, crt_rpm_inc, dyn_k, k, num_bearings, node_bearings, mbb, gbb,cbb )
  
  
  I = eye(length(mbb));      % identity matrix of (4n+4)*(4n+4) order or (4m X 4m)
  O = zeros(length(mbb));    % null matrix of (4n+4)*(4n+4) order

  a=1;
  clear A b evec1 evalu;

  for dyn_k = dyn_k
      
      b(1,a) = dyn_k;
      clear k_temp;
      k_temp = k;
      
      for num_b=1:num_bearings
        % Bearing nodes have already been located. Bearing stiffness addition in global stiffness matrix:
        k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) = k_temp(4*node_bearings(num_b)-3,4*node_bearings(num_b)-3) + dyn_k;
        k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) = k_temp(4*node_bearings(num_b)-2,4*node_bearings(num_b)-2) + dyn_k;
      end
      
      c=1;
      clear d arad A ;

      for arpm = 0:crt_rpm_inc:crt_rpm_range
        arad = (arpm*2*pi)/60;   % rotational velocity of rotor in rad/s

        d(1,c) = arpm;

        A = [ O I ; (-mbb\k_temp) (-mbb\(cbb+arad*gbb)) ];
        [evec1,evalu] = eig(A);    
        evalud = diag(evalu);
        evaludhz = evalud/(2*pi);    

        firstcolumns = 1:2:(2*length(mbb));
        evaludhz(firstcolumns) = []; % 4mX1
        [evalorder,indexhz] = sort(abs((evaludhz)));

        for cnt = 1:length(evaludhz)
          evalhzr(cnt,c) = (evaludhz(indexhz(cnt)));
        end

        damped_freq (:,c) = -imag (evalhzr(1:6,c));
        c = c+1;
      end
         
      damped_freq = damped_freq.*60;
      
      for pt = 1:length(damped_freq)
        linesdnf_1f(pt,:) =  [ d(pt) damped_freq(2,pt) ] ;
        linesdnf_2f(pt,:) =  [ d(pt) damped_freq(4,pt) ] ;
      end
    x1 = [d(1) d(1) ; d(end) d(end)];
    mode1_intersect_f(a,:) = intersectPolylines (linesdnf_1f,x1);
    mode2_intersect_f(a,:) = intersectPolylines (linesdnf_2f,x1);
          
  a=a+1;
  end
  
  figure
  subplot(1,1,1)
  loglog (b, mode1_intersect_f(:,1), ('-'),'LineWidth',1); hold on;
  loglog (b, mode2_intersect_f(:,1), ('-'),'LineWidth',1); hold on;
  hold off;

  title(['Critical Speed Map']);
  xlabel('Bearing Stiffness (N/m)')
  ylabel('Forward Synchronous Critical Speed (RPM)')

  grid on;
  set(gca,'box','off');    
  set(gcf,'color','w');
  set(gca,'fontsize',9);   

  disp('Press Enter to continue to the next selected functionality.'); pause
end