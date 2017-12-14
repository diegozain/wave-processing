function [x,z,dx,dz,t,dt,r,dr] = w_xztr(epsilon_w,fo,x,z,t)

% -------------------------------------------------------
%
%   example
%
% -------------------------------------------------------
% % known space interval of target space
%
% x = [0 20];
% z = [x(1) 3];
% t = total_time;
%
% % get x,z,t that will go into synthetic-wave,
% % and suggested r to perform experiment.
%
% [x,z,dx,dz,t,dt,r,dr] = w_xztr(epsilon_w,fo,x,z,t);

  % the speed of ---**LIGHT··))
  %
  eo = 8.854187817e-12;
  muo = 4*pi*1e-7;
  c = 1/sqrt(muo*eo);

  e_max = max(epsilon_w(:));
  e_max = max(1,e_max);
  e_min = min(epsilon_w(:));
  e_min = min(1,e_min);

  % the max/min speed of our target model
  %
  vel_min = c/sqrt( e_max );
  vel_max = c/sqrt( e_min );

  % smallest wavelength
  %
  fmax = 2.2*fo;
  l_min = vel_min / fmax;
  
  % dr, dx, dz
  %
  dr = l_min / 2;
  dx = l_min / 10;
  dz = l_min / 10;
  
  % r, x, z
  %
  r = x(1):dr:x(2);
  x = r(1):dx:r(end);
  z = r(1):dz:z(2);
  
  % dt
  %
  courant_factor = 0.9;
  dt = 1/(vel_max * sqrt((1/dx^2)+(1/dz^2))); 
  dt = courant_factor * dt;
  
  % t
  %
  t = 0 : (dt*1e+9) : t-dt; 

end