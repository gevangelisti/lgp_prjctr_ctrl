function  robot = plnr_so_ro( n )

% planar  create an n-link planar robot.
% planar(n)  creates an all-revolute n-link planar robot with identical
% links.  For efficient use with Simulink, this function stores a copy of
% the last robot it created, and returns the copy if the arguments match.

% persistent  last_robot;% last_n;

% if isempty(last_robot) ~= 0 %&& last_n == n
%   robot = last_robot;
%   return
% end

robot.NB = n;
robot.parent = [0:n-1];

l_n = 1/n;
m_n = 1/n;
r_n = 0;%0.1;%1/n;

for i = 1:n
  robot.jtype{i} = 'r';
  if i == 1
    robot.Xtree{i} = plnr( 0, [0 0]);
  else
    robot.Xtree{i} = plnr( 0, [l_n 0]);
  end
  robot.I{i} = mcI( m_n, [l_n/2 0], m_n*(1/12*l_n^2 + 1/4*r_n^2) );
end

robot.gravity = 9.81*[1;0];

robot.appearance.base = ...
  { 'tiles', [-3 3; -3 3; -0.2 -0.2], 0.25, ...
  'box', [-0.05 -0.075 -0.05; 0.05 0.075 -0.01] };
%   { 'box', [-0.2 -0.3 -0.2; 0.2 0.3 -0.07] };

for i = 1:n
  robot.appearance.body{i} = ...
    { 'box', [0 -0.01 -0.01; l_n 0.01 0.01], ...
      'cyl', [0 0 -0.01; 0 0 0.01], 0.01 };
%     { 'box', [0 -0.07 -0.04; l_n 0.07 0.04], ...
%       'cyl', [0 0 -0.07; 0 0 0.07], 0.1 };
end

% last_robot = robot;
% %last_n = n;
