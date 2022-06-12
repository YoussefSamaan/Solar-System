%--------------------------------------------------
%
%                  Solar System
%
%--------------------------------------------------

clc
clear
close all

disp("Type the number of the test run you want to try: \n1. Normal/ real conditions \n2. To show how collisions look like when they happen. \n3. To show how the first 4 planets interact with the sun if they got really close to it. \n4. To show how the first 3 planets would interact with the sun if they were close to the sun but not as close to trial no. 3 (Similar interaction to what would happen if alpha particles were shot onto a thin gold foil) and the forth planet heading directly towards the sun.");
n = input("Test: ");
if n == 1
  load 'Test1.mat'
elseif n == 2
  load 'Test2.mat'
elseif n == 3
  load 'Test3.mat'
elseif n == 4
  load 'Test4.mat'
endif

% constants
G = 6.67408*10^(-11); % (N*m^2/Kg^2)

%equations
Ac = @(m,r) G.*m./r.^2; %gravitational acceleration
distv = @(d,v,t) d+v*t; %Final distance from the velocity
spda = @(v,a,t) v+a*t; %Final speed from the acceleration
dist = @(x1,y1,x2,y2) sqrt((x2-x1)^2 + (y2-y1)^2); % distane formula (pythagorian equation)
theta = @(x1,y1,x2,y2) atan2((y2-y1),(x2-x1)); % calculating the angle

%matrixes
Masses = [Sun_mass;
          Mercury_mass;
          Venus_mass; 
          Earth_mass;
          Mars_mass;
          Jupiter_mass; 
          Saturn_mass; 
          Uranus_mass; 
          Neptune_mass];
       
positions_of_planets = [Sun_radius,    0;
                       Mercury_radius, 0;
                       Venus_radius,   0;
                       Earth_radius,   0;
                       Mars_radius,    0;
                       Jupiter_radius, 0;
                       Saturn_radius,  0;
                       Uranus_radius,  0;
                       Neptune_radius, 0];
                       

Velocities_of_planets = [0,     Sun_Vi;
                         0, Mercury_Vi;
                         0,   Venus_Vi;
                         0,   Earth_Vi;
                         0,    Mars_Vi;
                         0, Jupiter_Vi;
                         0,  Saturn_Vi;
                         0,  Uranus_Vi;
                         0, Neptune_Vi];
                         
Disp_properties = [[1 1 0],5;
                  [0 0 0],1;
                  [1 1 1],2;
                  [0 0 1],2;
                  [1 0 0],1;
                  [1 1 0],4;
                  [1 0 1],4;
                  [0 1 1],3;
                  [0 0 1],3];% color and the size of the particles
    
Relative_positions = zeros(length(positions_of_planets));
Angles = zeros(length(positions_of_planets));
one = ones(length(positions_of_planets),1); % ones colomn vector 
forces_on_planets = zeros(length(positions_of_planets));
forces_on_planetsy = zeros(length(positions_of_planets));
forces_on_planetsx = zeros(length(positions_of_planets));
total_forces_y = zeros(length(positions_of_planets),1);
total_forces_x = zeros(length(positions_of_planets),1);
%{  
----------------> for the gif <--------------------
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'solar_system.gif';
%}
crashed = false;
    
% Animation
for i = 0:1000

  %Time (1 day) in seconds
  t = 24*60*60;
  
  for z = 1:length(positions_of_planets)
    for j = 1:length(positions_of_planets)
      
      %calculating the distance beetween the velocities
      Relative_positions(z,j) = dist(positions_of_planets(z,1),positions_of_planets(z,2),positions_of_planets(j,1),positions_of_planets(j,2));
      
      % check if a collision has occured
      if (Relative_positions(z,j) < (Planets_radii(z,1) + Planets_radii(j,1))) == (Relative_positions(z,j)~=0)
        
        % changing the velocity and the mass of one of the planet
        Velocities_of_planets(z,1) = (Velocities_of_planets(z,1)*Masses(z,1)+Velocities_of_planets(j,1)*Masses(j,1))/(Masses(z,1)+Masses(j,1));
        Velocities_of_planets(z,2) = (Velocities_of_planets(z,2)*Masses(z,1)+Velocities_of_planets(j,2)*Masses(j,1))/(Masses(z,1)+Masses(j,1));
        Masses(z,1) = Masses(z,1) + Masses(j,1);
        
        crashed = true; 
        
        % resizing the matrixes
        Masses(j,:) = [];
        positions_of_planets(j,:) = [];
        Planets_radii(j,:) = [];
        Velocities_of_planets(j,:) = [];
        Relative_positions(:,j) = [];
        Angles(:,j) = [];
        one(j,:) = [];
        forces_on_planets(:,j) = [];
        forces_on_planetsy(:,j) = [];
        forces_on_planetsx(:,j) = [];
        total_forces_y(j,:) = [];
        total_forces_x(j,:) = [];
        Disp_properties(j,:) = [];
        
      endif
      
      % to quit the loop if a collision has occured
      if crashed == true
        break
      endif
      
      % calculating the angles between the planets
      Angles(z,j) = theta(positions_of_planets(z,1),positions_of_planets(z,2),positions_of_planets(j,1),positions_of_planets(j,2));
      
      % calculating the acceleration
      if Relative_positions(z,j) ~= 0
        forces_on_planets(z,j) = Ac(Masses(j),Relative_positions(z,j));
      endif
      
    endfor
    
    % to quit the loop if a collision has occured
    if crashed == true
      break
    endif
    
  endfor
  
  % to stop breaking the loop and go to the next iteration
  if crashed == true
     crashed = false;
     continue
  endif
  
  % calculating the forces in the x and y directions
  forces_on_planetsy = forces_on_planets.*sin(Angles);
  forces_on_planetsx = forces_on_planets.*cos(Angles);
  
  % Finding the sum of the forces in the x and y directions
  total_forces_y = forces_on_planetsy*one;
  total_forces_x = forces_on_planetsx*one;
  
  for i = 1:length(positions_of_planets)
    
    % finding the new velocities in the x and y directions
    Velocities_of_planets(i,1) = spda(Velocities_of_planets(i,1),total_forces_x(i),t);
    Velocities_of_planets(i,2) = spda(Velocities_of_planets(i,2),total_forces_y(i),t);
    
    % finding the new position in the x and y directions
    positions_of_planets(i,1) = distv(positions_of_planets(i,1),Velocities_of_planets(i,1),t);
    positions_of_planets(i,2) = distv(positions_of_planets(i,2),Velocities_of_planets(i,2),t);
    
    % plotting the planets      
    plot(positions_of_planets(i,1),positions_of_planets(i,2),"Marker","o","MarkerFaceColor",[Disp_properties(i,1:3)],"MarkerSize",Disp_properties(i,4))
    hold on
  endfor
  
  hold off
  xlim([-5e12 5e12])
  ylim([-5e12 5e12])
  pause(1/120)
%{
  -------------------> for the gif <--------------------
  drawnow
  % Capture the plot as an image 
  frame = getframe(h); 
  im = frame2im(frame); 
  [imind,cm] = rgb2ind(im); 
  % Write to the GIF File 
  if i == 0 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.01,'Quality',100); 
  else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01,'Quality',100); 
  end 
  
%}
endfor

