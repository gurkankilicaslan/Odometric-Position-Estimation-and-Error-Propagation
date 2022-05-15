close all; clear; clc; 

%%Simulation parameters
dt = 0.1; %step size
totalTime = 10; %Simulation time
t = 0:dt:totalTime; %Time span

%%Vehicle parameters
wheelR = 0.05; %radius of the wheel
lb = 0.1; %distance between wheel frame to vehicle frame
b = 2*lb;

%%Initial Conditions
x0 = 0; % initial x, you can change desirably
y0 = 0; % initial y, you can change desirably

alpha = -pi/2; 
beta = 0;


zeta0 = [x0;y0;0]; % initializing pose matrix

zeta(:,1) = zeta0; % initializing pose matrix

%%Loop starts here
    for i = 1:length(t)-1
        theta = zeta(3,i);  % current orientation in radius.
        RotationMatrix = [cos(theta), sin(theta), 0;
                 -sin(theta), cos(theta), 0;
                 0, 0, 1]; % rotation matrix
         
        %%inputs
        psi_dot_1 = 4; % right wheel angular velocity
        psi_dot_2 = 4; % left wheel angular velocity
        
        psi_dot = [psi_dot_1; psi_dot_2];  
        
          A = [-sin(alpha+beta), cos(alpha+beta), lb*cos(beta);
               -sin(alpha+beta), cos(alpha+beta), -lb*cos(beta);
               cos(alpha+beta), sin(alpha+beta), lb*sin(beta);
               cos(alpha+beta), sin(alpha+beta), -lb*sin(beta);];
           % first two rows are rolling constraints
           % last two rows are no-sliding constraints
           
          B = [wheelR, 0;
               0, wheelR;
               0, 0;
               0, 0;];
           
          zeta_dot(:,i) = inv(RotationMatrix)*inv(transpose(A)*A)...
              *transpose(A)*B*psi_dot; % zeta_dot derivation
        
          zeta(:,i+1) = zeta(:,i) + dt * zeta_dot(:,i); %updating pose
          
          dtheta = zeta(3,i+1) - zeta(3,i);
          ds = dtheta*wheelR;
          ds_r = dtheta*(wheelR + b/2);
          ds_l = dtheta*(wheelR - b/2);
          
          
    end

%%Plotting Functions

%%Animation (mobile robot motion animation)
l = 0.3; %length of the mobile robot
w = 2*lb; %width of the mobile robot


fh = figure();
fh.WindowState = 'maximized';

Covariance_Matrix=[0 0 0; 0 0 0; 0 0 0];

for i = 1:length(t)-1 %animation starts here
    theta = zeta(3,i);
             
    dtheta = zeta(3,i+1) - zeta(3,i);
    ds = dtheta*wheelR;
    ds_r = psi_dot_1*wheelR*dt; %dtheta*(wheelR + b/2);
    ds_l = psi_dot_2*wheelR*dt; %dtheta*(wheelR - b/2);
    
    
    grid on
    %axis([-2 10 -2 4]), 
    axis equal
    plot(zeta(1,1:i),zeta(2,1:i),'b-','LineWidth',1.8) % plotting TurtleBot's motion
    set(gca,'fontsize',12)
    xlabel('x,[m]');
    ylabel('y,[m]');
    title("Odometric Position Estimation and Error Propagation")
    
    
    
    [x,y,theta,Covariance_Matrix] = Odometry_f(zeta(1,i),zeta(2,i),theta,zeta(3,i+1),zeta(3,i),ds_r,ds_l,Covariance_Matrix);
        
    hold on
    
    if mod(i,10)==0
        
        C_x = Covariance_Matrix(1,1).*cos(0:0.01:2*pi);
        C_y = Covariance_Matrix(2,2).*sin(0:0.01:2*pi);
        plot(x + C_x.*cos(theta) - C_y.*sin(theta), y + C_x.*sin(theta) + C_y.*cos(theta),'-.','LineWidth',1.5);
    end
        
    pause(0.00000001);
    
    
end %animation ends

function [x,y,theta,Cp] = Odometry_f(x,y,theta,first,second,ds_r,ds_l,Cp)

lb = 0.1; %distance between wheel frame to vehicle frame
b = 2*lb;

dtheta = first-second;
%ds = dtheta*wheelR;
% ds_r = dtheta*(wheelR + b/2);
% ds_l = dtheta*(wheelR - b/2);

ds = (ds_r + ds_l)/2;

dx = ds*cos(theta + dtheta/2);
dy = ds*sin(theta + dtheta/2);

k_r = 0.01; % you can change desirably
k_l = 0.01; % you can change desirably

%Calculation of Covariance Matrix
Covariance_Matrix = [k_r*abs(ds_r) 0;
                    0 k_l*abs(ds_l)];

%Forward Kinematics Differentiation
d_position = [1 0 -ds*sin(theta + dtheta/2); 0 1 ds*cos(theta + dtheta/2); 0 0 1];

d_wheel = [0.5*cos(theta+dtheta/2)-(ds/(2*b))*sin(theta + dtheta/2)     0.5*cos(theta+dtheta/2)+(ds/(2*b))*sin(theta + dtheta/2);
    0.5*sin(theta+dtheta/2)+(ds/(2*b))*cos(theta + dtheta/2)     0.5*sin(theta+dtheta/2)-(ds/(2*b))*cos(theta + dtheta/2);
    1/b     -1/b];


x = x + dx;
y = y + dy;
theta = theta + dtheta;

Cp = d_position*Cp*d_position' + d_wheel*Covariance_Matrix*d_wheel';
end
