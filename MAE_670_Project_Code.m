% MAE670 NonLinear Control Class Project
% Instructor Professor Dr. TarnRaj Singh
% Submitted by Swapneel Bhavesh Mehta and Ankur Soni
% University at Buffalo

clc;
clear all;
close all;

tspan=[0, 20];
x_0 = [0 10 pi/4 0 0 0 0];
options=odeset('RelTol',1e-8,'AbsTol',1e-8);
[t2, x] = ode45(@mobile_robot, tspan, x_0, options);

q_0=[0,10,pi/4,0,0,0,0,0,0,0];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t1,q]=ode45(@dynamics,tspan,q_0,options);

figure(1)
plot(t1,q(:,1),t2,x(:,1))
legend('Uncontrolled','Controlled')
xlabel('Time')
ylabel('Position x')

figure(2)
plot(t1,q(:,2),t2,x(:,2))
legend('Uncontrolled','Controlled')
xlabel('Time')
ylabel('Position y')

figure(3)
plot(t1,q(:,3),t2,x(:,3))
legend('Uncontrolled','Controlled')
xlabel('Time')
ylabel('Orientation angle of vehicle (phi)')

figure(4)
plot(t1,q(:,4),t2,x(:,4))
legend('Uncontrolled','Controlled')
xlabel('Time')
ylabel('Anguar position of right wheel')

figure(5)
plot(t1,q(:,5),t2,x(:,5))
legend('Uncontrolled','Controlled')
xlabel('Time')
ylabel('Anguar position of left wheel')

figure(6)
plot(t1,q(:,9),t2,x(:,6))
legend('Uncontrolled','Controlled')
xlabel('Time')
ylabel('Anguar veloctiy of right wheel')

figure(7)
plot(t1,q(:,10),t2,x(:,7))
legend('Uncontrolled','Controlled')
xlabel('Time')
ylabel('Anguar velocity of left wheel')

function [X_dot] = mobile_robot(t, x)
l=0.5; % l=r/2b
d_1=0.75; d_2=1;
Alpha=2; Beta=5; % Feedback gains 
Xref_g=5; Yref_g=0; % Referenece desired point

% States defined 
Xg = x(1);  
Yg = x(2);
phi = x(3);
theta_r     = x(4); 
theta_l     = x(5); 
theta_r_dot = x(6); 
theta_l_dot = x(7); 

 % Defining desired path for robot
y_goal = t;      
x_goal = t;
traj_goal =[x_goal y_goal]';
traj_dot_goal = [1; 1]; 
traj_2_dot_goal = [0; 0]; 

Xref = Xg+Xref_g*cos(phi)-Yref_g*sin(phi);
Yref = Yg+Xref_g*sin(phi)+Yref_g*cos(phi);
ref_point = [Xref;Yref];

S = [l*(d_1*cos(phi)-d_2*sin(phi)) l*(d_1*cos(phi)+d_2*sin(phi)); l*(d_1*sin(phi)+d_2*cos(phi)) l*(d_1*sin(phi)-d_2*cos(phi));l -l;1 0;0 1];
 
ang_rate_wheels = [theta_r_dot;theta_l_dot];
cap_phi = [ l*(d_1*cos(phi)-d_2*sin(phi))+(-Xref_g*sin(phi)-Yref_g*cos(phi))*l,l*(d_1*cos(phi)+d_2*sin(phi))-(-Xref_g*sin(phi)-Yref_g*cos(phi))*l;
                l*(d_1*sin(phi)+d_2*cos(phi))+(Xref_g*cos(phi)-Yref_g*sin(phi))*l,l*(d_1*sin(phi)-d_2*cos(phi))-(Xref_g*cos(phi)-Yref_g*sin(phi))*l ]; 

cap_phi_dot =[ l*(-d_1*sin(phi)-d_2*cos(phi))+(-Xref_g*cos(phi)+Yref_g*sin(phi))*l,l*(-d_1*sin(phi)+d_2*cos(phi))-(-Xref_g*cos(phi)+Yref_g*sin(phi))*l;
             l*(d_1*cos(phi)-d_2*sin(phi))+(-Xref_g*sin(phi)-Yref_g*cos(phi))*l,l*(d_1*cos(phi)+d_2*sin(phi))-(-Xref_g*sin(phi)-Yref_g*cos(phi))*l ];
               
xy_dot  = cap_phi*ang_rate_wheels;
xy2_dot = traj_2_dot_goal+Beta*(traj_dot_goal - xy_dot)+Alpha*(traj_goal - ref_point);
u = (cap_phi)\(xy2_dot-cap_phi_dot*ang_rate_wheels);
q_dot = S*ang_rate_wheels;
X_dot = [q_dot(1);q_dot(2);q_dot(3);q_dot(4);q_dot(5); 0; 0] + [0 0;0 0;0 0;0 0;0 0;1 0;0 1]*[u(1);u(2)];
end

function q_dotdot= dynamics(~,q)

m=3.5;
m_c=2.5;
d=1.2;
I_1=12;
I_2=1.5;
TAU1=1.5; % Some constant values for torque input, No control
TAU2=2.5;
lambda1=1;
lambda2=1;
lambda3=1;
b=1;
r=0.5;

x=q(1);
y=q(4);
phi=q(3);
thetar=q(4);
thetal=q(5);
xdot=q(6);
ydot=q(7);
phidot=q(8);
thetardot=q(9);
thetaldot=q(10);

MI=[m                  0            -m_c*d*sin(phi) 0  0;
   0                 m              m_c*d*cos(phi) 0  0;
   -m_c*d*sin(phi)   m_c*d*cos(phi)  I_1             0  0;
   0                 0               0           I_2  0;
   0                 0               0           0  I_2];

V=[-m_c*d*((phidot)^2)*cos(phi);
    -m_c*d*((phidot)^2)*sin(phi);
               0;
               0;
               0];
 
 E=[0 0;
    0 0;
    0 0;
    1 0;
    0 1];

A=[-sin(phi) cos(phi) -d 0 0;
    -cos(phi) -sin(phi) -b r 0;
    -cos(phi) -sin(phi) b 0 r];

lambda=[lambda1;lambda2;lambda3];
TAU=[TAU1;TAU2];

vec=MI\(-V+E*TAU-((A)'*lambda));

x_dotdot=vec(1);
y_dotdot=vec(2);
phi_dotdot=vec(3);
thetar_dotdot=vec(4);
thetal_dotdot=vec(5);
q_dotdot=[xdot;ydot;phidot;thetardot;thetaldot;x_dotdot;y_dotdot;phi_dotdot;thetar_dotdot;thetal_dotdot];
end