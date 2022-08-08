%-------------------------------------------------------------------------%
%               UFSC - Federal University of Santa Catarina               %
%               Graduate Program in Mechanical Engineering                %
%                                                                         %
%     Programmer:                                                         %
%       Davi Klein                                                        %
%                                                                         %
%   Version: 1.0                                              08/09/2022  %
%=========================================================================%
%                          Program Descriprion                            %
%=========================================================================%
%	Main file, responsible for all the routines regarding the simulation  %
%   of a SPS Stewart-Gough Platform                                       %
%-------------------------------------------------------------------------%
clear
close all
clc
format short eng
%-------------------------------------------------------------------------%
% Parametros de simulacao %
%-------------------------------------------------------------------------%
step = 0.01;
Ts = 2.5;
N = 14.5370e+003;%round(Ts/step); 
w = 3.0;
to = 0;
syms to;
global e
e = zeros(6,1);
global aux
aux = 0;
%-------------------------------------------------------------------------%
%                        Structural Parameters                            %
%-------------------------------------------------------------------------%

e1 = 0.5; % distance between joint Ai and the cilinder center of mass (m)
e2 = 0.5; % distance between joint Ai and the piston center of mass (m)
m1 = [0.1,0.1,0.1,0.1,0.1,0.1]; % cilinder mass (Kg)
m2 = [0.1,0.1,0.1,0.1,0.1,0.1]; % piston mass (Kg)
I1ii = zeros(3,3,6); % cilinder moment of inertia (kg.m^2)
I2ii = zeros(3,3,6); % piston moment of inertia (kg.m^2)
for i = 1:6
    I1ii(:,:,i) = [0.00625 0 0; 0 0.00625 0; 0 0 0];
    I2ii(:,:,i) = [0.00625 0 0; 0 0.00625 0; 0 0 0];
end

mp = 1.5; % moving platform mass (Kg)
Ipb = [0.08 0 0; 0 0.08 0; 0 0 0.08]; %moving platform moment of inertia (kg.m^2)

% fixed platform joint coordinates (ai) for a generic configuration
a = zeros(3,6); 
a_abs = 2; 
alfa_o = 0; % angle between the first joint and the x axis of the reference system (°)
aj = 60; % angle between joints (°)
alfa = [alfa_o alfa_o+aj alfa_o+120 alfa_o+120+aj alfa_o-120 alfa_o-120+aj];
for i = 1:6
    alfa(i) = deg2rad(alfa(i));
end

for i = 1:3
    for j = 1:6
        if i == 1
            a(i,j) = a_abs*cos(alfa(1,j));
        elseif i == 2
            a(i,j) = a_abs*sin(alfa(1,j));
        else
            a(i,j) = 0;
        end
    end
end
a = [-2.120 -2.380 -2.380 -2.120 0 0; 1.374 1.224 -1.224...
     -1.374 -0.15 0.15; 0 0 0 0 0 0]; %Tsai parameters
 
% Moving platform joint coordinates (bi)

% u = 0.5; % scaling factor between joints ai and bi
% b = u*a;

b = [0.170 -0.600 -0.600 0.170 0.430 0.430; 0.595 0.150... 
     -0.150 -0.595 -0.445 0.445; -0.4 -0.4 -0.4 -0.4 -0.4 -0.4]; %Tsai parameters

%-------------------------------------------------------------------------%
%                         Simulation Routines                             %
%-------------------------------------------------------------------------%

L = zeros(1,6);
s = zeros(3,6,N);
sx = zeros(3,3,6,N);
su = zeros(3,6,N);

P = [-1.5+(0.2*sin(w*to)); 0.2*sin(w*to);1.0+(0.2*sin(w*to))];
dP = diff(P);
ddP = diff(dP);
sdP = size(dP);
sddP = size(ddP);

if sdP(1) < 3
    dP = zeros(3,N);
end

if sddP(1) < 3
    ddP = zeros(3,N);
end

if sddP(2) < N
    for i = 1:N
        ddP(:,i) = ddP(:,1);
    end
end

ea = zeros(3,N);
dea = diff(ea);
ddea = diff(dea);
sdea = size(dea);
sddea = size(ddea);

if sdea(1) < 3
    dea = zeros(3,N);
end

if sddea(1) < 3
    ddea = zeros(3,N);
end

to = linspace(0,2.0944,N);

P = subs(P);
dP = subs(dP);
ddP = subs(ddP);
ea = subs(ea);
dea = subs(dea);
ddea = subs(ddea);

X = [P;ea];
dX = [dP;dea];

fe = zeros(3,1);
ne = zeros(3,1);

RBA = zeros(3,3,N);

tau = zeros(6,N);
taune = zeros(6,N);

RCA = zeros(3,6,6,N);

wt = w*to;

for i = 1:N
    RBA(:,:,i) = [cos(ea(3,i))*cos(ea(2,i)) cos(ea(3,i))*sin(ea(2,i))*sin(ea(1,i))-sin(ea(3,i))*cos(ea(1,i)) cos(ea(3,i))*sin(ea(2,1))*cos(ea(1,1))+sin(ea(3,1)*sin(ea(1,1)));...
                  sin(ea(3,i))*cos(ea(2,i)) sin(ea(3,i))*sin(ea(2,i))*sin(ea(1,i))+cos(ea(3,i))*cos(ea(1,i)) sin(ea(3,i))*sin(ea(2,1))*cos(ea(1,1))-cos(ea(3,1))*sin(ea(1,1));...
                 -sin(ea(2,i)) cos(ea(2,i))*sin(ea(1,i)) cos(ea(2,i))*cos(ea(1,i))];
    
    % Inverse Kinematics %
    [L(i,:),s(:,:,i)] = Inv_Kin(P(:,i),a,b,RBA(:,:,i));
    
    % Inverse Dynamics %
    [tau(:,i)] = Din_Inv_v2(L(i,:),s(:,:,i),RBA(:,:,i),dP(:,i),ddP(:,i),dea(:,i),ddea(:,i),b,e1,e2,m1,m2,I1ii,I2ii,mp,Ipb,fe,ne);
end

% ode system initial conditions
X0 = [-1.3 ; 0.2 ; 1.2 ; 0 ; 0 ; 0] ; 
dX0 = [0 ; 0 ; 0 ; 0 ; 0 ; 0] ;

% Decentralized PD Control %
[t,X_new] = ode45(@(t,X_new) Control(t,X_new,a,b,I1ii,I2ii,m1,m2,e1,e2,mp,Ipb), [0 2.0944], [X0, dX0],odeset('OutputFcn','odeplot','OutputSel',[3]));
te = X(:,1:3)' - X_new(1:3,:);
re = X(:,4:6)' - X_new(4:6,:);
%-------------------------------------------------------------------------%
%                           Results Plottings                             %
%-------------------------------------------------------------------------%

%------------------------- Inverse Kinematics ----------------------------%

l1 = L(:,1);
l2 = L(:,2);
l3 = L(:,3);
l4 = L(:,4);
l5 = L(:,5);
l6 = L(:,6);

figure(1)
plot(wt,l1,wt,l2,wt,l3,wt,l4,wt,l5,wt,l6)
grid on
title('Actuator Length (m) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('Actuator Length (m)')
legend('l1','l2','l3','l4','l5','l6')

%-------------------------- Inverse Dynamics -----------------------------%

t1 = tau(1,:);
t2 = tau(2,:);
t3 = tau(3,:);
t4 = tau(4,:);
t5 = tau(5,:);
t6 = tau(6,:);

figure(2)
plot(wt,t1,wt,t2,wt,t3,wt,t4,wt,t5,wt,t6)
title('Input Forces (N) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('Input Forces (N)')
legend('\tau_1','\tau_2','\tau_3','\tau_4','\tau_5','\tau_6')
grid on

%----------------------- Decentralized PD Control -------------------------%

figure(3)
plot(w*t,X_new(:,1),wt,P(1,:))
grid on
title('X axis translation (m) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('X axis translation (m)')
legend('Executed trajectory', 'Reference trajectory')

figure(4)
plot(w*t,X_new(:,2),wt,P(2,:))
grid on
title('Y axis translation (m) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('Y axis translation (m)')
legend('Executed trajectory', 'Reference trajectory')

figure(5)
plot(w*t,X_new(:,3),wt,P(3,:))
grid on
title('Z axis translation (m) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('Z axis translation (m)')
legend('Executed trajectory', 'Reference trajectory')

figure(6)
plot(w*t,X_new(:,4),wt,ea(1,:))
grid on
title('X axis rotation (rad) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('X axis rotation (rad)')
legend('Executed trajectory', 'Reference trajectory')

figure(7)
plot(w*t,X_new(:,5),wt,ea(2,:))
grid on
title('Y axis rotation (rad) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('Y axis rotation (rad)')
legend('Executed trajectory', 'Reference trajectory')

figure(8)
plot(w*t,X_new(:,6),wt,ea(3,:))
grid on
title('Z axis rotation (rad) vs. Dimensionless Time (\omega t)')
xlabel('Dimensionless Time (\omega t)')
ylabel('Z axis rotation (rad)')
legend('Executed trajectory', 'Reference trajectory')

figure(9)
plot(w*t,te(1,:),w*t,te(2,:),w*t,te(3,:))
grid on
title('Translational tracking error')
xlabel('Dimensionless Time (\omega t)')
ylabel('Tracking error')

figure(10)
plot(w*t,re(1,:),w*t,re(2,:),w*t,re(3,:))
grid on
title('Rotational tracking error')
xlabel('Dimensionless Time (\omega t)')
ylabel('Tracking error')