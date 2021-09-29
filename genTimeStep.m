%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate Time step suitable for corresponding material
%           properties
%           input: material properties
%           output: time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear
format long;
%% Input parameters
global E Rho nu v r  m
E = 1e8; % Young's Modules
Rho = 3.2; % Density
nu = 0.30; % Poisson's ratio
v = 100; % maximum velocity
r = 0.3; % average radius
m = 4/3*pi*r^3*Rho; % particles mass
%% Estimate time step, model 1
T = 5.84*(Rho*(1-nu^2)/E)^(2/5)*r*v^(-1/5);
dT = T/20;
fprintf('模型一估计时间步长为  %d\n', dT);
%% Estimate time step, model 2
T = 2.94*(15*m/16/E/sqrt(r*v))^(0.4);
dT = T/20;
fprintf('模型二估计时间步长为  %d\n', dT);
%% Compare to Linear Model
Kn = 2.0*E*sqrt(r*r/2)/(3*(1-nu^2));
fprintf('等效线性刚度为  %d\n', Kn);

MaxAcc = 2.0*E*sqrt(r*r)/(3*(1-nu^2))*r/m;
fprintf('最大容许加速度为  %d\n', Kn);

%% Numerical simulation
[Time,Pos]=ode45('contact_model',[0 T*2],[0 v]);
plotyy(Time,real(Pos(:,1)),Time,real(Pos(:,2)));

%% Result check
MaxOverlap = max(real(Pos(:,1)))/r;
fprintf('最大重叠量  %.4f%%\n', MaxOverlap*100);