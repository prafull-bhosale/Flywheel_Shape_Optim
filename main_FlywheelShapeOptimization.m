% ----------------------- FlywheelShapeOptimization_code --------------------------

clearvars; clc; close all;
syms u z(u) t1 t2 t3 t4 t5 t6 t7 t8
global n k U
n=8; %number of control points
k=4; %for cubic B spline
S=n-k+1;
h=0.01; %step size for FDM
ll=0;
ul=S;
rho=7250;
nu=.3;
omega=65.45;
Mmax=115;
sigma_allow=6.4e6;
UArr=0.5:1:S;%[0.5 1.5 2.5 3.5 4.5];

R1=0.06;
R2=0.5;
DeltaR=(R2-R1)/(n-1);
point_r=R1:DeltaR:R2;
%point_r = sym('r',[1;n]);
point_t = sym('t',[1;n]);

%B spline implementation starts
for i=1:n+k
    x_i(i)=x(i-1);
end
x_i;

for j=1:S
    U=UArr(j);
    radj=0;
    thkj=0;
    for i=0:n-1
        N_ik=simplify(N(i,k));
        radj=radj+N_ik*point_r(i+1);
        thkj=thkj+N_ik*point_t(i+1);
    end
    rad(j)=radj;
    thk(j)=thkj;
end
rad
thk
Drad=diff(rad,u);
D2rad=diff(rad,u,2);
Dthk=diff(thk,u);
for j=1:S
    thk_Func{j}=matlabFunction(thk(j),'Vars',{point_t,u});
    rad_Func{j}=matlabFunction(rad(j));
    Drad_Func{j}=matlabFunction(Drad(j));
    D2rad_Func{j}=matlabFunction(D2rad(j));
    Dthk_Func{j}=matlabFunction(Dthk(j),'Vars',{point_t,u});
end

disp('---------------------------------------');

%Volume, Mass KE calculations starts
Ek=0;
M=0;
cnt1=1;
for j=1:S    
    Vk=2*pi*int(thk(j)*rad(j)*Drad(j),u,j-1,j); % volume of kth segment
    mk=rho*Vk; % mass of kth segment
    ek=pi*rho*omega^2*int(thk(j)*rad(j)^3*Drad(j),u,j-1,j); % KE of kth segment
    M=M+mk;
    Ek=Ek+ek;
    
    %stress calculations
    %requirements for FDM starts
    C(j)=(rad(j)^2)*Drad(j);
    D(j)=rad(j)*(Drad(j)^2)-(rad(j)^2)*D2rad(j)-((rad(j)^2)/thk(j))*Drad(j)*Dthk(j);
    E(j)=nu*(rad(j)/thk(j))*Dthk(j)*(Drad(j)^2)-(Drad(j)^3);
    F(j)=-((3+nu)*rho*omega^2*thk(j)*(rad(j)^3)*(Drad(j)^3));
    %requirements for FDM ends
end

M
Ek
mass=matlabFunction(M,'Vars',{point_t});
KE=matlabFunction(Ek,'Vars',{point_t});
C=simplify(C);
D=simplify(D);
E=simplify(E);
F=simplify(F);
%requirements for FDM starts
for j=1:S
    C_Func{j}=matlabFunction(C(j),'Vars',{point_t,u});
    D_Func{j}=matlabFunction(D(j),'Vars',{point_t,u});
    E_Func{j}=matlabFunction(E(j),'Vars',{point_t,u});
    F_Func{j}=matlabFunction(F(j),'Vars',{point_t,u});
end

%results from research paper
%x_rp=0.02*ones(1,n)%Original Varshney
x_rp=[0.0189 0.01 0.01 0.01 0.01 0.01 0.0246 0.06];% Result from research paper obtained using Jaya algo
disp('From Research paper');
fprintf('Mass = %f',mass(x_rp));
fprintf('Kinetic Energy = %f', KE(x_rp));

disp('---- JAYA algorithm starts ----------');
pop = 1000;               % Population size
var = 8;                 % Number of design variables
maxGen = 1000;            % Maximum number of iterations
mini = 0.010*ones(1,var);  % Lower Bound of Variables
maxi = 0.060*ones(1,var);   % Upper Bound of Variables
%objective = @myobj;      % Cost Function
FunctionTolerance=1e-6;
MaxStallGenerations=50;
penaltyParam=1e8;

[xopt,fopt] = JAYA_Algorithm(@myobj,pop,var,[],maxGen,mini,maxi,FunctionTolerance,MaxStallGenerations,{mass KE rho omega Mmax sigma_allow ll h ul C_Func D_Func E_Func F_Func thk_Func rad_Func Drad_Func penaltyParam});
fprintf('Mass = %f',mass(xopt));
fprintf('Kinetic Energy = %f',KE(xopt));

[thk_sol, rad_sol, sigma_r, sigma_t, sigma_vonMises] = Flywheel_stresses(xopt,ll, h, ul, C_Func, D_Func, E_Func, F_Func, thk_Func, rad_Func, Drad_Func, rho, omega);
%convert m to mm
point_r=point_r.*1e3;
xopt=xopt.*1e3;
rad_sol=rad_sol.*1e3;
thk_sol=thk_sol.*1e3;

figure(2)
plot(rad_sol,sigma_r,'DisplayName','Radial Stress','LineWidth',2);
hold on
plot(rad_sol,sigma_t,'DisplayName','Tangential Stress','LineWidth',2);
plot(rad_sol,sigma_vonMises,'DisplayName','Von-Mises Stress','LineWidth',2);
xlabel('Radius (mm)')
ylabel('Stresses (Pa)');
legend
hold off

figure(3)
plot(rad_sol,thk_sol,'b','DisplayName','Optimized shape of the flywheel','LineWidth',2);
hold on
scatter([point_r point_r],[xopt -xopt],[],[0 0 0],'filled','DisplayName','Control Points');
hold on
plot(rad_sol,-thk_sol,'b','LineWidth',2);
hold off
xlabel('Radius (mm)')
ylabel('Thickness (mm)');
legend

function [f]=myobj(x,sp_params)
rho=sp_params{3};
omega=sp_params{4};
Mmax=sp_params{5};
sigma_allow=sp_params{6};
ll=sp_params{7};
h=sp_params{8};
ul=sp_params{9};
C_Func=sp_params{10};
D_Func=sp_params{11};
E_Func=sp_params{12};
F_Func=sp_params{13};
thk_Func=sp_params{14};
rad_Func=sp_params{15};
Drad_Func=sp_params{16};
penaltyParam=sp_params{17};

[r,c]=size(x);
Z=zeros(r,1);
mass=sp_params{1};
KE=sp_params{2};
for i=1:r   
    z=-KE(x(i,:));
    g1=mass(x(i,:))-Mmax;
    [~, ~, ~, ~, sigma_vonMises] = Flywheel_stresses(x(i,:),ll, h, ul, C_Func, D_Func, E_Func, F_Func, thk_Func, rad_Func, Drad_Func, rho, omega);
    g2=max(sigma_vonMises)-sigma_allow;
    if(g1>0)
        p1=penaltyParam;
    else
        p1=0;
    end
    if(g2>0)
        p2=penaltyParam;
    else
        p2=0;
    end
%     p1=penaltyParam*((max(0,g1))^2); % penalty if constraint 1 is violated
%     p2=penaltyParam*((min(0,g2))^2); % penalty if constraint 2 is violated
    Z(i)=z+p1+p2; % penalized objective function value
end
f=Z;
end

function x=x(i)
global n k
%if 0<=i && i<=n+k
if i<k
    x=0;
elseif k<=i && i<=n
    x=i-k+1;
else
    x=n-k+1;
end
% else
%     t=0;
% end
end

function N_ik=N(i,k)
global U
syms u
if k==1
    if x(i)<=U && U<=x(i+1)
        N_ik=1;
    else
        N_ik=0;
    end
else
    %disp('N start:------------');
    xi_k=x(i);
    Ni_k_1=N(i,k-1);
    Ni_k_1(isinf(Ni_k_1))=0;
    
    xipk_1_k=x(i+k-1);
    xipk_k=x(i+k);
    Nip1_k_1=N(i+1,k-1);
    Nip1_k_1(isinf(Nip1_k_1))=0;
    
    xip1_k=x(i+1);
    
    num1=1/(xipk_1_k-xi_k);
    num1(isinf(num1))=0;
    num2=1/(xipk_k-xip1_k);
    num2(isinf(num2))=0;
    
    N_ik=Ni_k_1*(u-xi_k)*num1+Nip1_k_1*(xipk_k-u)*num2;
    %disp('N end:-----------');
end
end


