% ----------- Flywheel_stresses_code ---------------------------------

function [thk_sol, rad_sol, sigma_r, sigma_t, sigma_vonMises] = Flywheel_stresses(point_tArr,ll, h, ul, ...
    C_Func, D_Func, E_Func, F_Func, thk_Func, rad_Func, Drad_Func, rho, omega)
%testing...
%point_tArr=[0.0189 0.01 0.01 0.01 0.01 0.01 0.0246 0.06];

%solution using Finite difference method
u_seg=ll:h:ul;

for j=1:ul
    st=(1/h)*(j-1)+1;en=(1/h)*j;
    C_sol(st:en)=C_Func{j}(point_tArr,u_seg(st:en));
    D_sol(st:en)=D_Func{j}(point_tArr,u_seg(st:en));
    E_sol(st:en)=E_Func{j}(point_tArr,u_seg(st:en));
    F_sol(st:en)=F_Func{j}(point_tArr,u_seg(st:en));
    thk_sol(st:en)=thk_Func{j}(point_tArr,u_seg(st:en));
    rad_sol(st:en)=rad_Func{j}(u_seg(st:en));
    Drad_sol(st:en)=Drad_Func{j}(u_seg(st:en));    
end

u_seg_size=length(u_seg);
C_sol(u_seg_size)=C_Func{j}(point_tArr,u_seg(u_seg_size));
D_sol(u_seg_size)=D_Func{j}(point_tArr,u_seg(u_seg_size));
E_sol(u_seg_size)=E_Func{j}(point_tArr,u_seg(u_seg_size));
F_sol(u_seg_size)=F_Func{j}(point_tArr,u_seg(u_seg_size));
thk_sol(u_seg_size)=thk_Func{j}(point_tArr,u_seg(u_seg_size));
rad_sol(u_seg_size)=rad_Func{j}(u_seg(u_seg_size));
Drad_sol(u_seg_size)=Drad_Func{j}(u_seg(u_seg_size));
    
flag1 = 1; p1 = 0; flag2 = 1; p2 = 0;
z_sol = twopoint(u_seg,C_sol,D_sol,E_sol,F_sol,flag1,flag2,p1,p2)'; %to make column vector
for i=2:u_seg_size-1
    Dz(i)=(z_sol(i+1)-z_sol(i-1))./(2*h);%using central difference approximation for Dz at intermediate points
end

%using forward difference approximation for Dz at end points
Dz(1)=(z_sol(2)-z_sol(1))/h;
Dz(u_seg_size)=(z_sol(u_seg_size)-z_sol(u_seg_size-1))/h;

sigma_r=z_sol./(thk_sol.*rad_sol);
sigma_t=(1./thk_sol).*(Dz./Drad_sol+rho*omega^2*rad_sol.^2.*thk_sol);
sigma_vonMises=(sigma_r.^2+sigma_t.^2-sigma_r.*sigma_t).^(0.5);
end
