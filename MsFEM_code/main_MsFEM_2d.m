% This program use SIPDG method to solve the elliptic equation 
%     
%      - \nabla (a \nabla u) = f
%     
%
%=========================================================
%% initialize
%=========================================================
% clear all
format long
%=========================================================
%% input parameter
%=========================================================

plot_solution_fine = 1;
plot_solution_MS = 1;
compute_MS_basis = 1;
solve_fine = 1;
compute_error = 1;

% Input the num of Corase block , num of fine block in one Corase block 
Nx=32; Ny=Nx;       %num of corase block 
nx=32; ny=nx;       %num of fine block in one corase block


%=========================================================
%% input medium
%=========================================================
[X1,Y1] = meshgrid(1/nx/Nx/2:1/nx/Nx:1,1/ny/Ny/2:1/ny/Ny:1);


epsilon = 1/100;
% source_a = @(x,y) (1.1 + sin(pi*x/epsilon  ).*sin(pi*(1+y/epsilon)) + sin(pi*(1+(x+y)/epsilon)  ).^2);

%source_a = @(x,y) (4 + sin(pi*(sqrt(y/epsilon)-x/epsilon)))./(2 + sin(pi*(y/epsilon-2*x/epsilon)));
% source_a = @(x,y) (5 + sin(pi*(sqrt(y/epsilon)-x/epsilon)))./(2 + sin(pi*(y-x.^2).*(y/epsilon-2*x/epsilon)));
source_a = @(x,y) (2 + 1.8*sin(2*pi*x/epsilon))./(2+1.8*cos(2*pi*y/epsilon))...
    +(2+1.8*sin(2*pi*y/epsilon))./(2+1.8*sin(2*pi*x/epsilon));

%source_a = @(x,y) (1+0*x);

a = source_a(X1,Y1);


%=========================================================
%% forming matrix
%=========================================================

disp('Forming IPDG matrix')
% Global_DA is the fine-scale stiffness martix 
% Global_M  is the fine-scale mass martix
% boundary  is the boundary index

[Global_DA,Global_M,boundary] = finematrix_2d(a,nx,ny,Nx,Ny);

%=========================================================
%% forming RHS
%=========================================================
% input a source function (fun_F) in this part (or you can directly input a vector (F))

disp('Forming RHS')
%fun_F = @(x,y) 2*pi^2*sin(pi*x).*sin(pi*y);
fun_F = @(x,y) -1 + 0*x;

[X,Y] = meshgrid(1/nx/Nx/2:1/nx/Nx:1,1/ny/Ny/2:1/ny/Ny:1);
F = fun_F(X,Y);
f = form_Source(F,nx,ny,Nx,Ny);
%f = -ones(Ny*ny+1,Nx*nx+1);

%=========================================================
%% forming MS basis
%=========================================================
if compute_MS_basis
    disp('Forming MS basis')
    [loc_basis] = MsFEM_2d_basis(Global_DA,nx,ny,Nx,Ny);
%     [loc_basis] = bilinear_partition(nx,ny,Nx,Ny);
    save basis_data_MsFEM loc_basis
else
    
    load basis_data_MsFEM loc_basis
    
end
%=========================================================
%% forming MS Matrix
%=========================================================
disp('Forming MS Matrix')

interior_idx_fine = 1:(nx*Nx+1)*(ny*Ny+1);
interior_idx_fine(boundary) = [];

interior_idx_coarse = 1:(Nx+1)*(Ny+1);
interior_idx_coarse([1:Ny+1:(Nx+1)*(Ny+1), Ny+1:Ny+1:(Nx+1)*(Ny+1),  2:Ny  ,(Ny+1)*Nx + (2:Ny)]) = [];

Global_DA = Global_DA(interior_idx_fine,interior_idx_fine);
Global_M = Global_M(interior_idx_fine,interior_idx_fine);
f(boundary) = [];
loc_basis = loc_basis(interior_idx_fine,interior_idx_coarse);

MS_A = loc_basis'*Global_DA*loc_basis;
MS_f = loc_basis'*f;

%=========================================================
%% solving solution
%=========================================================
if solve_fine
    disp('Solving fine solution')

    Global_U = Global_DA\f;
end

disp('Solving MS solution')

MS_CU = MS_A\MS_f;
MS_U = loc_basis*MS_CU;

%=========================================================
%% visualization
%=========================================================
MS_sol = reshape(MS_U,ny*Ny-1,nx*Nx-1);
Global_sol = reshape(Global_U,ny*Ny-1,nx*Nx-1);
%=========================================================
%% compute error
%=========================================================
if compute_error
DG_error = sqrt(((MS_U(:)-Global_U(:))'*Global_DA*(MS_U(:)-Global_U(:)))/((Global_U(:))'*Global_DA*(Global_U(:))));
L2_error = sqrt(((MS_U(:)-Global_U(:))'*Global_M*(MS_U(:)-Global_U(:)))/((Global_U(:))'*Global_M*(Global_U(:))));

fprintf('The relative Energy Error for MsFEM is  %2.8f. \n',DG_error)
fprintf('The relative L2 Error for MsFEM is      %2.8f. \n',L2_error)
end
%=========================================================
