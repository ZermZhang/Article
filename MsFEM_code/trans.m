%改变matlab中的局部基函数的存储形式以便适应fortran程序的需求

function result = trans()
%=======================================
%% initialize
%=======================================

% read the data from the old files and full it.
load basis_data_MsFEM.mat  loc_basis
%loc_basis = full(loc_basis);      it's too big to full.
[m,n] = size(loc_basis);
fprintf('m=%d,n=%d\n',m,n);
Nx = sqrt(n) - 1;
Ny = Nx;
nx = (sqrt(m) - 1)/Nx;
ny=nx;
fprintf('Nx=%d,nx=%d\n',Nx,nx)

%=======================================
%% build the target files and init the variables
%========================================

% open the files contains the result
fp1=fopen('nbase1','w+');
fp2=fopen('nbase2','w+');
fp3=fopen('nbase3','w+');
fp4=fopen('nbase4','w+');

% init the varables needed
tmp_basis1=zeros(1,(nx+1)*(ny+1));
tmp_basis2=zeros(1,(nx+1)*(ny+1));
tmp_basis3=zeros(1,(nx+1)*(ny+1));
tmp_basis4=zeros(1,(nx+1)*(ny+1));

% init the index information
idx = reshape(1:(Ny*ny+1)*(Nx*nx+1),Ny*ny+1,Nx*nx+1);

%========================================
%% modify the form of the base_function information
%========================================

for i = 1:Nx
       for j = 1:Ny
           loc_idx = idx((j-1)*(ny)+(1:ny+1),(i-1)*(nx)+(1:nx+1));
           loc=reshape(loc_idx,(nx+1)*(ny+1),1);
           tmp_basis1=loc_basis(loc,(i-1)*(Ny+1)+j     );
           tmp_basis3=loc_basis(loc,(i-1)*(Ny+1)+j+1);
           tmp_basis2=loc_basis(loc,(i   )*(Ny+1)+j     );
           tmp_basis4=loc_basis(loc,(i   )*(Ny+1)+j+1);
           
           fprintf(fp1,'%d ',(Nx+1)*(Ny+1));
           fprintf(fp2,'%d ',(Nx+1)*(Ny+1));
           fprintf(fp3,'%d ',(Nx+1)*(Ny+1));
           fprintf(fp4,'%d ',(Nx+1)*(Ny+1));
           for m = 1:nx+1
               for l = 1:ny+1
                   fprintf(fp1,'%f ',full(tmp_basis1(((l-1)*(ny+1)+m),1)));
                   fprintf(fp2,'%f ',full(tmp_basis2(((l-1)*(ny+1)+m),1)));
                   fprintf(fp3,'%f ',full(tmp_basis3(((l-1)*(ny+1)+m),1)));
                   fprintf(fp4,'%f ',full(tmp_basis4(((l-1)*(ny+1)+m),1)));
               end
           end
           fprintf(fp1,'\n');
           fprintf(fp2,'\n');
           fprintf(fp3,'\n');
           fprintf(fp4,'\n');
       end
end

%========================================
%% close the files
%========================================
fclose(fp1');
fclose(fp2);
fclose(fp3);
fclose(fp4);
%========================================