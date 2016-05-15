function result = trans()
%改变matlab中的局部基函数的存储形式以便适应fortran程序的需求

%% modify the data
%=======================================

% read the data from the old files and full it.
load basis_data_MsFEM.mat  loc_basis
loc_basis = full(loc_basis);
[m,n] = size(loc_basis);
fprintf('%d,%d\n',m,n);
%=======================================
%% build the target files
%========================================

fopen('nbase1','w+');
fopen('nbase2','w+');
fopen('nbase3','w+');
fopen('nbase4','w+');

%========================================
%% modify the form of the base_function information
%========================================


%========================================
