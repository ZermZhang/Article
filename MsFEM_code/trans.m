function result = trans()
%�ı�matlab�еľֲ��������Ĵ洢��ʽ�Ա���Ӧfortran���������

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
