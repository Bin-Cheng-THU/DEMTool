function [Vol,Den,N_p,N_f,N_e,p_data,f_data,e_data,len_data,ne_data,nf_data]=loadpoly(c)
%
%specified type: triangle facet polyhydron
%vectors all represented in body coordinates, default unit: m
disp('Start preloading...');
fid = fopen(c);%unit: m
tmp = fscanf(fid,'%f',[1,2]);
Vol = tmp(1); % m^3
Den = tmp(2); % kg/m^3
tmp = fscanf(fid,'%f %f %f',[3,1]);
tmp = fscanf(fid,'%d %d %d',[3,1]);
N_p = tmp(1);
N_f = tmp(2);
N_e = tmp(3);
p_data = fscanf(fid,'%f %f %f',[3,N_p]);
p_data = [(1:N_p).',p_data'];
tmp = fscanf(fid,'%d %d %d %f %f %f',[6,N_f]);
tmp = tmp.';
f_data = [(1:N_f).',tmp(:,1:3)+1];
nf_data = [(1:N_f)',tmp(:,4:6)];
tmp = fscanf(fid,'%d %d %d %d %f %f %f %f %f %f %f',[11,N_e]);
tmp = tmp.';
e_data = [(1:N_e).',tmp(:,1:4)+1];
len_data = [(1:N_e)',tmp(:,5)];
ne_data = [(1:N_e)',e_data(:,4),tmp(:,6:8),e_data(:,5),tmp(:,9:11)];
fclose(fid);
disp('Preloading completed.')
end