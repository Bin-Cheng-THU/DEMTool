% Analysis of Surface Motion by using Spherical Harmonics
clc
clear
% basic parameters --------------------------------------------------------
G = 6.674184e-11; % G constant, unit: N(m/kg)^2 or m^3/kg/s^2
Period = 0;%2.2593*3600; % Autorotation period, unit: s
Radius = 450;%Equatorial radius; unit: m
% options -----------------------------------------------------------------
flag_genSH = 1;         % calculate the Spherical Harmonics
flag_loadSH = 1;        % load the Spherical Harmonics
flag_SHshape = 1;       % show the shape of Spherical Harmonics
flag_triSH = 1; 
%----------------------------------------------------------------------------------------------------------------------------------------------------
%%% preload the polyhedron geometric model from TXT file with data format ICQ
% [Vol,Den,N_p,N_f,N_e,p_data,f_data,e_data,len_data,ne_data,nf_data] = loadpoly('Apophis99942.txt');
fid = fopen('Rock1.bt');
tmp = fscanf(fid,'%d %d',[2,1]);
N_p = tmp(1);
N_f = tmp(2);
p_data = fscanf(fid,'%f %f %f',[3,N_p]);
p_data = [(1:N_p).',p_data'];
tmp = fscanf(fid,'%d %d %d',[3,N_f]);
tmp = tmp.';
f_data = [(1:N_f).',tmp(:,1:3)];
DT = delaunayTriangulation(p_data(:,2:4));
[~, Vol] = convexHull(DT);
disp('Volume: ');
disp(Vol);
%----------------------------------------------------------------------------------------------------------------------------------------------------
%%% Normalization
L_unit = (3*Vol/4/pi)^(1/3); % unit length in m
T_unit = Period; % unit time in s
p_data(:,2:4) = p_data(:,2:4)/L_unit;
% len_data(:,2) = len_data(:,2)/L_unit;
%
%----------------------------------------------------------------------------------------------------------------------------------------------------
%%% calculate the Spherical Harmonics and save as SH.bt
if flag_genSH
    MSH2COEF(p_data(:,2:4),5);
end
%----------------------------------------------------------------------------------------------------------------------------------------------------
%%% load the Spherical Harmonics 
if flag_loadSH
    CSALL = load('SH.bt');
    lmax=CSALL(end,1);
end
%----------------------------------------------------------------------------------------------------------------------------------------------------
%%% draw the contour surface
if flag_SHshape
    %..................SHMshape(CSALL,lmax,degree,Nmesh) 
    [XMSH,YMSH,ZMSH] = SHMshape(CSALL,lmax,5,500);
    %
    figure;
    p=patch('Vertices',p_data(:,2:4),'Faces',f_data(:,2:4)+1,'FaceVertexCData',[0.7,0.4,0.4],'FaceColor','flat');
    set(p,'edgecolor','none');
    hold on
    
    p1=surf(XMSH,YMSH,ZMSH,'FaceColor','flat');    
    grid on
    material dull
    %camlight;
    %lighting gouraud
    view(3);
    axis tight
    axis equal
%     xlim([-1.5,1.5]);
%     ylim([-1.5,1.5]);
%     zlim([-1.5,1.5]);
end
%
%----------------------------------------------------------------------------------------------------------------------------------------------------
if flag_triSH
    % generate mesh
    addpath('./spheretri/');  
    [Vertices, Faces] = spheretri(500);
    tmp = Faces(:,1);
    Faces(:,1) = Faces(:,2);
    Faces(:,2) = tmp;
    Vertices1 = Vertices*0;
    [Azimuth,Elevation,tmp] = cart2sph(Vertices(:,1),Vertices(:,2),Vertices(:,3));
    Theta = Elevation+pi/2;
    Phi = Azimuth+pi;
    for i=1:length(Azimuth)
        rvec = SHfunvec(Theta(i),Phi(i),CSALL,lmax,8); 
        Vertices1(i,:) = rvec.'; 
    end
    %
    figure;
    p=patch('Vertices',Vertices1,'Faces',Faces,'FaceVertexCData',[0.7,0.4,0.4],'FaceColor','flat');
    hold on
    T = Faces;
    Xb = Vertices1;
    TR = triangulation(T,Xb);
    trisurf(T,Xb(:,1),Xb(:,2),Xb(:,3), ...
        'FaceColor', 'cyan', 'faceAlpha', 0.8);
    axis equal;
    hold on;
    P = incenter(TR);
    fn = faceNormal(TR);
    quiver3(P(:,1),P(:,2),P(:,3), ...
        fn(:,1),fn(:,2),fn(:,3),0.5, 'color','r');
    %
    material dull
    %camlight;
    %lighting gouraud
    view(3);
    axis tight
    axis equal
    DT = delaunayTriangulation(Vertices1);
    [~, VolTri] = convexHull(DT);
    L_unit_Tri = (3*VolTri/4/pi)^(1/3);
    Vertices1 = Vertices1 * Radius / L_unit_Tri;
    DT = delaunayTriangulation(Vertices1);
    [~, VolTri] = convexHull(DT);
    disp('New Volume:');
    disp(VolTri);
    %
%     fid = fopen('Bennu.vtk','w+');
%     fprintf(fid,'# vtk DataFile Version 3.0\n');
%     fprintf(fid,'Asteroid\n');
%     fprintf(fid,'ASCII\n');
%     fprintf(fid,'DATASET POLYDATA\n');
%     fprintf(fid,'POLYGONS %d %d\n',length(Faces),length(Faces)*4);
%     for ii = 1:length(Faces)
%         fprintf(fid,'%d %d %d %d\n',3,Faces(ii,:)-1);
%     end
%     fprintf(fid,'POINTS %d float\n',length(Vertices1));
%     for ii = 1:length(Vertices1)
%         fprintf(fid,'%f %f %f\n',Vertices1(ii,:));
%     end    
%     
%     writeBT(Vertices1,Faces,'Bennu.bt');
end











