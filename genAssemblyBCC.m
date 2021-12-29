%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate an Assembly system with a BBC configuration
%           output: input files
%
%           description: particle generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% Parameters
Xrange = [-20,20];
Yrange = [-20,20];
Zrange = [0,100];
spacing = 2.0;
Rrange = [1 2];
Srange = 0.1;
Vrange = 0.01;
Rotrange = [-180 180]/180*pi;
inFile = 'pointCell.txt';
outFileCSV = 'configure.csv';
outFileP = 'input_points.txt';
outFileA = 'input_assembly.txt';

%% Particles
x = Xrange(1)+spacing:spacing:Xrange(2)-spacing;
y = Yrange(1)+spacing:spacing:Yrange(2)-spacing;
z = Zrange(1)+spacing:spacing:Zrange(2)-spacing;
pointx = repmat(x',length(y)*length(z),1);
pointy = repmat(y,length(x),length(z));
pointy = reshape(pointy,[length(x)*length(y)*length(z),1]);
pointz = repmat(z,length(x)*length(y),1);
pointz = reshape(pointz,[length(x)*length(y)*length(z),1]);
point = [pointx,pointy,pointz];

radius = rand(size(point,1),1)*(Rrange(2)-Rrange(1)) + Rrange(1);
shift = rand(size(point,1),3)*Srange*2-Srange;
point = point + shift;
velocity = rand(size(point,1),3)*Vrange*2-Vrange;
rotate = rand(size(point,1),3)*(Rotrange(2)-Rotrange(1)) + Rotrange(1);

%% Assembly
fid = fopen(inFile,'r');
numCell = fscanf(fid,'%d',[1,1]);
densityCell = fscanf(fid,'%f',[1,1]);
tmp = fscanf(fid,'%f %f %f %f',[1,4]);
massCell = tmp(1);
inertiaCell(1:3) = tmp(2:4);
pointCell = zeros(numCell,4);
for ii = 1:numCell
    tmp = fscanf(fid,'%f %f %f %f',[1,4]);
    pointCell(ii,:) = tmp;
end
fclose('all');

%% Generator
dataPoint = zeros(size(point,1)*numCell,10);
dataCell = zeros(size(point,1),18);
dataData = zeros(size(point,1)*numCell,6);
for ii=1:size(point,1)
    tmp_mat = angle2dcm(rotate(ii,1),rotate(ii,2),rotate(ii,3),'XYZ');
    tmp_quat = angle2quat(rotate(ii,1),rotate(ii,2),rotate(ii,3),'XYZ');
    for jj=1:numCell
        tmp_point = pointCell(jj,1:3)*radius(ii);
        tmp_point = (tmp_mat'*tmp_point')' + point(ii,:);
        tmp_vel = velocity(ii,:);
        tmp_omega = zeros(1,3);
        tmp_r = pointCell(jj,4)*radius(ii);
        dataPoint((ii-1)*numCell+jj,:) = [tmp_point,tmp_vel,tmp_omega,tmp_r];
        dataData((ii-1)*numCell+jj,:) = [massCell*radius(ii)^3 0.0 pointCell(jj,1:4)*radius(ii)];
    end
    dataCell(ii,:) = [numCell massCell*radius(ii)^3 inertiaCell*radius(ii)^5 point(ii,:) velocity(ii,:) zeros(1,3) tmp_quat(2:4) tmp_quat(1)];
end

%% write csv/point file/assembly file
TagCSV = 1;
TagPoint = 1;
TagAssembly = 1;

if (TagCSV == 1)
    disp('Start writting CSV file...');
    filename = outFileCSV;
    fid = fopen(filename,'w');
    fprintf(fid,'X,Y,Z,U,V,W,W1,W2,W3,R,Cell\r\n');
    for ii = 1:size(dataPoint,1)
        fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d\r\n',dataPoint(ii,:),floor((ii-1)/numCell));
    end
    fclose(fid);
    disp('Writting completed.')
end

if (TagPoint == 1)
    disp('Start writing point file...');
    save(outFileP,'-ascii','dataData');
    disp('Writeing completed.');
end

if (TagAssembly == 1)
    disp('Start writing assembly file...');
    filename = outFileA;
    fid = fopen(filename,'w');
    fprintf(fid,'%d\r\n',size(dataCell,1));
    for ii = 1:size(dataCell,1)
        fprintf(fid,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\r\n',dataCell(ii,:));
    end
    fclose(fid);
    disp('Writeing completed.');
end

fclose('all');