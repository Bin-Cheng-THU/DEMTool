%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           convert POV-Ray files to Movie
%           input: POV-Ray files
%           output: Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% parameters
fps = 10;
N = 193;
folderPath = 'C:\Users\chengbin\Desktop\genWallPovRay\genWallPovRay\Data\';
aviobj = VideoWriter('out-yorp-x.avi'); %��ʼ��һ��avi�ļ�
aviobj.FrameRate = fps;
open(aviobj);

for ii=1000:(1000+N-1)
    filename = strcat(folderPath,'Scene');
    filename = strcat(filename,num2str(ii));
    filename = strcat(filename,'.png');
    
    frames=imread(filename); %avi��Ҫ��ɫ����
    writeVideo(aviobj,frames); %һ֡һ֡��д��avi
end
close(aviobj); %����������д��avi