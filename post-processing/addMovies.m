%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           combine Movie
%           input: several movies
%           output: Movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% parameters
moviePath = 'C:\Users\chengbin\Desktop\';

%%
filename = strcat(moviePath,'out1.avi');
aviobj1 = VideoReader(filename);

filename = strcat(moviePath,'out2.avi');
aviobj2 = VideoReader(filename);

%%
fps = aviobj1.FrameRate;
% N = aviobj1.NumberOfFrames + aviobj2.NumberOfFrames;

aviobj = VideoWriter('out.avi'); %��ʼ��һ��avi�ļ�
fps = aviobj.FrameRate;
open(aviobj);

while hasFrame(aviobj1)
    frames = readFrame(aviobj1); %avi��Ҫ��ɫ����
    writeVideo(aviobj,frames); %һ֡һ֡��д��avi
end

while hasFrame(aviobj2)
    frames = readFrame(aviobj2); %avi��Ҫ��ɫ����
    writeVideo(aviobj,frames); %һ֡һ֡��д��avi
end
close(aviobj); %����������д��avi
