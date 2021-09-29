%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           generate TriMesh of points
%           input:  points data
%           output: bt file
%
%           description: TriMesh calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
format long;
%% Parameters
inFile = 'points.txt';
outFile = 'trimesh.bt';
%% Load data
load('points.txt');

%% Triangulation
DT = delaunayTriangulation(points);
[T,Xb] = freeBoundary(DT);
TR = triangulation(T,Xb);
trisurf(T,Xb(:,1),Xb(:,2),Xb(:,3), ...
     'FaceColor', 'cyan', 'faceAlpha', 0.8);
axis equal;
hold on;

