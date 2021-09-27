%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% Mesh grid points coordinates in unit 1
% 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
function MSH2COEF(p_data,lmax)
%
N=length(p_data(:,1));
L=(lmax+1)^2;
if L>N
    disp('Warning: Max degree beyond mesh info covers!');
    return;
end
TPR=zeros(N,3); %Mesh point coordinates Theta-Phi-R
for i=1:N
    pt=p_data(i,:);
    TPR(i,1)=pi/2-atan2(pt(3),sqrt(pt(1)^2+pt(2)^2));
    TPR(i,2)=atan2(pt(2),pt(1));
    if TPR(i,2)<0
        TPR(i,2)=TPR(i,2)+2*pi;
    end
    TPR(i,3)=norm(pt,2);
end
NLM=zeros(1,L);
for l=0:lmax
    NLM(l^2+1)=sqrt((2*l+1)/4/pi);
    if l>0
        for m=1:l
            NLM(l^2+1+m)=sqrt((2*l+1)/4/pi*factorial(l-m)/factorial(l+m));
            NLM(l^2+l+1+m)=NLM(l^2+1+m);
        end
    end
end
PLM=zeros(N,L);
for l=0:1:lmax
    plm=legendre(l,cos(TPR(:,1)));
    PLM(:,l^2+1)=plm(1,:).';
    if l>0
        PLM(:,l^2+2:l^2+1+l)=plm(2:l+1,:).';
        PLM(:,l^2+2+l:l^2+1+2*l)=PLM(:,l^2+2:l^2+1+l);
    end
end
QCS=ones(N,L);
sqr2=sqrt(2);
for l=1:lmax
    tmp=TPR(:,2)*(1:l);
    QCS(:,l^2+2:l^2+1+l)=sqr2*cos(tmp);
    QCS(:,l^2+2+l:l^2+1+2*l)=sqr2*sin(tmp);
end
MF=QCS.*(ones(N,1)*NLM).*PLM;
clear NLM PLM QCS tmp plm MeshPoint
LFT=MF.'*MF;
RGT=MF.'*TPR(:,3);
RES=LFT\RGT;
residual=norm(LFT*RES-RGT);
%
disp('The condition number of Left Matrix: ');
disp(cond(LFT));
disp('The residual: ');
disp(residual);
disp('The mean magnitude of SHcoef: ');
disp(norm(RES)/L);
%
%clear LFT RGT MF TPR
%write to reconding file
REC=zeros((lmax+1)*(lmax+2)/2,4);
for l=0:lmax
    REC((l^2+l)/2+1,1)=l;
    REC((l^2+l)/2+1,3)=RES(l^2+1);
    if l>0
        for m=1:l
            REC((l^2+l)/2+1+m,1)=l;
            REC((l^2+l)/2+1+m,2)=m;
            REC((l^2+l)/2+1+m,3)=RES(l^2+1+m);
            REC((l^2+l)/2+1+m,4)=RES(l^2+l+1+m);
        end
    end
end
fid=fopen('SH.bt','wt');
fmt='%8d %8d %28.16e %28.16e\n';
for i=1:length(REC(:,1))
    fprintf(fid,fmt,REC(i,:));
end
fclose(fid);
end
