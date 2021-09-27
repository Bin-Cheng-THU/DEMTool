%
% gen the mesh shape of Spherical Harmonics
%
%
function [XMSH,YMSH,ZMSH] = SHMshape(CSALL,lmax,degree,Nmesh)
%
N=Nmesh;
%
if degree>lmax
    disp('degree overflow!');
    return
end
len=(degree+1)*(degree+2)/2;
CS=CSALL(1:len(1),3:4);
theta=linspace(0,pi,N);
phi=linspace(0,2*pi,N);
PS=zeros(len,N);
for k=0:1:degree
    p2=legendre(k,cos(theta));
    PS(k*(k+1)/2+1:(k+1)*(k+2)/2,:)=p2;
end
coss=zeros(degree,N);
sins=zeros(degree,N);
for k=1:1:degree
    coss(k,:)=cos(k*phi);
    sins(k,:)=sin(k*phi);
end
KS=zeros(1,length(CS(:,1)));
id=0;
for k=0:1:degree
    for m=0:1:k
        id=id+1;
        KS(id)=sqrt((2*k+1)/4/pi*factorial(k-m)/factorial(k+m));
    end
end
R=zeros(N,N);
id=0;
for k=0:1:degree
    for m=0:1:k
        id=id+1;
        M1=ones(N,1)*PS(id,:);
        if m==0
            tmp=M1*KS(id)*CS(id,1);
        else            
            M2=coss(m,:)'*ones(1,N);
            M3=sins(m,:)'*ones(1,N);
            tmp=sqrt(2)*M1.*(M2*CS(id,1)+M3*CS(id,2))*KS(id);
        end
        R=R+tmp;     
    end
end
%
[TH,PH]=meshgrid(theta,phi);
% [X,Y,Z]=sph2cart(TH,PH,R);
XMSH=zeros(N,N);
YMSH=zeros(N,N);
ZMSH=zeros(N,N);
for i=1:1:N
    for j=1:1:N
        XMSH(i,j)=sin(TH(i,j))*cos(PH(i,j))*R(i,j);
        YMSH(i,j)=sin(TH(i,j))*sin(PH(i,j))*R(i,j);
        ZMSH(i,j)=cos(TH(i,j))*R(i,j);
    end
end
%
end
%