function R = SHfun(theta,phi,CSALL,lmax,degree)
%
% theta = 0...pi, scalar
% phi = 0...2*pi, scalar
%
if degree>lmax
    disp('degree overflow!');
    return
end
len = (degree+1)*(degree+2)/2;
CS = CSALL(1:len(1),3:4);
%
PS = zeros(len,1);
for k = 0:1:degree
    p2 = legendre(k,cos(theta));
    PS(k*(k+1)/2+1:(k+1)*(k+2)/2,:) = p2;
end
%
ks = 1:1:degree;
coss = cos(ks*phi);
sins = sin(ks*phi);
%
R = 0;
id = 0;
for k = 0:1:degree
    for m = 0:1:k
        id = id+1;
        KS = sqrt((2*k+1)/4/pi*factorial(k-m)/factorial(k+m));
        if m==0
            tmp=PS(id)*KS*CS(id,1);
        else            
            tmp=sqrt(2)*PS(id)*(coss(m)*CS(id,1)+sins(m)*CS(id,2))*KS;
        end
        R=R+tmp;     
    end
end
%
end