%
% gen the mesh shape of Spherical Harmonics
%
function rvec = SHfunvec(theta,phi,CSALL,lmax,degree)
%
% theta = 0...pi, scalar
% phi = 0...2*pi, scalar
%
R = SHfun(theta,phi,CSALL,lmax,degree);
%
nvec = [sin(theta)*cos(phi);
        sin(theta)*sin(phi);
        cos(theta)];
%
rvec = nvec * R;
%
end
%
%
%
