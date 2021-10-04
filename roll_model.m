function Dy = roll_model(t,y)

global E nu r e m beta
dn = r*0.1;
%% Model
Kn = -2*E*sqrt(r*dn)/(3*(1-nu^2));
lnCOR = log(e);
Cn = 2*sqrt(5/6)*lnCOR/sqrt(lnCOR^2+pi^2)*sqrt(m*E/(1-nu^2))*r^(1/4)*dn^(1/4);
Kr = 0.25*Kn*(beta*r)^2;
Cr = 0.25*Cn*(beta*r)^2;

Dy = zeros(2,1);
Dy(1) = y(2);
Dy(2) = (Kr*y(1)+Cr*y(2))/(0.4*m*r^2);
end