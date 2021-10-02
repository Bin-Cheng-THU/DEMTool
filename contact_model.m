function Dy = contact_model(t,y)

global E nu r e m

%% Model
Kn = -2*E*sqrt(r*y(1))/(3*(1-nu^2));
lnCOR = log(e);
Cn = 2*sqrt(5/6)*lnCOR/sqrt(lnCOR^2+pi^2)*sqrt(m*E/(1-nu^2))*r^(1/4)*y(1)^(1/4);
Dy = zeros(2,1);
Dy(1) = y(2);
if y(1)>0
    Dy(2) = (Kn*y(1)+Cn*y(2))/m;
else
    Dy(2) = 0.0;
end
end