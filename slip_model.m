function Dy = slip_model(t,y)

global E nu r e m

dn = r*0.001;
%% Model
Ks = -2*E/(1+nu)/(2-nu)*sqrt(r*dn);
lnCOR = log(e);
Cs = 2*sqrt(5/6)*lnCOR/sqrt(lnCOR^2+pi^2)*sqrt(2*m*E/(1+nu)/(2-nu))*r^(1/4)*dn^(1/4);

Dy = zeros(2,1);
Dy(1) = y(2);
Dy(2) = (Ks*y(1)+Cs*y(2))/m;

end 