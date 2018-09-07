function dyds = EPA(~,y)

global t E Ep fy P Qc yc I

A = Ep*t^3/12;
B = 1/3*(E-Ep)*(fy/E)^3;

if 0 <= y(1) && y(1) < pi/2
    C = -P*sin(y(1))-Qc*sin(y(1))/sin(yc);
elseif pi/2 <= y(1) <= yc
    C = -Qc*sin(y(1))/sin(yc);
end

if y(1) == 0
    Tc = 0;
else
    Tc = -Qc/tan(y(1));
end

dyds = zeros(2,1);    % a column vector

if 0 <= y(1) && y(1) < pi/2
    dyds(1) = sqrt(2/E/I)*sqrt(P*cos(y(1))+Tc*(1-cos(y(1)-yc))-Qc*sin(y(1)-yc));
elseif pi/2 <= y(1) <= yc
    dyds(1) = 0;
end

dyds(2) = C/(A+B/dyds(1)^3);

