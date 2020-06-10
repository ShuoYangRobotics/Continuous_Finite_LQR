function [A, B] = getLinTraj2(t,tRef, xRef,yRef, uRef)

%Linearizes the system about the nominal trajectory. This function is used
%to turn a non-linear trajectory tracking problem into a time-varying
%linear problem.
if (abs(t - tRef(end)) < 1e-3)
    x = xRef(end);
    y = yRef(end);
    u = uRef(end);    
else
    x = interp1(tRef, xRef, t);
    y = interp1(tRef, yRef, t);
    u = interp1(tRef, uRef, t);
end
[A,B] = getLinSys(x,y,u);

end