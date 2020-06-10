function Soln = trackingLqr(tRef,xRef,yRef,uRef,Q,R,F,xN,tol)
% Soln = trajectoryLqr(t,linSys,Q,R,F,tol)
%
% This function is used to solve the finite-horizon continuout-time linear
% quadratic regulator problem.
%
% INPUTS:
%   t = monotonically increasing vector of times at which you would like
%       the gain matrix to be calculated
%   xRef,
%   yRef,
%   uRef,
%   linSys = function handle for time-varying linear dynamics
%       [A, B] = linSys(t)
%   Q = state cost matrix 
%   R = input cost matrix 
%   F = final state cost matrix 
%   tol = accuracy of riccati equation propagation
%
% OUTPUTS:
%   Soln = struct array with solution at each point in t
%   Soln(i).t = t(i);
%   Soln(i).K = gain matrix at t(i)
%   Soln(i).S = Cost to go
%   Soln(i).E = close-loop eigen-values for system at t(i)
%
% NOTES:
%
%   J = (x-xN)'F(x-xN) + Integral {x'Qx + u'Ru} dt
%
% See Also LQR

nState = size(Q,1);
nInput = size(R,1);

userFun = @(t,z)rhs(t,z,tRef, xRef, yRef, uRef,Q,R,nState);
Sxx_N = reshape(F,nState*nState,1);
sxx_N = -F*(xN-[xRef(end);yRef(end)]);
s0_N  = (xN-[xRef(end);yRef(end)])'*F*(xN-[xRef(end);yRef(end)]);
costate = [Sxx_N;sxx_N;s0_N];
tSpan = [tRef(end),tRef(1)];

options = odeset();
options.RelTol = 1e-1;
options.AbsTol = 1e-1;
sol = ode45(userFun,tSpan,costate);
costate = deval(sol,tRef);

nSoln = length(tRef);
Soln(nSoln).t = 0;
Soln(nSoln).K = zeros(nState,nInput);
Soln(nSoln).uff = zeros(nInput,1);
Soln(nSoln).Sxx = zeros(nState,nState);
Soln(nSoln).sx = zeros(nState,1);
Soln(nSoln).s0 = 0;

for idx=1:nSoln
    i = nSoln-idx+1;
    costateNow = costate(:,i);
    tNow = tRef(i);
    Sxx = reshape(costateNow(1:nState*nState),nState,nState);
    sx = costateNow(nState*nState+1:nState*nState+nState);
    s0 = costateNow(end);
    [A,B] = getLinTraj2(tNow, tRef, xRef,yRef,uRef);
    Soln(i).t = tNow;
    Soln(i).K = -R\B'*Sxx;
    Soln(i).uff = -R\B'*sx;
    Soln(i).Sxx = Sxx;
    Soln(i).sx = sx;
    Soln(i).s0 = s0;
end
% the last control is 0
% Soln(nSoln).K = [0 0];
% Soln(nSoln).uff = 0;
end

function dz = rhs(t,z,tRef, xRef, yRef, uRef,Q,R,nState)

Sxx = reshape(z(1:nState*nState),nState,nState);
sx = z(nState*nState+1:nState*nState+nState);
s0 = z(end);

[A,B] = getLinTraj2(t, tRef, xRef,yRef,uRef);  % this is get on xRef
dSxx = -(Q-Sxx*B*(R\B')*Sxx+Sxx*A+A'*Sxx); 
dsx  = -((A'-Sxx*B*(R\B'))*sx);
ds0  = -(-sx'*B*(R\B')*sx);
dP = [reshape(dSxx,nState*nState,1);dsx;ds0];
dz = reshape(dP,nState*nState+nState+1,1);
end

% function dz = rhs(t,z,tRef, xRef, yRef, uRef,Q,R,nState)
% 
% Sxx = reshape(z(1:nState*nState),nState,nState);
% sx = z(nState*nState+1:nState*nState+nState);
% s0 = z(end);
% d1 = interp1(tRef, xRef, t);
% d2 = interp1(tRef, yRef, t);
% xd = [d1;d2];
% ud = interp1(tRef, uRef, t);
% [A,B] = getLinTraj2(t, tRef, xRef,yRef,uRef);  % this is get on xRef
% dSxx = -(Q-Sxx*B*(R\B')*Sxx+Sxx*A+A'*Sxx); 
% dsx  = -(-Q*xd+(A'-Sxx*B*(R\B'))*sx+Sxx*B*ud);
% ds0  = -(xd'*Q*xd-sx'*B*(R\B')*sx+2*sx'*B*ud);
% dP = [reshape(dSxx,nState*nState,1);dsx;ds0];
% dz = reshape(dP,nState*nState+nState+1,1);
% end


