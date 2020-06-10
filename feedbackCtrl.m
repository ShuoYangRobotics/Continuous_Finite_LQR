function u = feedbackCtrl(t,z,tSol,xSol,ySol,uff,K1,K2)

cur_uff = interp1(tSol, uff, t);
cur_K = [interp1(tSol, K1, t);interp1(tSol, K2, t)];
cur_ref = [interp1(tSol, xSol, t);interp1(tSol, ySol, t)];

u = cur_uff + cur_K(1,:).*(z(1,:)-cur_ref(1,:)) ...
    + cur_K(2,:).*(z(2,:)-cur_ref(2,:));

end