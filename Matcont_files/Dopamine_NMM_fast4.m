function out = Dopamine_NMM_fast4
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = [];
out{6} = [];
out{7} = [];
out{8} = [];
out{9} = [];

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
Tu=50;
gaa=1;
Ea=0;
Eg=-80;
Eu=-90;
Sja=15;
Ja=15;
b=1;
sgnum=2*(1/(1+exp(-par_K*(kmrgd(2)-par_Vr))))-1;
dydt=[(par_gammaa*kmrgd(1)*(2*kmrgd(2)-par_Vr-par_Vt)-kmrgd(1)*(par_Md+par_B)*gaa*kmrgd(4)+par_Delta*par_gammaa*par_gammaa*(kmrgd(2)-par_Vr)*sgnum/(par_C*pi))/par_C;
(par_gammaa*(kmrgd(2)-par_Vr)*(kmrgd(2)-par_Vt)-kmrgd(3)+par_Iext+(par_Md+par_B)*gaa*kmrgd(4)*(Ea-kmrgd(2))-par_C*pi*kmrgd(1)*(par_C*pi*kmrgd(1)/par_gammaa+par_Delta*sgnum))/par_C;
(b*(kmrgd(2)-par_Vr)-(kmrgd(3)-Eu))/Tu+par_ujump*kmrgd(1);
-kmrgd(4)/par_Ta+Sja*kmrgd(1)+Ja*par_cexc;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(Dopamine_NMM_fast4);
y0=[0,0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',[],'HessiansP',[]);
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
jac=[ -(kmrgd(4)*(par_B + par_Md) + par_gammaa*(par_Vr - 2*kmrgd(2) + par_Vt))/par_C , (2*kmrgd(1)*par_gammaa + (par_Delta*par_gammaa^2*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1))/(par_C*pi) + (2*par_K*par_Delta*par_gammaa^2*exp(-par_K*(kmrgd(2) - par_Vr))*(kmrgd(2) - par_Vr))/(par_C*pi*(exp(-par_K*(kmrgd(2) - par_Vr)) + 1)^2))/par_C , 0 , -(kmrgd(1)*(par_B + par_Md))/par_C ; -(par_C*pi*(par_Delta*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1) + (kmrgd(1)*par_C*pi)/par_gammaa) + (kmrgd(1)*par_C^2*pi^2)/par_gammaa)/par_C , -(kmrgd(4)*(par_B + par_Md) - par_gammaa*(kmrgd(2) - par_Vr) - par_gammaa*(kmrgd(2) - par_Vt) + (2*kmrgd(1)*par_C*par_K*par_Delta*pi*exp(-par_K*(kmrgd(2) - par_Vr)))/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1)^2)/par_C , -1/par_C , -(kmrgd(2)*(par_B + par_Md))/par_C ; par_ujump , 1/50 , -1/50 , 0 ; 15 , 0 , 0 , -1/par_Ta ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
jacp=[ (kmrgd(1)*par_gammaa*(par_Vr - 2*kmrgd(2) + par_Vt) + kmrgd(1)*kmrgd(4)*(par_B + par_Md) - (par_Delta*par_gammaa^2*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1)*(kmrgd(2) - par_Vr))/(par_C*pi))/par_C^2 - (par_Delta*par_gammaa^2*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1)*(kmrgd(2) - par_Vr))/(par_C^3*pi) , 0 , 0 , -(kmrgd(1)*kmrgd(4))/par_C , -(kmrgd(1)*kmrgd(4))/par_C , -(kmrgd(1)*(par_Vr - 2*kmrgd(2) + par_Vt) - (2*par_Delta*par_gammaa*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1)*(kmrgd(2) - par_Vr))/(par_C*pi))/par_C , -(kmrgd(1)*par_gammaa + (par_Delta*par_gammaa^2*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1))/(par_C*pi) + (2*par_K*par_Delta*par_gammaa^2*exp(-par_K*(kmrgd(2) - par_Vr))*(kmrgd(2) - par_Vr))/(par_C*pi*(exp(-par_K*(kmrgd(2) - par_Vr)) + 1)^2))/par_C , -(kmrgd(1)*par_gammaa)/par_C , (par_gammaa^2*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1)*(kmrgd(2) - par_Vr))/(par_C^2*pi) , 0 , 0 , (2*par_Delta*par_gammaa^2*exp(-par_K*(kmrgd(2) - par_Vr))*(kmrgd(2) - par_Vr)^2)/(par_C^2*pi*(exp(-par_K*(kmrgd(2) - par_Vr)) + 1)^2) ; (kmrgd(3) - par_Iext + kmrgd(2)*kmrgd(4)*(par_B + par_Md) - par_gammaa*(kmrgd(2) - par_Vr)*(kmrgd(2) - par_Vt) + kmrgd(1)*par_C*pi*(par_Delta*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1) + (kmrgd(1)*par_C*pi)/par_gammaa))/par_C^2 - (kmrgd(1)*pi*(par_Delta*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1) + (kmrgd(1)*par_C*pi)/par_gammaa) + (kmrgd(1)^2*par_C*pi^2)/par_gammaa)/par_C , 1/par_C , 0 , -(kmrgd(2)*kmrgd(4))/par_C , -(kmrgd(2)*kmrgd(4))/par_C , ((kmrgd(2) - par_Vr)*(kmrgd(2) - par_Vt) + (kmrgd(1)^2*par_C^2*pi^2)/par_gammaa^2)/par_C , -(par_gammaa*(kmrgd(2) - par_Vt) - (2*kmrgd(1)*par_C*par_K*par_Delta*pi*exp(-par_K*(kmrgd(2) - par_Vr)))/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1)^2)/par_C , -(par_gammaa*(kmrgd(2) - par_Vr))/par_C , -kmrgd(1)*pi*(2/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1) - 1) , 0 , 0 , -(2*kmrgd(1)*par_Delta*pi*exp(-par_K*(kmrgd(2) - par_Vr))*(kmrgd(2) - par_Vr))/(exp(-par_K*(kmrgd(2) - par_Vr)) + 1)^2 ; 0 , 0 , kmrgd(1) , 0 , 0 , 0 , -1/50 , 0 , 0 , 0 , 0 , 0 ; 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 15 , kmrgd(4)/par_Ta^2 , 0 ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,par_C,par_Iext,par_ujump,par_Md,par_B,par_gammaa,par_Vr,par_Vt,par_Delta,par_cexc,par_Ta,par_K)
