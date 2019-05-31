function [amp, rise, decay, hfwidth, rise10_90, tc_dcy] = analyze_EPSC(t,cur)
%
%   assume we are only sent time vs current data of a single EPSC.

%
%  Find the 'easy' stats: amplitude and half-width.
%
[mn, mnI] = min(cur);
[mx, mxI] = max(cur);
Ishift = cur(1)-cur';
amp = abs(cur(mxI)-cur(mnI));
hfmx = (mn+mx)/2;
hfidx = find(cur<hfmx);
hfbds = [hfidx(1) hfidx(end)];
hfwidth = abs(t(hfbds(2))-t(hfbds(1)));

%
%  now fit exponential.  This is hard because the data are negative and
%  near zero, so linearization doesn't work well. 
%

% Only fit to the part of the curve where the EPSC is 'active':
% assume EPSC starts when when abs(dI/dt) > 1e-5 and ends when abs(dI/dt) <
% 1e-5 for the last time. 
%
[dI]=get_dVdt(t,cur);
[EPidx]=find(abs(dI)>1e-5);
EPst  = EPidx(1);
EPend = EPidx(end);

rise = t(mnI) - t(EPst);
[mn10,idx10] = min(find(abs(Ishift)>=0.1*amp));
[mn90,idx90] = min(find(abs(Ishift)>=0.9*amp));
rise10_90 = t(mn90)-t(mn10);

decay = t(EPend)-t(mnI);
EPbds = [EPst mnI EPend];

% now:  decay time constant as fit to exponential.
tidx = mnI:EPend;
yuse = (find(isfinite(log(Ishift(tidx))) == 1 ));
% paramEstsLin = [ones(size(t(tidx))), t(tidx)] \ log(y(tidx));
% paramEstsLin = [ones(size(t(yuse))), t(yuse)] \ log(Ishift(yuse))';
% paramEstsLin(1) = exp(paramEstsLin(1))
% tc_dcy = 1 / paramEstsLin(2);
% % yoffst = 10;
% % pfit = polyfit(t(yuse),log(Ishift(yuse)+yoffst)',1);
% % tc_dcy = 1 / pfit(1);
% % figure; plot(t,Ishift+10,t(tidx),exp(pfit(2))*exp(t(tidx)*pfit(1)));

% http://www.mathworks.com/products/statistics/demos.html?file=/products/de
% mos/shipping/stats/xform2lineardemo.html
modelFun = @(p,x) p(1)*exp(p(2)*x);

% paramEstsLin = [ones(size(x)), x] \ log(y);
% paramEstsLin(1) = exp(paramEstsLin(1))
% 
% paramEsts = nlinfit(x, y, modelFun, paramEstsLin)

paramEstsLin = [ones(size(t(tidx))) t(tidx)] \ log(Ishift(tidx));
paramEstsLin(1) = exp(paramEstsLin(1));
linFit = myExp(paramEstsLin,t(tidx));

paramEsts = nlinfit(t(tidx), Ishift(tidx), modelFun, paramEstsLin);
iFit = modelFun(paramEsts,t(tidx));
doPlot = 0;
if( doPlot)
    close;
    plot(t,Ishift,'.',t(tidx),iFit,'r-',t(tidx),linFit,'k-');
    legend('data','exp fit','linear fit');
end;

tc_dcy = -1 / paramEsts(2);
% fprintf('');


% fprintf('\tamp %g\n\trise %g VS %g 10-90 pcnt\ndecay %g VS %g\thalf width %g\n',amp,rise,rise10_90,decay,tc_dcy,hfwidth);

if( ~isreal(tc_dcy) )
    tc_dcy = real(tc_dcy);
    fprintf('imaginary!\t');
end;
rise = rise10_90;
decay = tc_dcy;

return;


function [yout] = myExp(p,x)
    yout = p(1).*exp(p(2).*x);
    return;
    
    
function [dVdt] = get_dVdt(tvec, vvec)
% function [dVdt,dVall] = get_dVdt(tvec, vvec)
%
%   get_dVdt    use the central difference formula to estimate the
%               derivative of the second vector with respect to the first. 
%
%   INPUT       TVEC, vector of times
%               VVEC, vector of membrane potential, for example.
%   OUTPUT      dVdt, the derivative of VVEC with respect to TVEC.
%
%   Christina Weaver, christina.weaver@mssm.edu, July 2005
%

for i = 2:length(tvec)-1 
    dVdt(i) = (vvec(i+1)-vvec(i-1))/(tvec(i+1)-tvec(i-1));
end;
