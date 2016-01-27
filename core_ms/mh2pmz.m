function pmz = mh2pmz(pmass,pcharge)
% MH2PMZ calculate pmz from pmass (MH+) and pcharge.
%	See also PMZ2MH.

pmz = (pmass+(pcharge-1)*1.0072764668)./pcharge;
% pmz = (pmass+(pcharge-1)*CONSTS.mProton)./pcharge; % MUCH slower