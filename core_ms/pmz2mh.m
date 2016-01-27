function pmass = pmz2mh(pmz,pcharge)
% PMZ2MH calculate pmass (MH+) from pmz and pcharge.
%	See also MH2PMZ.

pmz = pmz(:);
pcharge = pcharge(:);

pmass = pmz.*pcharge-(pcharge-1)*1.0072764668;

%%  benchmarking various pmass calcs
%%% Results:
%   method   t1     t2
%   1        0.85   1.93 seconds
%   2        1.01   3.25 seconds
%   3        1.04   3.5  seconds
%   4        1.55   5.4  seconds
%   5        2.2    5.8  seconds
%   6        1.25   5.6  seconds
% Timings are each avg of 3 reps with params:
%	t1: n=100000; nreps=10000;
%	t2: n=10000; nreps=100000;
% Method 1 is fastest, but it uses its own definition of proton mass. The others
%	are better in that they refer to CONSTS.mProton

%%% timing code:
%{
pmz = rand(n,1)*2000;
pcharge = randi(4,n,1);
for i = 1:10, t = pmz2mh(pmz,pcharge); end
tic
for i = 1:nreps, t = pmz2mh(pmz,pcharge); end
toc
%}

%%% Methods
%{
%%% 1: 0.85, 1.93 seconds
pmass = pmz.*pcharge-(pcharge-1)*1.0072764668;

%%% 2: 1.01, 3.25 seconds
% define global 'mProton' outside of this function:
% global mProton; mProton = CONSTS.mProton;
pmass = pmz.*pcharge-(pcharge-1)*mProton;

%%% 3: 1.04, 4.0 seconds
% define global 'consts' outside of this function:
%	global consts; consts = CONSTS;
global consts
pmass = pmz.*pcharge-(pcharge-1)*consts.mProton;

%%% 4: 1.55, 5.4 seconds
pmass = pmz.*pcharge-(pcharge-1)*CONSTS.mProton;

%%% 5: 2.2, 5.8 seconds
persistent mProton
if isempty(mProton), mProton = CONSTS.mProton; end
pmass = pmz.*pcharge-(pcharge-1)*mProton;

%%% 6: 2.2, 7.6  seconds
% setappdata(0,'mProton',CONSTS.mProton)
mProt = getappdata(0,'mProton');
pmass = pmz.*pcharge-(pcharge-1)*mProt;

%}