function VARopt = VARoption
% =======================================================================
% Optional inputs for VAR analysis
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

VARopt.nsteps    = 40;
VARopt.impact    = 0;
VARopt.shut      = 0;
VARopt.ident     = 'oir';
VARopt.S         = [];
VARopt.IV        = [];
VARopt.invA      = [];
VARopt.Fcomp     = [];
VARopt.maxEig    = [];
VARopt.ndraws    = 100;
VARopt.pctg      = 90;
VARopt.method    = 'bs';
VARopt.pick      = 0;
VARopt.quality   = 0;
VARopt.suptitle  = 0;
VARopt.firstdate = [];
VARopt.frequency = 'q';
VARopt.IV_fstage = [];
VARopt.IV_fitted = [];
VARopt.IV_F      = [];
VARopt.IV_R2     = [];

