function [IRF, VARopt] = VARir(VAR,VARopt)
% =======================================================================
% Compute IRFs for a VAR model estimated with VARmodel. Three
% identification schemes can be specified: zero short-run restrictions,
% zero long run restrictions, and sign restrictions
% =======================================================================
% [IRF VARopt] = VARir(VAR,VARopt)
% -----------------------------------------------------------------------
% INPUT
%   - VAR: structure, result of VARmodel function
%   - VARopt: options of the IRFs (see VARoption)
% ----------------------------------------------------------------------- 
% OUTPUT
%   - IRF(t,j,k): matrix with 't' steps, containing the IRF of 'j' variable 
%       to 'k' shock
%   - VARopt: options of the IRFs (see VARoption)
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

% Note. This code follows the notation as in the lecture notes available at
% https://sites.google.com/site/ambropo/MatlabCodes

% Specifically if Y, u, and e are [NxT] matrices, the following is true:
% Reduced form VAR    -->  Y = F*Y(-1) + u    (F => Ft' from VARmodel function) 
% Structural form VAR --> AY = B*Y(-1) + e     

% Where:
%     F = invA*B  -->  B = A*F  -->  B = IRF.invA\transpose(VAR.Ft)
%     u = invA*e  -->  e = A*u  -->  e = IRF.invA\transpose(VAR.residuals);

% Impulse responses:
% IRF(1) = invA*e          where "e" is impulse in the code
% IRF(j) = H(j)*IRF(1)     where H(j)=H^j and for j=2,...,nsteps



%% Check inputs
%===============================================
if ~exist('VARopt','var')
    error('You need to provide options (see VARoption)');
end

nsteps = VARopt.nsteps;
impact = VARopt.impact;
shut   = VARopt.shut;
ident  = VARopt.ident;
S      = VARopt.S;
IV     = VARopt.IV;


%% Retrieve parameters and preallocate variables
%===============================================
det   = VAR.det;
nvar  = VAR.nvar;
nlag  = VAR.nlag;
Ft    = VAR.Ft;
sigma = VAR.sigma;

IRF = nan(nsteps,nvar,nvar);

%% Compute the impulse response
%==============================

for mm=1:nvar

    % Compute the companion matrix
    F = Ft';
    Fcomp = [F(:,1+det:nvar*nlag+det); eye(nvar*(nlag-1)) zeros(nvar*(nlag-1),nvar)];
    Fcomp_eye = eye(size(Fcomp,1));
    if shut~=0
        Fcomp(shut,:) = 0;
    end
      
    % Initialize the impulse response vector
    response = zeros(nvar, nsteps);
    
    % Create the impulse vector
    impulse = zeros(nvar,1); 
    
    % Get the matrix invA containing the structural impulses
    if strcmp(ident,'oir')
        [out, chol_flag] = chol(sigma);
        if chol_flag~=0; error('VCV is not positive definite'); end
        invA = out';
    elseif strcmp(ident,'bq')
        Finf_big = inv(eye(length(Fcomp))-Fcomp); % from the companion
        Finf = Finf_big(1:nvar,1:nvar);
        D  = chol(Finf*sigma*Finf')'; % identification: u2 has no effect on y1 in the long run
        invA = Finf\D;
    elseif strcmp(ident,'sr')
        [out, chol_flag] = chol(sigma);
        if chol_flag~=0; error('VCV is not positive definite'); end
        if isempty(S); error('Rotation matrix is not provided'); end
        invA = (out')*(S');
    else
    end

    if ~strcmp(ident,'iv')
        % Set the size of the shock
        if impact==0
            impulse(mm,1) = 1; % one stdev shock
        elseif impact==1
            impulse(mm,1) = 1/invA(mm,mm); % unitary shock
        else
            error('Impact must be either 0 or 1');
        end

        % First period impulse response (=impulse vector)
        response(:,1) = invA*impulse;

    elseif strcmp(ident,'iv')
        % Recover residuals (order matter!)
        up = VAR.residuals(:,1);     % shock to be identified
        uq = VAR.residuals(:,2:end); % other shocks

        % Make sample of IV comparable with up and uq
        [aux, fo, lo] = CommonSample([up IV(VAR.nlag+1:end)]);
        p = aux(:,1);
        q = uq(end-length(p)+1:end,:);
        Z = aux(:,2:end);

        % Run first stage regression and fitted
        OLS = OLSmodel(p,Z,1);
        p_hat = OLS.yhat;

        % First period impulse response (=impulse vector)
        response(1,1) = 1;
        for ii=2:nvar
            second = OLSmodel(q(:,ii-1),p_hat,0);
            response(ii,1) = second.beta;
        end
    end
    
    % Shut down the response of "shut"
    if shut~=0
        response(shut,1) = 0;
    end

    % Make it comparable with companion
    impulse_big  = [response(:,1)' zeros(1, nvar*(nlag-1))]';
    
    % Recursive computation of impulse response
    for kk = 2:nsteps
        Fcomp_eye = Fcomp_eye * Fcomp; % this is the multiplier Fcomp^n
        response_big   = Fcomp_eye * impulse_big;
        response(:,kk) = response_big(1:nvar);
    end
    IRF(:,:,mm) = response';
end

% Update VARopt
VARopt.Fcomp = Fcomp;
VARopt.maxEig = max(abs(eig(Fcomp)));
if ~strcmp(ident,'iv')
    VARopt.invA      = invA;
else
    VARopt.IVfstage = OLS;
    aux = [nan(VAR.nlag,1); nan(fo,1); p_hat];
    VARopt.IVfitted = aux;
    VARopt.IVF = OLS.F;
    VARopt.IVR2 = OLS.rsqr;
end

