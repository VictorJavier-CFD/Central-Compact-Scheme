%==========================================================================
%
% CCS:  Computes the first-order derivative of a function f = f(xi) using 
% a Central Compact Scheme of m-th order accuracy. A central compact 
% scheme solves a linear system in the form of
%                            df     1
%                          M --- = --- Q f
%                            dxi   Dxi
% where M and Q are banded matrices. Since the Central Compact Scheme
% implemented only uses 4th-, 6th-order accuracy, M is a tridiagonal 
% matrix. 
%
% In case of convergence issues (e.g. matrix M is not diagonal dominant), 
% use a successive over-relaxation method. If the linear system to solve 
% is written as a_{ij} x_j = b_i, then this iterative method starts from
% a guess value, i.e. x_i^{0} = 0, and then apply the formula:   
%      a_{ii}                                   1 - SOR     
%      ------ x_i^k + sum(a_{ij} x_j^k) = b_i + ------- a_{ii} x_i^{k-1}
%       SOR                                       SOR
% where 0 < SOR <= 1 until convergence. The stopping criterion is based on
% reduce the residual so that resid_CCS <= residMIN_CCS. The residual is 
% defined as resid_CCS := ||b - A x^k||_2 / ||b - A x^0||_2 being || ||_2 
% the L2-norm. 
%
%
% Copyright © 2022 Llorente Lázaro, Víctor Javier
% Last Update: July 25, 2022
% Website: https://sites.google.com/view/vjllorente
% Contact: victor.javier.llorente@gmail.com
%
%--------------------------------------------------------------------------
%
% ------------
% Description:
% ------------
%
%  [ mcf(1) mef(1)                                                   ] [   vecder(1)    ]       [ qm(1,1) qm(1,2) qm(1,3) qm(1,4) qm(1,5) qm(1,6)                                                                 ] [   vecfun(1)    ]
%  [ mwf(2) mcf(2) mef(2)                                            ] [   vecder(2)    ]       [ qm(2,1) qm(2,2) qm(2,3) qm(2,4) qm(2,5)                                                                         ] [   vecfun(2)    ]
%  [        mwf(3) mcf(3) mef(3)                                     ] [                ]    1  [ qm(3,1) qm(3,2) qm(3,3) qm(3,4) qm(3,5)                                                                         ] [                ]
%  [               ...    ...    ...                                 ] [      ...       ] = --- [         ...     ...     ...     ...     ...                                                                     ] [      ...       ]
%  [                      ...    ...         ...                     ] [                ]   Dxi [                 ...     ...     ...     ...     ...                                                             ] [                ]
%  [                             mwf(nnod-1) mcf(nnod-1) mef(nnod-1) ] [ vecder(nnod-1) ]       [                         qm(nnod-1,nnod-4) qm(nnod-1,nnod-3) qm(nnod-1,nnod-2) qm(nnod-1,nnod-1) qm(nnod-1,nnod) ] [ vecfun(nnod-1) ]
%  [                                         mwf(nnod)   mcf(nnod)   ] [  vecder(nnod)  ]       [         qm(nnod,nnod-5) qm(nnod,nnod-4)   qm(nnod,nnod-3)   qm(nnod,nnod-2)   qm(nnod,nnod-1)   qm(nnod,nnod)   ] [  vecfun(nnod)  ]
%                                                                                           [   bcf(1)    ]
%                                                                                           [   bcf(2)    ]
%                                                                                           [             ]
%                                                                                         = [     ...     ]
%                                                                                           [             ]
%                                                                                           [ bcf(nnod-1) ]
%                                                                                           [  bcf(nnod)  ]
%
% ----------------
% MATLAB function: 
% ----------------
%   vecder = CCS( vecfun, Dxi, mth )
%   Dependencies: 
%     1. TDMA.m
%
% ------
% INPUT:
% ------
%   Double, (nnod)x1 array      :: vecfun       - function to derive
%   Integer                     :: Dxi          - interval size
%   Integer                     :: mth          - CCS order of accuracy
%
% -------
% OUTPUT:
% -------
%   Double, (nnod)x1 array      :: vecder       - first-order derivative
%
% ----------------
% LOCAL VARIABLES:
% ----------------
%   Double, (nnod)x1 array      :: mcf          - diagonal of M matrix
%   Double, (nnod-1)x1 array    :: mwf, mef     - lower and upper diagonal of M matrix (Note: mwf(1) = mef(nnod) = 0)
%   Double, (nnod)x(nnod) array :: qm           - Q matrix
%   Double, (nnod)x1 array      :: bcf          - discrete source
%   Double                      :: SOR_CCS      - successive over-relaxation parameter
%   Double                      :: residMIN_CCS - minimal residual of SOR
%   Integer                     :: iterMAX      - maximum iteration of SOR
%
% -----------
% REFERENCES:
% -----------
%   1. Lele, S.K. (1992) Compact finite difference schemes with 
%      spectral-like resolution. Journal of Computational Physics, 103(1), 
%      16-42.
%
%==========================================================================
function vecder = CCS( vecfun, Dxi, mth )

    %% determines number of nodes
    nnod = length( vecfun );
    
    %% allocates
    qm = zeros(nnod,nnod);
    mcf = zeros(nnod,1); mef = mcf; mwf = mcf;
    vecder = zeros(nnod,1);
    
    %% build up
    switch mth
        % Fourth order of accuracy
        case 4
            % Matrix M
            % upper diagonal
            mef(1) = 3; mef(2:nnod-1) = 1 / 4;
            % principal diagonal
            mcf = ones(nnod,1);
            % lower diagonal
            mwf(2:nnod-1) = mef(2:nnod-1); mwf(nnod) = mef(1);
            % Matrix Q
            % left boundary
            qm(1,1) = -17 / 6; qm(1,2) = 3 / 2; qm(1,3) = 3 / 2; qm(1,4) = -1 / 6;
            % inner
            for i = 2:nnod-1
                qm(i,i - 1) = -3 / 4; qm(i,i + 1) = -qm(i,i - 1);
            end
            % right boundary
            qm(nnod,nnod - 3) = -qm(1,4); qm(nnod,nnod - 2) = -qm(1,3); qm(nnod,nnod - 1) = -qm(1,2); qm(nnod,nnod) = -qm(1,1);
        % Sixth order of accuracy    
        case 6
            % Matrix M
            % upper diagonal
            mef(1) = 5; mef(2) = 3 / 4; mef(3:nnod-2) = 1 / 3; mef(nnod-1) = 1 / 8;
            % principal diagonal
            mcf = ones(nnod,1);
            % lower diagonal
            mwf(2) = mef(nnod-1); mwf(3:nnod-2) = mef(3:nnod-2); mwf(nnod-1) = mef(2); mwf(nnod) = mef(1);
            % Matrix Q
            % left boundary
            qm(1,1) = - 197 / 6; qm(1,2) = - 5 / 12; qm(1,3) = 5; qm(1,4) = - 5 / 3; qm(1,5) = 5 / 12; qm(1,6) = - 1 / 20;
            % near left boundary
            qm(2,1) = - 43 / 96; qm(2,2) = - 5 / 6; qm(2,3) = 9 / 8; qm(2,4) = 1 / 6; qm(2,5) = - 1 / 96;
            % inner
            for i = 3:nnod - 2
                qm(i,i - 2) = - 1 / 36; qm(i,i - 1) = - 7 / 9; qm(i,i + 1) = -qm(i,i - 1); qm(i,i + 2) = -qm(i,i - 2);
            end
            % near right boundary
            qm(nnod-1,nnod-4) = -qm(2,5); qm(nnod-1,nnod-3) = -qm(2,4); qm(nnod-1,nnod-2) = -qm(2,3); qm(nnod-1,nnod-1) = -qm(2,2); qm(nnod-1,nnod) = -qm(2,1);
            % right boundary
            qm(nnod,nnod - 5) = -qm(1,6); qm(nnod,nnod - 4) = -qm(1,5); qm(nnod,nnod - 3) = -qm(1,4); qm(nnod,nnod - 2) = -qm(1,3); qm(nnod,nnod - 1) = -qm(1,2); qm(nnod,nnod) = -qm(1,1);
    end
    
    %% discrete source
    bcf = ( 1 / Dxi ) * qm * vecfun;
    
    %% successive over-relaxation (SOR)
    SOR_CCS = 1.0;
    iterMAX = 1;
    residMIN_CCS = 1.0e-6;
    
    %% TDMA solver with SOR
    mcf_mod = mcf / SOR_CCS;
    matM = diag(mcf,0) + diag(mwf(2:nnod),-1) + diag(mef(1:nnod-1),1);
    for iter = 1:iterMAX
        bcf_mod = bcf + ( ( 1.0 - SOR_CCS ) / SOR_CCS ) * (mcf .* vecder);
        vecder = TDMA( mwf, mcf_mod, mef, bcf_mod );
        resid_CCS = norm( bcf - matM * vecder, 2 ) / norm( bcf, 2 );
        if resid_CCS <= residMIN_CCS
            break
        end
    end

end