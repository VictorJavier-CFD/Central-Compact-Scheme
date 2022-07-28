%**************************************************************************
%  Function to calculate a derivative by CCS    
%    y   -> df/dx vector 
%    var -> f vector 
%    nn  -> number of nodes
%    hh  -> interval length
%    th  -> order of CCS
%**************************************************************************

function y = CCS(var,nn,hh,th)
  %% Matrix
  Ms = sparse(nn,nn); Qs = sparse(nn,nn);
  %% Builder
  switch th
    case 4
      % Left boundary
      i = 1; 
      C = i; E = C + 1; EE = C + 2; EEE = C + 3;
      Ms(C,C) = 1; Ms(C,E) = 3;
      Qs(C,C) = -17/6; Qs(C,E) = 3/2; Qs(C,EE) = 3/2; Qs(C,EEE) = -1/6; 
      % Inner
      for i = 2:nn-1
        C = i; W = C - 1; E = C + 1;
        Ms(C,W) = 1/4; Ms(C,C) = 1; Ms(C,E) = 1/4;
        Qs(C,W) = -3/4; Qs(C,C) = 0; Qs(C,E) = 3/4;
      end
      % Right boundary
      i = nn; 
      C = i; W = C - 1; WW = C - 2; WWW = C - 3;
      Ms(C,W) = 3; Ms(C,C) = 1;
      Qs(C,WWW) = 1/6; Qs(C,WW) = -3/2; Qs(C,W) = -3/2; Qs(C,C) = 17/6;       
    case 6
      % Left boundary
      i = 1;
      C = i; E = C + 1; EE = C + 2; EEE = C + 3; EEEE = C + 4; EEEEE = C + 5;
      Ms(C,C) = 1; Ms(C,E) = 5;
      Qs(C,C) = -197/60; Qs(C,E) = -5/12; Qs(C,EE) = 5; Qs(C,EEE) = -5/3; Qs(C,EEEE) = 5/12; Qs(C,EEEEE) = -1/20; 
      % Near left boundary
      i = 2;
      C = i; W = C - 1; E = C + 1; EE = C + 2; EEE = C + 3;
      Ms(C,W) = 1/8; Ms(C,C) = 1; Ms(C,E) = 3/4;
      Qs(C,W) = -43/96; Qs(C,C) = -5/6; Qs(C,E) = 9/8; Qs(C,EE) = 1/6; Qs(C,EEE) = -1/96;
      % Inner
      for i = 3:nn-2
        C = i; W = C - 1; WW = C - 2; E = C + 1; EE = C + 2;
        Ms(C,W) = 1/3; Ms(C,C) = 1; Ms(C,E) = 1/3;
        Qs(C,WW) = -1/36; Qs(C,W) = -7/9; Qs(C,C) = 0; Qs(C,E) = 7/9; Qs(C,EE) = 1/36; 
      end
      % Near right boundary
      i = nn-1;
      C = i; W = C - 1; WW = C - 2; WWW = C - 3; E = C + 1;
      Ms(C,W) = 3/4; Ms(C,C) = 1; Ms(C,E) = 1/8;
      Qs(C,WWW) = 1/96; Qs(C,WW) = -1/6; Qs(C,W) = -9/8; Qs(C,C) = 5/6; Qs(C,E) = 43/96; 
      % Right boundary
      i = nn;
      C = i; W = C - 1; WW = C - 2; WWW = C - 3; WWWW = C - 4; WWWWW = C - 5;
      Ms(C,W) = 5; Ms(C,C) = 1;
      Qs(C,WWWWW) = 1/20; Qs(C,WWWW) = -5/12; Qs(C,WWW) = 5/3; Qs(C,WW) = -5; Qs(C,W) = 5/12; Qs(C,C) = 197/60;
  end
  %% Solver
  M = full(Ms); Q = full(Qs);
  y = linsolve(M,Q*var'/hh)';
end