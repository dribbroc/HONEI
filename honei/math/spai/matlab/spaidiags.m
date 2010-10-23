function M = spaidiags(A,varargin)
% M = spaidiags(A,ld,ud,vb)
%
%     Computes the SPAI preconditioner (a left preconditioner)
%     parameters:
%
%       A  (required): input matrix (square, real, and sparse)
%       ld (optional): no. of lower (below main) diagonals
%                      default is 0
%       ud (optional): no. of upper (above main) diagonals
%                      default is 0
%       vb (optional): verbose mode
%                      vb = 1 turns it on (default)
%                      vb = 0 turns it off
%
% Note: This constructs a right preconditioner. To construct a left
%       preconditioner use transpose(A).


if (length(varargin) < 1) ld = 0;
else ld = varargin{1}; end

if (length(varargin) < 2) ud = 0;
else ud = varargin{2}; end

if (length(varargin) < 3) vb = 1;
else vb = varargin{3}; end

%       sc (optional): sparsity control:
%			 0 = adaptive,
%			 1 = specified tau,
%			 2 = fixed diagonals (-ud & -ld)
%			     (main diagonal always included)
%                      default is 0

Mfull = spai_full(A, 0, 0, 5000, 0, 1, vb, 2, ld, ud, 0);
M = spconvert(Mfull);
