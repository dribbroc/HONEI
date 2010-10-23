function M = spaitau(A,varargin)
% M = spaitau(A,ta,vb)
%
%     Computes the SPAI preconditioner (a left preconditioner)
%     parameters:
%
%       A  (required): input matrix (square, real, and sparse)
%       ta (optional): (tau) only those entries Mij included for which
%                      |Aij| > (1-ta)*max_{j}|Aij|        (for -sc = 1)
%                      i.e. ta=0: main diagonal, ta=1: full pattern of A
%                      Main diagonal always included
%                      default is 0
%       vb (optional): verbose mode
%                      vb = 1 turns it on (default)
%                      vb = 0 turns it off
%
% Note: This constructs a right preconditioner. To construct a left
%       preconditioner use transpose(A).

if (length(varargin) < 1) ta = .4;
else ta = varargin{1}; end

if (length(varargin) < 2) vb = 1;
else vb = varargin{2}; end

%       sc (optional): sparsity control:
%			 0 = adaptive,
%			 1 = specified tau,
%			 2 = fixed diagonals (-ud & -ld)
%			     (main diagonal always included)
%                      default is 0

Mfull = spai_full(A, 0, 0, 5000, 0, 1, vb, 1, 0, 0, ta);
M = spconvert(Mfull);
