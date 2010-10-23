function M = spai(A,varargin)
% M = spai(A,ep,ns,mn,bs,mb,vb)
%
%     Computes the SPAI preconditioner (a left preconditioner)
%     parameters:
%
%       A  (required): input matrix (square, real, and sparse)
%       ep (optional): epsilon
%                      default is ep = .6
%       ns (optional): number of improvement steps allowed per column
%                      default is ns = 5
%       mn (optional): maximum number of new nonzeros accepted per step
%                      default is mn = 5
%       bs (optional): block size
%                      bs > 0 specifies the constant-block method
%                      bs = 0 specifies the variable-block method
%                      default is bs = 1
%       mb (optional): maximum buffer size
%                      default is 100000
%       vb (optional): verbose mode
%                      vb = 1 turns it on (default)
%                      vb = 0 turns it off
%
% Note: This constructs a right preconditioner. To construct a left
%       preconditioner use transpose(A).


if (length(varargin) < 1) ep = .6;
else ep = varargin{1}; end

if (length(varargin) < 2) ns = 5;
else ns = varargin{2}; end

if (length(varargin) < 3) mn = 5;
else mn = varargin{3}; end

if (length(varargin) < 4) bs = 1;
else bs = varargin{4}; end

if (length(varargin) < 5) mb = 100000;
else mb = varargin{5}; end

if (length(varargin) < 6) vb = 1;
else vb = varargin{6}; end

Mfull = spai_full(A,ep,ns,mb,mn,bs,vb,0,0,0,0);
M = spconvert(Mfull);
