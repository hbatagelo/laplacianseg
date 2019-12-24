function [A, b, labels, cindex, xsol] = LClinear_hard(W, maskconstraints)
% ------------------------------------------------------------------------
% Description: Builds the system of linear equations ('hard' version)
% Usage: [A, b, labels] = LClinear_hard(W, maskconstraints);
% 
% Input(s):
%   W: Matrix of weights (sparse/adjacency matrix)
%   maskconstraints:  Matrix of seeds (binary m x n x 2 mask)
%
% Output(s):
%   A: Coefficient matrix
%   b: Matrix containing the right-side of the linear system of equations
%   labels: Seeded labels
%   cindex: complemetary label indices
%   xsol: pre-solution with fixed 'hard' values 
% ------------------------------------------------------------------------

% Initializations
Eps = 1.0;
n = size(W,1);
nRegions = size(maskconstraints,3);
Constrb = zeros(n,nRegions);
%d = [1 1];
labels = [0 1]; %[xB xF]

%%% Dealing with the indices of the seeded regions
indc = 0;
%diag = 0;
for k=1:nRegions
    mask = logical(maskconstraints(:,:,k));
    indb = find(mask);
    Constrb(indb,k) = 1;
    %indc = [indc; indb];
    %diag = [diag; d(k)*ones(length(indb),1)];
end

% Building the graph Laplacian matrix
DW = spdiags(sum(W,2), 0, n, n);
L = Eps*(DW - W);

% Right-side vector b
b = labels(2)*Constrb(:,1) + labels(1)*Constrb(:,2);

% Building the hard constrainted system of linear equations
index = 1:n;
index(~logical(Constrb(:,1) + Constrb(:,2))) = [];
cindex = 1:n;
cindex(index) = [];

% Matrix decompoisiton
Bt = L(cindex, index);
B = transpose(Bt);
xm = b(index);
Lu = L(cindex, cindex);
Lm = L(index, index);

% LC hard system [A|b]
A = (Bt*B + Lu^2);
b = -(Bt*Lm + Lu*Bt)*xm;

% Assigning the constrained pixels to the solution
xsol = zeros(size(xm));
xsol(index) = xm;

end