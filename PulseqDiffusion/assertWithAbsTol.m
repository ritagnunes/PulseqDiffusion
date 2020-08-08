function assertWithAbsTol(actVal,expVal,varargin)
% Helper function to assert equality within an absolute tolerance.
% Takes two values and an optional message and compares
% them within an absolute tolerance of 1e-6.
tol = 1e-6;
N = numel(actVal);
dif = reshape( abs(actVal - expVal),[N 1]);
tf = sum(dif) <= tol;
assert(tf, varargin{:});
end  
