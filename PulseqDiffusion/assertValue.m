function assertValue(actVal,expVal)
% Helper function to assert equality when testing the code
% Takes two images and compares them 

tf = isequal(actVal,expVal);
assert(tf, 'Recon images do not match!');
end
