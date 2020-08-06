function assertFiles(indexpred,acqpred)
if not(exist('index','file'))
    assert(0,'File index not found');
else
    index=load('index');
end

if not(exist('acqp','file'))
    assert(0,'File acqp not found');
else
    acqp=load('acqp');
end
 
tol = 1e-6;
tf = norm(index-indexpred) <= tol;
assert(tf, 'index values do not match prediction');

tf = norm(acqp-acqpred) <= tol;
assert(tf, 'acqp values do not match prediction');

delete('index','acqp')
end      
