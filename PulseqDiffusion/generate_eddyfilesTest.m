%% Test 3 directions, Nb0s=1, Nb=1, Ndir=3, 10 slices
dimpars=struct('Nb0s',1,'Nb',1,'Ndir',3,'Nsl',10)
generate_eddyfiles(dimpars)
indexpred = ones(4,1);
acqpred = [0 1 0 0.05];
assertFiles(indexpred,acqpred);

%% Test 3 directions, Nb0s=1, Nb=5, Ndir=3, 10 slices
dimpars=struct('Nb0s',1,'Nb',5,'Ndir',3,'Nsl',10)
generate_eddyfiles(dimpars)
indexpred = ones(16,1);
acqpred = [0 1 0 0.05];
assertFiles(indexpred,acqpred);


%% Test 6 directions, Nb0s=1, Nb=5, Ndir=6, 20 slices
dimpars=struct('Nb0s',1,'Nb',5,'Ndir',6,'Nsl',20)
generate_eddyfiles(dimpars)
indexpred = ones(31,1);
acqpred = [0 1 0 0.05];
assertFiles(indexpred,acqpred);

%% Test 12 directions, Nb0s=1, Nb=1, Ndir=12, 4 slices
dimpars=struct('Nb0s',1,'Nb',1,'Ndir',12,'Nsl',4)
generate_eddyfiles(dimpars)
indexpred = ones(13,1);
acqpred = [0 1 0 0.05];
assertFiles(indexpred,acqpred);


%% Test 60 directions, Nb0s=3, Nb=1, Ndir=60, 20 slices
dimpars=struct('Nb0s',3,'Nb',1,'Ndir',60,'Nsl',20)
generate_eddyfiles(dimpars)
indexpred = ones(63,1);
acqpred = [0 1 0 0.05];
assertFiles(indexpred,acqpred);
