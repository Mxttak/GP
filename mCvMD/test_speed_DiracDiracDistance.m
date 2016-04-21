density1 = rand(100,2000);
density2 = rand(100,2000);

mean1 = mean(density1,2);
mean2 = mean(density2,2);

bmax = 1e5;
dist = ClassDiracDiracDistance(bmax);

tic
dist1 = dist.ComputeDistance(density1,density2);
t1 = toc;

tic
dist2 = DiracDiracDistance(density1,mean1,density2,mean2,bmax);
t2 = toc;