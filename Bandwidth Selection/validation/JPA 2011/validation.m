format long;
for i=1:100
name = strcat('validation_ce_', int2str(i),'.mat');   
load(name)
nameBand=strcat('CV_ce_M_', int2str(i),'.csv');
nameBandExp=strcat('CV_ce_exp_M_', int2str(i),'.csv');
bands=zeros(100,729);
bandsExp=zeros(100,729);
for k=1:729
for j=1:100
bands(j,k)     =  Ker_LSCV_OUT(h(:,k), X(j,:), X, Y, Z, n, q, d);
bandsExp(j,k)     =  Ker_LSCV_OUT(hexp(:,k), X(j,:), X, Yexp, Z, n, q, d);
end
end
disp(i)
dlmwrite(nameBand, bands,'delimiter','\t','precision', 22)
dlmwrite(nameBandExp, bandsExp,'delimiter','\t','precision', 22)
end