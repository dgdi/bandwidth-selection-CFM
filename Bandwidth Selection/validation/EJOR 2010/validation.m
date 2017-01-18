format long;
for i=1:100
name = strcat('validation_ban_', int2str(i),'.mat');   
load(name)
nameBandUni=strcat('CV_ban_Uni_M_', int2str(i),'.csv');
nameBandMul=strcat('CV_ban_Mul_M_', int2str(i),'.csv');
bandsUni=zeros(100,9);
bandsMul=zeros(100,27);
for k=1:81
for j=1:100
bandsUni(j,k)     =  Ker_LSCV_OUT(hUni(:,k), X(j,:), X, YUni, Z(:,1), n, 2, 1);
end
end
for l=1:729
for m=1:100
bandsMul(m,l)     =  Ker_LSCV_OUT(hMul(:,l), X(m,:), X, YMul, Z, n, 2, 2);
end
end
disp(i)
dlmwrite(nameBandUni, bandsUni,'delimiter','\t','precision', 22)
dlmwrite(nameBandMul, bandsMul,'delimiter','\t','precision', 22)
end