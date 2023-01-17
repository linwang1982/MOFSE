function SA = KNN_fun(A, kn)
n = size(A,1);
% SA = A - eye(n);
SA = A;
[~,I] = sort(SA,2);
for i = 1:n
    SA(i,I(i,1:(n-kn))) = 0;
end
