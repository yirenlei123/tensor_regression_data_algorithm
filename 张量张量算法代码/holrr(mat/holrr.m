function W = holrr(X, Y, rank, gamma)
% Noemie Jaquier X = trainx; Y = trainy;gamma = lambda;

d0 = size(X,2);

U = cell(1,length(rank));
Ut = cell(1,length(rank));

%YTmp = tens2mat(Y,1);
YTmp = classical_mode_unfolding(Y,1);
UTmp = (X'*X + gamma.*eye(d0)) \ (X' * (YTmp * YTmp') * X);
[V,D] = eig(UTmp);
[~, ind] = sort(diag(D),'descend');
U{1,1} = V(:,ind(1:rank(1)));

for i = 2:length(rank)
	%YTmp = tens2mat(Y,i);
    YTmp = classical_mode_unfolding(Y,i);
	[V,D] = eig(YTmp*YTmp');
	[~, ind] = sort(diag(D),'descend');
	U{1,i} = V(:,ind(1:rank(i)));
	Ut{1,i} = U{1,i}';
end

Ut{1,1} = (U{1,1}' * (X'*X + gamma.*eye(d0)) * U{1,1}) \ (U{1,1}'*X');

%G = tmprod(Y,Ut,1:length(rank));
G = tensorprod(Y,Ut{1},1,2);
for i = 2:length(rank)
    G = tensorprod(G,Ut{i},1,2);
end

%W = tmprod(G,U,1:length(rank));
W = tensorprod(G,U{1},1,2);
for i = 2:length(rank)
    W = tensorprod(W,U{i},1,2);
end

end
