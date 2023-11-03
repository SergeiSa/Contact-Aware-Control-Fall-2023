function suit = svd_suit(M, tol)

suit.self = M;

[U, S, V] = svd(M);
suit.U = U;
suit.S = S;
suit.V = V;

if min(size(S)) > 1
    s = diag(S);
else
    s = S(1);
end
suit.s = s;

%matlab-style automatic tolerance
if nargin < 2 
    suit.tol = max(size(M)) * eps(norm(s, inf));
else
    suit.tol = tol;
end

%rank
suit.rank = sum(s > suit.tol);

%pinv
V(:,(suit.rank+1):end) = [];
U(:,(suit.rank+1):end) = [];
s((suit.rank+1):end) = [];
s = 1./s(:);
suit.pinv = (V.*s.')*U';

%null
V = suit.V;
suit.null = V(:,(suit.rank+1):end);

%row_space
suit.row_space = V(:,1:suit.rank);

%orth (column space)
U = suit.U;
U(:, (suit.rank+1):end) = [];
suit.orth = U;

%left_null
U = suit.U;
U(:, 1:suit.rank) = [];
suit.left_null = U;
end
