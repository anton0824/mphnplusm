% get the steady-state probabilities based on pi_i=pi_i*R{i}.

function pii = piscalar(R,U0,n, trunc)

  pii = zeros(1, n+trunc);
  dim = size(U0);
  M = -U0;
  M(:,dim(1)) = ones(dim(1),1);
  vec = zeros(1,dim(1));
  vec(1,dim(1)) = 1;
  pi = vec/M;
  pii(1) = sum(pi);

  for i = 2:n
    pi = pi*R{i};
    pii(i) = sum(pi);
  end
  for i = n+1:n+trunc
      pi = pi*R{n};
      pii(i) = sum(pi);
  end
  total = sum(pii);
  Tmp = inv(eye(size(R{n}))-R{n});
  total = total + sum(pi*Tmp);
  pii = pii/total;