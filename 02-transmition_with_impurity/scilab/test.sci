

function H = GeneraH(N,eps,gamma)

  H = zeros(N,N)


  for i =2:N;
    H(i,i-1) = gamma;
    H(i,i)   = eps/2;
  end
  H = H + H' 
  H(1,1) = eps


end
