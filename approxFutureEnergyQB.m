function [w] = approxFutureEnergyQB(A,N,B,C,eta,d,verbose)
%  Calculates a polynomial approximation to the future energy function
%  for a quadratic system.
%
%  w = approxFutureEnergyQB(A,N,B,Q,C,eta,d) 
%
%  Computes a degree d polynomial approximation to the future energy function 
%
%          E^+(x) = 1/2 ( w{2}'*kron(x,x) + ... + w{d}'*kron(.. x) )
%
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + N*kron(x,x) + Q*kron(x,u)
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation for W2,
%
%    A'*W2 + W2*A - eta*W2*B*B'*W2 + C'*C = 0,
%  
%  and in the subsequent linear systems arising from the Future H-infinity 
%  Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%  Note that w2 = vec(W2).  Details are in Section III.B of reference [1].
%
%  Requires the following functions from the KroneckerTools repository:
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Reference: [1] Nonlinear balanced truncation: Part 1--Computing energy 
%             functions, by Kramer, Gugercin, and Borggaard, arXiv:2209.07645.
%             [2] Nonlinear balanced truncation for quadratic bilinear
%             systems, by Nick Corbin ... (in progress).
%
%             See Algorithm 1 in [1].
%
%  Part of the NLbalancing repository.
%%

  if (nargin<7)
    verbose = false;
  end

  n = size(A,1);   % A should be n-by-n
  m = size(B,2);   % B should be n-by-m
  p = size(C,1);   % C should be p-by-n

  Q = sparse(n,n*m); % Hardcoded zero Q for now
  Q(1,1) = 1; 
  
  % Create a vec function for readability
  vec = @(X) X(:);


  %% k=2 case
  R = eye(m)/eta;

  if ( eta>0 )
    [W2] = icare(A,B,(C.'*C),R);
    if ( isempty(W2) && verbose )
      warning('approxFutureEnergyQB: icare couldn''t find stabilizing solution')
    end
    
  elseif ( eta<0 )
    [W2] = icare(A,B,(C.'*C),R,'anti');
    
    if ( isempty(W2) && verbose )
      warning('approxFutureEnergyQB: icare couldn''t find stabilizing solution')
      warning('approxFutureEnergyQB: using the hamiltonian')
      [~,W2,~] = hamiltonian(A,B,C.'*C,R,true);
    end
    
  else % eta==0
    [W2] = lyap(A.',(C.'*C));
    
  end

  if (isempty(W2))
    error('approxFutureEnergyQB: Can''t find a stabilizing solution')
  end
  
  %  Check the residual of the Riccati/Lyapunov equation
  if (verbose)
    RES = A'*W2 + W2*A - eta*(W2*B)*(B'*W2) + C'*C ;
    fprintf('The residual of the Riccati equation is %g\n',norm(RES,'inf'));
  end

  %  Reshape the resulting quadratic coefficients
  w{2} = vec(W2);

  %% k=3 case
  if ( d>2 ) % set up the generalized Lyapunov solver
    W2BB = W2*(B*B.');
    Acell = cell(1,d);
    for i=1:d % Pre-computing the left-hand-sides of (46) in [1]
      Acell{i} = A.'-eta*W2BB; 
    end

    b = -LyapProduct(N.',w{2},2);
    [w{3}] = KroneckerSumSolver(Acell(1:3),b,3); % Solve Ax=b for k=3
    
    [w{3}] = kronMonomialSymmetrize(w{3},n,3); % Symmetrize w3
  end
  
  %% k>3 cases (up to d)
  BWk = cell(1,d-1); % Pre-compute B.'*Wi for all the i we need
  QWk = cell(1,d-1); % Pre-compute Q.'*Wi for all the i we need
  BWk{2} = B.'*W2; % Newly needed for QB work
  QWk{2} = Q.'*W2; % Newly needed for QB work
  Im = speye(m);
  if ( d>3 )

    for k=4:d
      BWk{k-1} = B.'*reshape(w{k-1},n,n^(k-2));
      QWk{k-1} = Q.'*reshape(w{k-1},n,n^(k-2));

      b = -LyapProduct(N.',w{k-1},k-1); % Pre-compue all the L(N') terms
      
      % Now add all the terms from the 'B' sum by looping through the i and j
      for i=3:(k+1)/2 % i+j=k+2
        j   = k+2-i;
        tmp = BWk{i}.'*BWk{j};
        b   = b + 0.25*eta*i*j*(vec(tmp) + vec(tmp.'));
      end

      if (mod(k,2)==0) % k+2 is even
        i   = (k+2)/2; 
        j   = i;
        tmp = BWk{i}.'*BWk{j};
        b   = b + 0.25*eta*i*j*vec(tmp);
      end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%
      % Now add the linear Q terms by iterating through the first sum
      for i=2:k-1 % i+j=k+1
        j   = k+1-i;
        
        tmp = kron(QWk{j},BWk{i}); 
        
        %
        b   = b + 0.25*eta*i*j*kron(speye(n^k),vec(Im).')*vec(tmp);
        % The first part " kron(speye(n^k),vec(Im).') " can probably be made quickly w/ a special function; it is some almost `rectangular diagonal` thing, not a permutation matrix tho; look to how perfect shuffle is created
        
      end
      
      %%
      % Now add the quadratic Q terms by iterating through the second sum
      for i=2:k-2 % i+j=k
        j   = k-i;
        
        tmp = kron(speye(n),vec(Im).') ...
        * kron(vec(QWk{j}).' , kron(QWk{i},Im)) ...
        * kron(speye(n^(j-1)),kron(perfectShuffle(n^(i-1),n*m),Im))...
        * kron(speye(n^(k-1)),vec(Im));
        %
        b   = b + 0.125*eta*i*j*vec(tmp);
      end

 
      [w{k}] = KroneckerSumSolver(Acell(1:k),b,k);

      [w{k}] = kronMonomialSymmetrize(w{k},n,k);
    end
  end

end
