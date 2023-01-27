function [v] = approxPastEnergyQB(A,N,B,Q,C,eta,d,verbose)
%  Calculates a polynomial approximation to the past energy function
%  for a quadratic bilinear system.
%
%  v = approxPastEnergy(A,N,B,Q,C,eta,d) 
%
%  Computes a degree d polynomial approximation to the past energy function 
%
%          E^-(x) = 1/2 ( v{2}'*kron(x,x) + ... + v{d}'*kron(.. x) )
%
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + N*kron(x,x) + Q*kron(x,u)
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation
%
%    A'*V2 + V2*A + V2*B*B'*V2 - eta*C'*C = 0.
%
%  and in the subsequent linear systems arising from the Past H-infinity 
%  Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%  Note that v{2} = vec(V2) = V2(:).  Details are in Section III.C of reference [1].
%
%  Requires functions from the KroneckerTools repository
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

  if (nargin<8)
    verbose = false;
  end

  n = size(A,1);   % A should be n-by-n
  m = size(B,2);   % B should be n-by-m
  p = size(C,1);   % C should be p-by-n
  
%   Q = sparse(n,n*m); % Hardcoded zero Q for now
%   Q(1,1) = 1; 

  % Create a vec function for readability
  vec = @(X) X(:);


  %% k=2 case
  R = eye(m);
  
  if ( eta~=0 )
    %  We multiply the ARE by -1 to put it into the standard form in icare,
    %  and know (-A,-B) is a controllable pair if (A,B) is.
    [V2] = icare(-A,-B,eta*(C.'*C),R);

    if ( isempty(V2) && verbose )
      warning('approxPastEnergy: icare couldn''t find stabilizing solution')
    end
    
    if ( isempty(V2) )
      if (verbose) 
        warning('approxPastEnergy: using matrix inverse')
      end
      [Yinf] = icare(A.',C.',B*B.',eye(p)/eta);
      if ( isempty(Yinf) )
        V2 = [];
      else
        V2 = inv(Yinf);  % yikes!!!
      end
    end

    if ( isempty(V2) )
      if (verbose), fprintf("trying the anti-solution\n"), end
      [V2] = icare(-A,B,eta*(C.'*C),R,'anti');
    end
    
    if ( isempty(V2) )
      if (verbose)
        warning('approxPastEnergy: icare couldn''t find stabilizing solution')
        fprintf('approxPastEnergy: using the hamiltonian\n')
      end
      [~,V2,~] = hamiltonian(-A,B,eta*(C.'*C),R,true);
      V2 = real(V2);
    end
     
    if ( isempty(V2) )
      error('Could not find a solution to the ARE, adjust the eta parameter')
    end
  
  else % eta==0
    %  This case is described in Section II.B of the paper and requires a
    %  matrix inverse to calculate E_c.
    [V2] = lyap(A,(B*B.'));
    V2 = inv(V2); % yikes!!!!!!!!

    %  To do: look at approximating this by [V2] = icare(-A,-B,eta*(C.'*C),R)
    %  with a small value of eta (and perhaps other choices for C.'*C)
    
  end
  
  %  Reshape the resulting quadratic coefficients
  v{2} = vec(V2);

  %% k=3 case
  BVk = cell(1,d-1);% Pre-compute B.'*Vi for all the i we need
  QVk = cell(1,d-1); % Pre-compute Q.'*Vi for all the i we need
  BVk{2} = B.'*V2; % Newly needed for QB work
  QVk{2} = Q.'*V2; % Newly needed for QB work
  Im = speye(m);
  if ( d>2 )
    % set up the generalized Lyapunov solver
    V2BB = V2*(B*B.');
    Acell = cell(1,d);
    for i=1:d
      Acell{i} = A.'+V2BB;
    end
    
    b = -LyapProduct(N.',v{2},2) ... 
        - 2*kron(speye(n^3),vec(Im).')*vec(kron(QVk{2},BVk{2})); % New QB term
    
    [v{3}] = KroneckerSumSolver(Acell(1:3),b,3);

    [v{3}] = kronMonomialSymmetrize(v{3},n,3);
  end
  
  %% k>3 case (up to d)
  if ( d>3 )
    for k=4:d
      BVk{k-1} = B.'*reshape(v{k-1},n,n^(k-2));
      QVk{k-1} = Q.'*reshape(v{k-1},n,n^(k-2));

      b = -LyapProduct(N.',v{k-1},k-1);

      for i=3:(k+1)/2
        j   = k+2-i;
        tmp = BVk{i}.'*BVk{j};
        b   = b - 0.25*i*j*(vec(tmp) + vec(tmp.'));
      end

      if (mod(k,2)==0) % k is even
        i   = (k+2)/2;
        j   = i;
        tmp = BVk{i}.'*BVk{j};
        b   = b - 0.25*i*j*vec(tmp);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%
      % Now add the linear Q terms by iterating through the first sum
      for i=2:k-1 % i+j=k+1
        j   = k+1-i;
        
        tmp = kron(QVk{j},BVk{i}); 
        
        %
        b   = b - 0.5*i*j*kron(speye(n^k),vec(Im).')*vec(tmp);
        % The first part " kron(speye(n^k),vec(Im).') " can probably be made quickly w/ a special function; it is some almost `rectangular diagonal` thing, not a permutation matrix tho; look to how perfect shuffle is created
        
      end
      
      %%
      % Now add the quadratic Q terms by iterating through the second sum
      for i=2:k-2 % i+j=k
        j   = k-i;
        
        tmp = kron(speye(n),vec(Im).') ...
        * kron(vec(QVk{j}).' , kron(QVk{i},Im)) ...
        * kron(speye(n^(j-1)),kron(perfectShuffle(n^(i-1),n*m),Im))...
        * kron(speye(n^(k-1)),vec(Im));
        %
        b   = b - 0.25*i*j*vec(tmp);
      end

      
      [v{k}] = KroneckerSumSolver(Acell(1:k),b,k);

      [v{k}] = kronMonomialSymmetrize(v{k},n,k);
    end
  end
  
end
