function [x] = KroneckerSumSolver(A,b,degree,M)
%KroneckerSumSolver Provides an efficient solution for a Kronecker sum system
%
%  Implements an N-Way version of the Bartels-Stewart algorithm for a 
%  special Kronecker sum system (the special case all matrices are the same)
%  with an additional
%
%  (kron(A{d},eye(n^(d-1))) + ... + kron(eye(n^(d-1),A{1}))+diag(M))x = b
%
%  where A is a cell array containing d matrices each of size (n,n) 
%  and b and x are size (n^d,1).  We assume that Schur factorizations have 
%  already been applied so that A{k} is upper triangular, for k=1:d.
%
%  This is consistent with KroneckerSumSolver that was developed for the
%  QQR problem.
%%

  n = size(A{1},1);

  
  %  As a simplification, we assume all A{i} are the same size.  We furthermore
  %  assume they are all the same (so we can use the kroneckerLeft function as is)
  %  Both of these can be easily relaxed if an application arises.
  [U,T] = schur(A{1},'complex');

  b = kroneckerLeft(U',b);

  L = length(A);
  for l=1:L
    A{l} = T;
  end
  
  if ( nargin<4 )
    UMU = spzeros(n);
  else
    UMU = U'*M*U;
  end

  
  switch degree

    case 2    %  The degree=2 case...

      X = zeros(n,n^2);
      B = reshape(b,n,n^2);

      for j3=n:-1:1
        j3Range = j3+1:n;
        for j2=n:-1:1
          j2Range = j2+1:n;

          At = A{1} + (A{2}(j2,j2)+A{3}(j3,j3))*eye(n) + UMU;
          rhs = B(:,j2+(j3-1)*n);

          if (~isempty(j2Range))
            rhs = rhs - X(:,j2Range+(j3-1)*n)*A{2}(j2,j2Range).';
          end

          if (~isempty(j3Range))
            rhs = rhs - X(:,j2+(j3Range-1)*n)*A{3}(j3,j3Range).';
          end

          X(:,j2+(j3-1)*n) = At\rhs;
        end
      end

      x = X(:);

    case 3    %  The degree=3 case...

      X = zeros(n,n^3);
      B = reshape(b,n,n^3);

      for j4=n:-1:1
        j4Range = j4+1:n;
        for j3=n:-1:1
          j3Range = j3+1:n;
          for j2=n:-1:1
            j2Range = j2+1:n;

            At = A{1} + (A{2}(j2,j2)+A{3}(j3,j3)+A{4}(j4,j4))*eye(n) + UMU;
            rhs = B(:,j2+(j3-1)*n+(j4-1)*n^2);

            if (~isempty(j2Range))
              rhs = rhs - X(:,j2Range+(j3-1)*n+(j4-1)*n^2)*A{2}(j2,j2Range).';
            end

            if (~isempty(j3Range))
              rhs = rhs - X(:,j2+(j3Range-1)*n+(j4-1)*n^2)*A{3}(j3,j3Range).';
            end

            if (~isempty(j4Range))
              rhs = rhs - X(:,j2+(j3-1)*n+(j4Range-1)*n^2)*A{4}(j4,j4Range).';
            end

            X(:,j2+(j3-1)*n+(j4-1)*n^2) = At\rhs;
          end
        end
      end  

      x = X(:);

    case 4 % the degree=4 case

      X = zeros(n,n^4);
      B = reshape(b,n,n^4);

      for j5=n:-1:1
        j5Range = j5+1:n;
        for j4=n:-1:1
          j4Range = j4+1:n;
          for j3=n:-1:1
            j3Range = j3+1:n;
            for j2=n:-1:1
              j2Range = j2+1:n;

              At = A{1} + (A{2}(j2,j2)+A{3}(j3,j3)+A{4}(j4,j4)+A{5}(j5,j5))*eye(n) + UMU;
              rhs = B(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3);

              if (~isempty(j2Range))
                rhs = rhs - X(:,j2Range+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3)*A{2}(j2,j2Range).';
              end

              if (~isempty(j3Range))
                rhs = rhs - X(:,j2+(j3Range-1)*n+(j4-1)*n^2+(j5-1)*n^3)*A{3}(j3,j3Range).';
              end

              if (~isempty(j4Range))
                rhs = rhs - X(:,j2+(j3-1)*n+(j4Range-1)*n^2+(j5-1)*n^3)*A{4}(j4,j4Range).';
              end

              if (~isempty(j5Range))
                rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5Range-1)*n^3)*A{5}(j5,j5Range).';
              end

              X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3) = At\rhs;
            end
          end
        end
      end   

      x = X(:);
      
    case 5 % the degree=5 case

      X = zeros(n,n^5);
      B = reshape(b,n,n^5);

      for j6=n:-1:1
        j6Range = j6+1:n;
        for j5=n:-1:1
          j5Range = j5+1:n;
          for j4=n:-1:1
            j4Range = j4+1:n;
            for j3=n:-1:1
              j3Range = j3+1:n;
              for j2=n:-1:1
                j2Range = j2+1:n;

                At = A{1} + (A{2}(j2,j2)+A{3}(j3,j3)+A{4}(j4,j4)+A{5}(j5,j5)+A{6}(j6,j6))*eye(n) + UMU;
                rhs = B(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4);

                if (~isempty(j2Range))
                  rhs = rhs - X(:,j2Range+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4)*A{2}(j2,j2Range).';
                end

                if (~isempty(j3Range))
                  rhs = rhs - X(:,j2+(j3Range-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4)*A{3}(j3,j3Range).';
                end

                if (~isempty(j4Range))
                  rhs = rhs - X(:,j2+(j3-1)*n+(j4Range-1)*n^2+(j5-1)*n^3+(j6-1)*n^4)*A{4}(j4,j4Range).';
                end

                if (~isempty(j5Range))
                  rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5Range-1)*n^3+(j6-1)*n^4)*A{5}(j5,j5Range).';
                end

                if (~isempty(j6Range))
                  rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6Range-1)*n^4)*A{6}(j6,j6Range).';
                end

                X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4) = At\rhs;
              end
            end
          end
        end   
      end
      
      x = X(:);

    case 6 % the degree=6 case

      X = zeros(n,n^6);
      B = reshape(b,n,n^6);

      for j7=n:-1:1
        j7Range = j7+1:n;
        for j6=n:-1:1
          j6Range = j6+1:n;
          for j5=n:-1:1
            j5Range = j5+1:n;
            for j4=n:-1:1
              j4Range = j4+1:n;
              for j3=n:-1:1
                j3Range = j3+1:n;
                for j2=n:-1:1
                  j2Range = j2+1:n;

                  At = A{1} + (A{2}(j2,j2)+A{3}(j3,j3)+A{4}(j4,j4)+A{5}(j5,j5)+A{6}(j6,j6)+A{7}(j7,j7))*eye(n) + UMU;
                  rhs = B(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5);

                  if (~isempty(j2Range))
                    rhs = rhs - X(:,j2Range+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5)*A{2}(j2,j2Range).';
                  end

                  if (~isempty(j3Range))
                    rhs = rhs - X(:,j2+(j3Range-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5)*A{3}(j3,j3Range).';
                  end

                  if (~isempty(j4Range))
                    rhs = rhs - X(:,j2+(j3-1)*n+(j4Range-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5)*A{4}(j4,j4Range).';
                  end

                  if (~isempty(j5Range))
                    rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5Range-1)*n^3+(j6-1)*n^4+(j7-1)*n^5)*A{5}(j5,j5Range).';
                  end

                  if (~isempty(j6Range))
                    rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6Range-1)*n^4+(j7-1)*n^5)*A{6}(j6,j6Range).';
                  end

                  if (~isempty(j7Range))
                    rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7Range-1)*n^5)*A{7}(j7,j7Range).';
                  end

                  X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5) = At\rhs;
                end
              end
            end
          end
        end   
      end
      
      x = X(:);

    case 7 % the degree=7 case

      X = zeros(n,n^7);
      B = reshape(b,n,n^7);

      for j8=n:-1:1
        j8Range = j8+1:n;
        for j7=n:-1:1
          j7Range = j7+1:n;
          for j6=n:-1:1
            j6Range = j6+1:n;
            for j5=n:-1:1
              j5Range = j5+1:n;
              for j4=n:-1:1
                j4Range = j4+1:n;
                for j3=n:-1:1
                  j3Range = j3+1:n;
                  for j2=n:-1:1
                    j2Range = j2+1:n;

                    At = A{1} + (A{2}(j2,j2)+A{3}(j3,j3)+A{4}(j4,j4)+A{5}(j5,j5)+A{6}(j6,j6)+A{7}(j7,j7)+A{8}(j8,j8))*eye(n) + UMU;
                    rhs = B(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5+(j8-1)*n^6);

                    if (~isempty(j2Range))
                      rhs = rhs - X(:,j2Range+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5+(j8-1)*n^6)*A{2}(j2,j2Range).';
                    end

                    if (~isempty(j3Range))
                      rhs = rhs - X(:,j2+(j3Range-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5+(j8-1)*n^6)*A{3}(j3,j3Range).';
                    end

                    if (~isempty(j4Range))
                      rhs = rhs - X(:,j2+(j3-1)*n+(j4Range-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5+(j8-1)*n^6)*A{4}(j4,j4Range).';
                    end

                    if (~isempty(j5Range))
                      rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5Range-1)*n^3+(j6-1)*n^4+(j7-1)*n^5+(j8-1)*n^6)*A{5}(j5,j5Range).';
                    end

                    if (~isempty(j6Range))
                      rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6Range-1)*n^4+(j7-1)*n^5+(j8-1)*n^6)*A{6}(j6,j6Range).';
                    end

                    if (~isempty(j7Range))
                      rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7Range-1)*n^5+(j8-1)*n^6)*A{7}(j7,j7Range).';
                    end

                    if (~isempty(j8Range))
                      rhs = rhs - X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5+(j8Range-1)*n^6)*A{8}(j8,j8Range).';
                    end

                    X(:,j2+(j3-1)*n+(j4-1)*n^2+(j5-1)*n^3+(j6-1)*n^4+(j7-1)*n^5+(j8-1)*n^6) = At\rhs;
                  end
                end
              end
            end
          end   
        end
      end
      
      x = X(:);

    otherwise
      error('not yet implemented')
  end

  x = kroneckerLeft(U,x);
end