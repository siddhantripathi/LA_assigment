function GaussElim(A,b,n)
    for i = 1 : n-1
        for j = i+1 : n
            factor = A(i,j) / A(j,j);
            for k = i : n
                A(i,k) = A(i,k) - (factor*A(i,j));
            end
            b(j) = b(j) - (factor*b(i));    
        end    
    end
    
    x(n)=b(n)/A(n,n);
    for i=n-1:-1:1
        sum=0;
        for j=i+1:n
            sum=sum+A(i,j)*x(j);
        end
        x(i)=(b(i)-sum)/A(i,i);
    end

    disp([A,b],"Reduced Matrix:");
    disp(x,"Solution(x):");
endfunction

function LU(A,n)
    L = eye(n, n);
    
    for i = 1 : n-1
        for j = i+1 : n
            factor=A(i,j)/A(i,i);
            L(i,j) = factor;
            for k=i:n
                A(j,k)=A(j,k)-factor*A(i,k);
            end
        end
    end
    
    disp(L,"L = ");
    disp(A,"U = ");
endfunction

function Inverse(A,n)
    aug = [A , eye(n,n)];
    for i = 1:n-1
        for j = i+1:n
            aug(j,i:2*n)=aug(j,i:2*n)-aug(j,i)/aug(i,i) * aug(i,i:2*n);
        end
    end
    
    for i = n:-1:2
        aug(1:i-1,:) = aug(1: i-1,:) - aug(1:i-1,i)/aug(i,i) * aug(i,:);
    end
    
    for i = 1:n
        aug(i,:) = aug(i,:)/aug(i,i);
    end
    
    B = aug(:,n+1:2*n);
    disp(B,"Inverse of the matrix = ");
endfunction

n = input("Enter the size of the matrix:");
A = x_matrix("Enter the elements of the coefficient matrix",zeros(n,n));
b = x_matrix("Enter the elements of the coefficient matrix",zeros(n,1));
disp([A b],"Augmented matrix:")
GaussElim(A,b,n)
LU(A,n)
Inverse(A,n)