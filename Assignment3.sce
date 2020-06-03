function col_span(a,n)
    disp("Column Span:");
    disp("Converting to row reduced form");
    for i = 1 : n - 1
        for j = i + 1 : n
            a(j,:) = a(j,:) - (a(j,i) / a(i,i)) * a(i,:);
        end
        disp(a);
    end
    for i = 1 : n
        if(a(i,i) <> 0)
            a(i,:) = a(i,:) / a(i,i);
        end
    end
    disp(a)
    for i = 1 : n
        for j = i : n
            if(a(i,j) <> 0)
                disp('has a pivot element ',j,'column');
                break
            end
        end
    end
endfunction

function fundamental_subspaces(a)
    disp("Fundamental Spaces:");
    [m,n] = size(a);
    disp(m,'m is ');
    disp(n,'n is ');
    [v,pivot] = rref(a);
    disp(rref(a));
    disp(v);
    r = length(pivot);
    disp(r,'rank is ')
    cs = a(:,pivot);
    disp(cs,'Column Space is ');
    ns = kernel(a);
    disp(ns,'Null Space is ');
    rs = v(1:r,:)';
    disp(rs,'Row Space is ')
    lns = kernel(a');
    disp(lns,'Left Null Space is ');
endfunction

a = x_matrix("Enter matrix:",zeros(3,3));
disp(a,"a = ");
col_span(a,3);
fundamental_subspaces(a);