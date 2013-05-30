function [Min, Max] = EECG(A, b, Size)
    x = zeros(b)
    r = b
    p = r
    D = eye(Size, Size)
    B0 = eye(Size, Size)
    B1 = eye(Size, Size)
    L = eye(Size, Size)
    L(1, 1) = p' * p
    D(1, 1) = norm(r)
    for i=1:Size
        rrp = (r' * r)
        alp = (r' * r) / (p' * A * p)
        r = r - alp * A * p
        bet = (r' * r) / rrp
        p = r + bet * p 
        if i<>Size then
            B0(i, i+1) = -1*bet
            B1(i+1, i) = -1*bet
            L(i+1, i+1) = p' * p
            D(i+1, i+1) = norm(r)
        end
    end
    T = inv(D)*B1*L*B0*inv(D)
    d = spec(T)
    disp(d)
//    [Min, Max] = SEE(d)
//    Min = min(d)
//    Max = max(d)
endfunction
