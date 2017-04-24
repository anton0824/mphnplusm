% get the G matrix for the repeating part by the quadratic convergent
% algorithm

function G = Gmatrix(A0,A1,A2,epsilon)

maxiterate = 100000;

i = 1;
D = 1000;
H = (-A1)^(-1)*A0;
L = (-A1)^(-1)*A2;
G = L;
T = H;

while( D > epsilon )
    if( i > maxiterate )
        disp ('The algorithm does not converge!');
        break;
    end

    U = H*L+L*H;
    M = H^2;
    S = size(U);
    E = eye(S(1));
    H = (E-U)^(-1)*M;
    M = L^2;
    L = (E-U)^(-1)*M;
    G = G+T*L;
    T = T*H;

    S = size(G);
    I = ones(S(1),1);
    error = I-G*I;

    diff = 0;
    for k = 1:S
        diff = max(diff,abs(error(k)));
    end
    D = diff;
end
