function [F, G] = solve_diophantin(A, B, D, delay) 
    % A*F +  q^-d * B*G = D
    syms q; 
    A_poly = poly2sym(A,q);
    B_poly = poly2sym(B,q);
    delay_poly = poly2sym([1, zeros(1, delay)], q);
    D_poly = poly2sym(D,q);
    
    deg_F = delay - 1;
    deg_G = polynomialDegree(A_poly) - 1;
    
    F = sym("f", [1, deg_F]);
    F = [1, F];
    G = sym("g", [1, deg_G+1]);
    
    F_poly = poly2sym(F,q);
    G_poly = poly2sym(G,q);
    
    LHS = A_poly * F_poly + delay_poly * B_poly * G_poly;
    RHS = D_poly;
    assert(polynomialDegree(LHS,q)==polynomialDegree(RHS,q), "LHS and RHS are not in same degree");
    
    LHS_coef = coeffs(LHS, q, 'all');
    RHS_coef = coeffs(RHS, q, 'all');
    
    
    if polynomialDegree(Bminus_poly * S_poly,q) < polynomialDegree(A_poly * Rprime_poly,q)
        results = double(vpa(struct2array(solve(LHS_coef(2:end)-RHS_coef(2:end)))));
    else
        results = double(vpa(struct2array(solve(LHS_coef-RHS_coef))));
    end
    
    Rprime_solved = results(1:length(Rprime)-1);
    S_solved = results(length(Rprime):length(Rprime) + length(S)-1);
    
    Rprime_solved = [1, Rprime_solved];
    Rprime_solved_poly = poly2sym(Rprime_solved,q);
    
    S_solved_poly = poly2sym(S_solved,q);
    R_solved_poly = Rprime_solved_poly * Bplus_poly;
    
    R_solved = double(vpa(coeffs(R_solved_poly, q, 'all')));
    

end