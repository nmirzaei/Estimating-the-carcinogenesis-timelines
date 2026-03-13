function y = timelines(PARAMS,N0,T,Age)

        TT = min(T,max(Age));

        x_ic1 = [1 PARAMS(1)*N0];
        [tp1,yp1] = ode15s(@(t,x,p) hazardfunc_sub_model_1(t,x,p,N0),0:TT, x_ic1, [], PARAMS);
    
        x_ic2 = [1 0 1 -PARAMS(2)];
        [tp2,yp2] = ode15s(@(t,x,p) hazardfunc_sub_model_2(t,x,p,N0),0:TT, x_ic2, [], PARAMS);    
    
        x_ic3 = [1 0 1 0 1 -PARAMS(3)];
        [tp3,yp3] = ode15s(@(t,x,p) hazardfunc_single_malignant_cells(t,x,p,N0),0:TT, x_ic3, [], PARAMS);   
    
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t_span  = 0:min(T,TT);
        dt = t_span(2) - t_span(1);
        P0_prime_1 = [0;yp1(2:min(T,TT)+1,2).*yp1(2:min(T,TT)+1,1)];
    
        P0_prime_2 = yp2(1:min(T,TT)+1,2).*yp2(1:min(T,TT)+1,1);
        P1_prime_2 = -yp2(1:min(T,TT)+1,4);
    
        P0_prime_3 = yp3(1:min(T,TT)+1,2).*yp3(1:min(T,TT)+1,1);
        P1_prime_3 = -yp3(1:min(T,TT)+1,4);
        P2_prime_3 = -yp3(1:min(T,TT)+1,6);
    
        % Normalize to get densities
        f_FSMC = P0_prime_1 / trapz(t_span, P0_prime_1); 
    
        f_SSMC_given_FSMC = P1_prime_2 / trapz(t_span, P1_prime_2);
        f_SSMC = conv(f_FSMC,f_SSMC_given_FSMC)*dt;
        f_SSMC = f_SSMC(1:length(t_span));
        f_SSMC = f_SSMC/trapz(t_span,f_SSMC);

        f_Malignant_given_SSMC = P2_prime_3 / trapz(t_span, P2_prime_3);
        f_Malignant = conv(f_SSMC,f_Malignant_given_SSMC)*dt;
        f_Malignant = f_Malignant(1:length(t_span));
        f_Malignant = f_Malignant/trapz(t_span,f_Malignant);

        E_FSMC = trapz(t_span, t_span .* f_FSMC');
        E_SSMC = trapz(t_span, t_span .* f_SSMC');
        E_Malignant = trapz(t_span, t_span .* f_Malignant');

        y = [E_FSMC E_SSMC E_Malignant];
end