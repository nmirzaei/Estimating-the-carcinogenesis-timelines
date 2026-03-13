function y = hazardfunc_sub_model_1_brc(t,x,p,N0,Menarche_Age)
    muN = p(1); mu1 = p(2); mu2 = p(3); alpha1=p(4); alpha2=p(5); beta1=p(6); beta2=p(7);
    y = zeros(2,1);

    Delay =  double((t-Menarche_Age) >= 0);

    y(1) = Delay*(-muN*N0*x(1));
    y(2) = 0; 
end
