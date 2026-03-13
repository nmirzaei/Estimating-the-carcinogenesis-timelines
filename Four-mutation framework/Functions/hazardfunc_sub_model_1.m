function y = hazardfunc_sub_model_1(t,x,p,N0)
    muN = p(1); mu1 = p(2); mu2 = p(3); mu3 = p(4); alpha1=p(5); alpha2=p(6); alpha3=p(7); beta1=p(8); beta2=p(9); beta3=p(10);
    y = zeros(2,1);

    y(1) = (-muN*N0*x(1));
    y(2) = 0; 
end
