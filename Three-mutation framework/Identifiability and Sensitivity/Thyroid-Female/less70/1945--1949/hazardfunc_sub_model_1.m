function y = hazardfunc_sub_model_1(t,x,p,N0)
    muN = p(1); mu1 = p(2); mu2 = p(3); alpha1=p(4); alpha2=p(5); beta1=p(6); beta2=p(7);
    y = zeros(2,1);

    y(1) = -muN*N0*x(1);
    y(2) = 0;
end
