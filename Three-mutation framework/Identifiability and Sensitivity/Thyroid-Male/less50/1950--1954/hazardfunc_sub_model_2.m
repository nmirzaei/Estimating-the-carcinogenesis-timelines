function y = hazardfunc_sub_model_2(t,x,p,N0)
    muN = p(1); mu1 = p(2); mu2 = p(3); alpha1=p(4); alpha2=p(5); beta1=p(6); beta2=p(7);
    y = zeros(4,1);


    y(1) = muN*N0*x(1)*(x(3)-1);
    y(2) = -muN*N0*x(4); 
    y(3) = beta1 - (alpha1+beta1+mu1)*x(3)+alpha1*x(3)^2;
    y(4) = -(alpha1+beta1+mu1)*x(4)+2*alpha1*x(3)*x(4);
end
