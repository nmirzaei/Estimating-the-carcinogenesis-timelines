function y = hazardfunc_sub_model_3(t,x,p,N0)
    muN = p(1); mu1 = p(2); mu2 = p(3); mu3 = p(4); alpha1=p(5); alpha2=p(6); alpha3=p(7); beta1=p(8); beta2=p(9); beta3=p(10);
    y = zeros(6,1);


    y(1) = (muN*N0*x(1)*(x(3)-1));
    y(2) = (-muN*N0*x(4)); 
    y(3) = (beta1 - (alpha1+beta1+mu1)*x(3)+mu1*x(3)*x(5)+alpha1*x(3)^2);
    y(4) = (-(alpha1+beta1+mu1)*x(4)+mu1*x(4)*x(5)+mu1*x(3)*x(6)+2*alpha1*x(3)*x(4));
    y(5) = (beta2 - (alpha2+beta2+mu2)*x(5)+alpha2*x(5)^2);
    y(6) = (-(alpha2+beta2+mu2)*x(6)+2*alpha2*x(5)*x(6));
end
