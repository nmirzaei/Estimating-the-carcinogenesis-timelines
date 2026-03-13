function y = hazardfunc_sub_model_1(t,x,p,N0)
    muN = p(1); mu1 = p(2); mu2 = p(3); alpha1=p(4); alpha2=p(5); beta1=p(6); beta2=p(7);
    y = zeros(2,1);



    % f= 1 - (alpha1-beta1)/(alpha1-beta1*exp((beta1-alpha1)*t));
    % fprime = beta1*((alpha1-beta1)^2)*exp((beta1-alpha1)*t)/((alpha1-beta1*exp((beta1-alpha1)*t))^2);
    


    y(1) = -muN*N0*x(1);
    y(2) = 0; %Cumulative hazard. To get hazard rate you take the derivative
end
