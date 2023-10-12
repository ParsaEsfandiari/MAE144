  % This code calculates the output of the diophantine method for the polynomials a & b & f
  % it also calculates the error of the reached answer with the desired f
  % but in this section we are looking for a strictly proper controller and we reach
  % the solution by adding poles to our desired f
  % Parsa Esfandiari, hw1, https://github.com/ParsaEsfandiari/
        
    b=RR_poly([-2,2,-5,5],1);
    a=RR_poly([-1,1,-3,3,-6,6],1);
    f=RR_poly([-1,-1,-3,-3,-6,-6,-60,-60 ,-80,-80],1); 
    [x,y]=RR_diophantine(a,b,f)
    test= trim(a*x+b*y)
    residual=norm(f-trim(a*x+b*y))
