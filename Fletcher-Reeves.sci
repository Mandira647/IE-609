function [x_opt,f_opt,k] = conjugate_gradient (fob,g_fob,x0,tol_grad); //tol_grad is the tolerance eg 0.00001

    c0 = feval(x0,g_fob); //evaluate gradient at initial point
    k = 0;   

    if norm(c0) < tol_grad

        x_opt = x0;  //optimum point  
        f_opt = feval(x_opt,fob); //cost function value
    else 
        d0= -c0;  //search direction
        alfa0 = equal_interval_line_search(x0,d0,fob,0.5,1e-6); //line search (step size)
        x1= x0+ alfa0*d0;  
        c1 = feval(x1,g_fob); 

        while norm(c1) > tol_grad
            beta = (norm(c1)/norm(c0))^2;
            d1= -c1+beta*d0;
            alfa1 = equal_interval_line_search(x1,d1,fob,0.5,1e-6);
            x2= x1+alfa1*d1;
            c0=c1;
            c1= feval(x2,g_fob);
            d0=d1;
            x1=x2;
            k=k+1;
        end
    x_opt = x1;    
    f_opt = feval(x_opt,fob); 
    end
endfunction
