function g=line_search(xl,xu)
    t = [xl:0.01:xu]
    tol = 10^(-12)
    r = (3-sqrt(5))/2
    c = xl+r*(xu-xl)
    d = xu-r*(xu-xl)
    fa = myfun(xu)
    fb = myfun(xl)
    fc = myfun(c)
    fd = myfun(d)
    while (xu-xl>tol)
        if (fc<fd) then
            xu = d
            fb = fd
            d = c
            fd = fc
            c = xl+r*(xu-xl)
            fc = myfun(c) 
        else
            xl = c
            fa = fc
            c = d
            fc = fd
            d = xu-r*(xu-xl)
            fd = myfun(d)
        end   
    end
    g=(xl+xu)/2
endfunction
