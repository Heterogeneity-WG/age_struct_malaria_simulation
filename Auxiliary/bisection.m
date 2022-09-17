function  [r, ierror, a,b,k, ao, bo]  = bisection(myfunction, afirst, bfirst, max_steps, y_tol, x_tol, lprint)
% use the method of bisection to find a root

% Check that that neither end-point is a root
% and if f(a) and f(b) have the same sign, warning.
format long

ierror=0;
a=afirst;
b=bfirst;
ao = a;
bo = b;
fa=myfunction(a);
fb=myfunction(b);

header = ' Iter                c               f(c)                    a                  b               f(a)              f(b)';
if lprint; disp(header) ;end

if fa == 0
    r = a;
    return;
elseif fb == 0
    r = b;
    return;
elseif fa*fb>0
    ierror=1;
    r = a; k = 0;
    warning( 'myfunction(a) and myfunction(b) do not have opposite signs' );
    fa
    fb
    keyboard
    return;
end

% We will iterate max_steps times and if a root was not
% found after max_steps iterations, a warning will be thrown.

for k = 1:max_steps
    c = (a + b)/2;
    fc = myfunction(c);   
    if lprint; fprintf('%5.0f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f \n',k,c,fc,a, b,fa,fb); end
    % Check if we found a root or whether or not
    % we continue with:
    %    [a, c] if f(a) and f(c) have opposite signs, or
    %    [c, b] if f(c) and f(b) have opposite signs.
    if abs(fc)<y_tol
        r = c; 
        if fc*fa<0
          bo = b;
          b = c; 
          fb = fc;
        else
          ao = a;
          a = c;
          fa = fc;
        end
        if lprint; fprintf('%5.0f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f \n',k,c,fc,a, b,fa,fb); end
        if lprint; disp('** < ytol **'); end
        return;
    elseif fc*fa<0
        bo = b;
        b = c; 
        fb = fc;
    else
        ao = a;
        a = c; 
        fa = fc;
    end
    
%     If |b - a| < eps_step, check convergence
    if abs(b-a)<x_tol
        if fa<fb
            r = a;
            if lprint; fprintf('%5.0f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f \n',k,c,fc,a, b,fa,fb); end
            if lprint; disp('** < xtol **'); end
            return;
        else
            r = b; 
            if lprint; fprintf('%5.0f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f \n',k,c,fc,a, b,fa,fb); end
            if lprint; disp('** < xtol **'); end
            return;
        end
    end
    
end
ierror=2;
if lprint; disp('** > max_steps **'); end
if fa<fb
    r = a;
    if lprint; fprintf('%5.0f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f \n',k,c,fc,a, b,fa,fb); end
    return;
else
    r = b;
    if lprint; fprintf('%5.0f %18.16f %18.16f %18.16f %18.16f %18.16f %18.16f \n',k,c,fc,a, b,fa,fb); end
    return;
end
end