function R3 = HestonCall(S,K,v,r,t,kp,et,sm,rh,ld,trunc,greek)
if nargin < 11, trunc = 100; greek = 1; end

I = sqrt(-1);
x = log(S) + r*t;

switch greek
    case 1 % Price
        inte = @(w) (exp(-I*w*x).*(K.^(1+I*w)).*Ffun(t,w,v,kp,ld,et,sm,rh)...
            ./(I*w-w.^2));
    case 2 % delta
        inte = @(w) (-I*w./S.*exp(-I*w*x).*(K.^(1+I*w))...
            .*Ffun(t,w,v,kp,ld,et,sm,rh)./(I*w-w.^2));
    case 3 % gamma
        inte = @(w) (-(w.^2)./(S^2).*exp(-I*w*x).*(K.^(1+I*w))...
            .*Ffun(t,w,v,kp,ld,et,sm,rh)./(I*w-w.^2));
    case 4 % vega
        inte = @(w) (exp(-I*w*x).*(K.^(1+I*w))...
            .*Ffun(t,w,v,kp,ld,et,sm,rh,1)./(I*w-w.^2));
end
R3t = exp(-r*t)*quadgk(inte,-trunc+2i,trunc+2i)/(2*pi);
R3 = real(R3t);
end


function R1 = Ffun(t,w,v,kp,ld,et,sm,rh,Indi)

if nargin < 9
Indi = 0;
end
        
I = sqrt(-1);

d = sqrt((rh*sm*I*w+kp+ld).^2+(w.^2-I*w).*(sm^2));
g = (kp+ld+rh*sm*I*w+d)./(kp+ld+rh*sm*I*w-d);
D = ((kp+ld+rh*sm*I*w+d)/(sm^2)).*(1-exp(d*t))./(1-g.*exp(d*t));
C = kp*et*((kp+ld+rh*sm*I*w+d)*t-2*log((1-g.*exp(d*t))./(1-g)))/(sm^2);

if Indi == 0
    R1 = exp(C+D.*v);
else R1 = D.*exp(C+D.*v);
end

end

        