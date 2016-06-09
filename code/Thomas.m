function x = thomas(alpha, beta, gamma, f)

bbb(1)=beta(1);
fff(1)=f(1);
n=length(f);

for i=2:n
    mult=alpha(i)/bbb(i-1);
    bbb(i)=beta(i)-mult*gamma(i-1);
    fff(i)=f(i)-mult*fff(i-1);
end
x(n)=fff(n)/bbb(n);
for i=n-1:-1:1
    x(i)=(fff(i)-gamma(i)*x(i+1))/bbb(i);
end
    
