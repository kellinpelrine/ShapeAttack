function f = u(c,n,gamma,b,epsilon)

%global gamma b epsilon

if gamma == 1
    f = log(c)-b*(n.^epsilon)/epsilon;
else
    f = (c.^(1-gamma))/(1-gamma)-b*(n.^epsilon)/epsilon;
end