function f = u2(c,n,gamma,b,epsilon,xgrid,ugrid)

%global gamma b epsilon

if gamma == 1
    f = interp1(xgrid,ugrid,c,'linear') - b*(n.^epsilon)/epsilon;
else
    f = (c.^(1-gamma))/(1-gamma)-b*(n.^epsilon)/epsilon;
end