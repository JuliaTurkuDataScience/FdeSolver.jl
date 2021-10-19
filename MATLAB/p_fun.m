function dy=p_fun(t,y,param)
p=param;
if t>1
    dy=1/gamma(2-p)*t^(1-p)-2/gamma(3-p)*(t-1)^(2-p);
else
    dy=1/gamma(2-p)*t^(1-p);
end

end
