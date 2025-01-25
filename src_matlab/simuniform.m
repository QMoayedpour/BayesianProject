function st = simuniform(m)
p = ones(m,1)/m;
st = sum(cumsum(p) < rand(1)) + 1;   
