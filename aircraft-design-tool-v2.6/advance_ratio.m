function [mu] = advance_ratio(n,r,v,aoa)
mu = v * cosd(aoa)/ n / r;
end