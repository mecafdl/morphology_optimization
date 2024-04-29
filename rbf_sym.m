function y = rbf_sym(x, a, b, c)

N = length(a);
y = 0.;
for i=1:N
    y = y + a(i) * exp(-b(i) * (x - c(i)).^2);
end

for i=1:N
    y = y - a(i) * exp(-b(i) * (x + c(i)).^2);
end



end