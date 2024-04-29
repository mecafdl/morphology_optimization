function y = rbf_norm(x, a, b, c)

N = length(a);
y = 0.;
y_sum = 0.;
for i=1:N
    y = y + a(i) * exp(-b(i) * (x - c(i)).^2);
    y_sum = y_sum +exp(-b(i) * (x - c(i)).^2);
end

y = y/y_sum;
end