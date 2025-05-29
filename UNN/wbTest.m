a3 = 1;
N = 100;

d1 = complex(randn(N,1),randn(N,1))/sqrt(2); 
d2 = complex(randn(N,1),randn(N,1))/sqrt(2); 

d3 = zeros(3*N-2,1);
for n = 1:N
    for m = 1:N
        for k = 1:N
            q = n + m - k - 1 + N;
            d3(q) = d3(q)+ 3/4*a3*d1(n)*d1(m)*conj(d2(k));
        end
    end
end

A = d1*d1.';
A = A(:,end:-1:1);
b = NaN(2*N-1,1);
for k = -N+1:N-1
    b(k+N) = sum(diag(A,k));
end
c = 3/4*a3*conv(b,conj(d2));
c = c(end:-1:1);

figure(1);
clf;
plot(real(d3));
hold on;
plot(real(c),'--');
q = floor(3*N/2) - 1;
plot([q q],[min(real(d3)) max(real(d3))],'r:');

