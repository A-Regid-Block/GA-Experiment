x = linspace(-6,6,2000);
y = linspace(-6,6,2000);

for i = 1:2000
    for j = 1:2000
        z(i,j) = adaptability(x(i),y(j));
    end
end

mesh(x,y,z);
grid on;