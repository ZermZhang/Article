%源项（系数）的表达式
function result = source_a(x,y)

if (x>(90/1024) && x<(934/1024))
    result = 100000;
else
    result = (2 + 1.8*sin(2*pi*x/epsilon))./(2+1.8*cos(2*pi*y/epsilon))...
    +(2+1.8*sin(2*pi*y/epsilon))./(2+1.8*sin(2*pi*x/epsilon));
end