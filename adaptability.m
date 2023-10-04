function [z] = adaptability(x,y)
if (x>=-3 && x<=-1.5 && y>=-3 && y<=3) || (x>=-3 && x<=3 && y>=-3 && y<=-1.5) || (x>=1.5 && x<=3 && y>=1.5 && y<=3)
    z = power(sin(pi*x/1.5),2)* power(sin(pi*y/1.5),2)*exp((x+y)/5/pi);
else
    z = -power(sin(pi*x/1.5),2)* power(sin(pi*y/1.5),2)*exp((x+y)/5/pi);
end
end

