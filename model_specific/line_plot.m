function [] = line_plot(hyp)
    
x = 0:1:10;
y = (-hyp(1)*x - hyp(3))./hyp(2);

plot(x,y,'r-','LineWidth',2);

end