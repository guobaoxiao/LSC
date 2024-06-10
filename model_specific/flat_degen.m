function r = flat_degen(A)
        
dis = pdist(A');
if min(dis)<eps
    r = 1;
else
    r = 0;
end

end
