function inliers =  x84(res, n)
% n e' il numero minimo di inliers che ritorna

theta = 2.5; 
% attenzione: il theta di default di X84 e' 3.5
% cosi' e' piu' restrittivo, ma in linea con la regola di selezione degli
% inliers in LMEDS

location = nanmedian(res(:));
scale = theta/0.6745 * nanmedian(abs(res(:)-location));
inliers = abs(res-location) <= scale;

fprintf('x84 scale: %f \n', scale+location);

if length(inliers) < n
    % ritorna i primi n
    [~, i] = sort(res);
    inliers = i(1:min(numel(i),n));
    inliers =  sort(inliers);
    
end
end