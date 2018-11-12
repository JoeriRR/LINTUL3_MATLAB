function [F] = GetF_gui_grid(m,scaling)
%generates the cost matrix E
%it assumes the fields are the same grid as the farm gui allocates the
%fields
nrows = round(sqrt(m));
if (mod(m,nrows) == 0)
    ncols = m/nrows;
else
    if(nrows^2<m)
        ncols = nrows+1;
    else
        ncols = nrows;
    end
end
E = zeros(m);
rowvec=[]; colvec=[];
for i = 1:nrows
    rowvec = [ones(1,ncols)*i rowvec];    
end
for i = 1:nrows
    colvec = [1:ncols colvec];    
end
midpoint = (max(rowvec)+min(rowvec))/2;
F = zeros(m,1);
for i = 1:m
   F(i) = sqrt((rowvec(i)-midpoint)^2+colvec(i)^2); 
end
F = scaling*F;
end

