function [E] = GetE_gui_grid(m,scaling)
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
index = 1;
rowvec=[]; colvec=[];
for i = 1:nrows
    rowvec = [ones(1,ncols)*i rowvec];    
end
for i = 1:nrows
    colvec = [1:ncols colvec];    
end
for i = 1:m   
   row_init = rowvec(i);
   col_init = colvec(i);
   for j = 1:m       
       row = rowvec(j);
       col = colvec(j);       
       E(i,j) = sqrt((row-row_init)^2+(col-col_init)^2);       
   end
end
E = scaling*E;
end

