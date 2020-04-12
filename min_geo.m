function [ L ] = min_geo( X,Y,Z ,V , start,ending,change_stop,perc_stop,max_iter)

[rows columns] = size(X);

D0=inf(max(rows,columns));
D0 = D0(1:rows,1:columns);
D0(start(1),start(2));
D0(start(1),start(2))=0;

[D iter]=fm(X,Y,Z,D0,V,change_stop,perc_stop,max_iter);

row = ending(1);
column = ending(2);
L(:,1) = [X(row,column),Y(row,column),Z(row,column)];
path_size = 2;


prev_r1 = 0 ;
prev_c1 = 0;
prev_r2 = 0 ;
prev_c2 = 0;

while (~(row == start(1) && column==start(2)))
    
    min_r = 0;
    min_c = 0;
    min_val = 9E99;
    for r=-1:1
        for c=-1:1
            if (r~=0 && c~=0)
              row_ind = row+r;
              col_ind = column+c;
              if ((row_ind>0) && (col_ind>0) &&( row_ind<rows+1) && (col_ind<columns+1) && (row_ind~=prev_r2) && (col_ind ~= prev_c2))
                  val = D(row_ind,col_ind);
                  if (val<min_val)
                      min_val = val;
                      min_r = row_ind;
                      min_c = col_ind;
                  end
              end
            end
        end
    end
    prev_r2 = prev_r1;
    prev_c2 =  prev_c1;
    prev_r1 = min_r;
    prev_c1 = min_c;
    row = min_r;
    column = min_c;
    if (min_r==0 && min_c==0)
        row = prev_r1-(prev_r2-prev_r1);
        column = prev_c1-(prev_c2-prev_c1);
    else
      L(:,path_size) = [X(min_r,min_c),Y(min_r,min_c),Z(min_r,min_c)];
      path_size = path_size +1;    
    end
end

end

