load surface.mat

[rows columns] = size(X);

max_iter = 50;
perc_stop = 10;
change_stop =  100;

mkdir jpg
for i=1:5
    V_num = ceil(rand()*4);
    switch (V_num)
        case 1 
            V = V1;
        case 2 
            V = V2;
        case 3 
            V = V3;
        case 4
            V = V4;
    end
  i0 = ceil(rand()*rows);
  j0 = ceil(rand()*columns);
  i1 = ceil(rand()*rows);
  j1 = ceil(rand()*columns);
  L= min_geo( X,Y,Z , V , [i0,j0],[i1,j1],change_stop,perc_stop,max_iter);
  fig1=figure;
  surface(X,Y,Z,V); 
  axis image; 
  shading interp;
  hold on; 
  plot3(L(1,:), L(2,:), L(3,:), 'k'); 
  hold off;
  title(sprintf('Minimum geodisics (%d,%d)->(%d,%d) V=V%d',i0,j0,i1,j1,V_num),'FontSize',14,'FontWeight','bold');
  print(fig1, '-djpeg',['jpg/',sprintf('min_geo%d.jpg',i)]);
end

