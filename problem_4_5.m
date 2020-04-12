load surface.mat

max_iter = 50;
perc_stop = 10;
change_stop = 100;

mkdir jpg
[D iter]=fm(X,Y,Z,D1,V1,change_stop,perc_stop,max_iter);
fig1=figure;
surf(X,Y,Z,D); 
axis image; 
shading interp;
title('Fast marching D1/V1','FontSize',14,'FontWeight','bold');
print(fig1, '-djpeg',['jpg/','fm1.jpg']);
imwrite(D/max(max(D)),'jpg/im1.jpg','JPG','Quality' ,95);

[D iter]=fm(X,Y,Z,D2,V2,change_stop,perc_stop,max_iter);
fig2=figure;
surf(X,Y,Z,D); 
axis image; 
shading interp;
title('Fast marching D2/V2','FontSize',14,'FontWeight','bold');
print(fig2, '-djpeg',['jpg/','fm2.jpg']);
imwrite(D/max(max(D)),'jpg/im2.jpg','JPG','Quality' ,95);

[D iter]=fm(X,Y,Z,D3,V3,change_stop,perc_stop,max_iter);
fig3=figure;
surf(X,Y,Z,D); 
axis image; 
shading interp;
title('Fast marching D3/V3','FontSize',14,'FontWeight','bold');
print(fig3, '-djpeg',['jpg/','fm3.jpg']);
imwrite(D/max(max(D)),'jpg/im3.jpg','JPG','Quality' ,95);

[D iter]=fm(X,Y,Z,D4,V4,change_stop,perc_stop,max_iter);
fig4=figure;
surf(X,Y,Z,D); 
axis image; 
shading interp;
title('Fast marching D4/V4','FontSize',14,'FontWeight','bold');
print(fig4, '-djpeg',['jpg/','fm4.jpg']);
imwrite(D/max(max(D)),'jpg/im4.jpg','JPG','Quality' ,95);

