function fi = mphinterp_uy(x,y)

load([pwd,'\uy.mat'],'f');

fi = f(x,y);