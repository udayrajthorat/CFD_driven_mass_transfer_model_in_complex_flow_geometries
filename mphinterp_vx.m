function fi = mphinterp_vx(x,y)

load([pwd,'\vx.mat'],'f');

fi = f(x,y);