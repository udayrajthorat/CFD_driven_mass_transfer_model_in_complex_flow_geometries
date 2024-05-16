function fi = mphinterp_vy(x,y)

load([pwd,'\vy.mat'],'f');

fi = f(x,y);