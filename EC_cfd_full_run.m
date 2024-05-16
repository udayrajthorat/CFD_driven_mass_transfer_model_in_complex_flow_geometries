%%
tic
close all
clear all
 
rho = 1000;
mu = 0.001;
U0 = 1; 

[x0,y0,delta,xv,yv,x_,y_,u_,v_,muT_,ux_,vx_,muTx_,uy_,vy_,muTy_,up1,yp1] = EC_cfd_full(rho,mu,U0);

save output0.mat

elapsed_time = toc;

disp(['Elapsed time: ' num2str(elapsed_time) ' seconds']);

return

%%

load output0.mat

pseudo_periodicity = 0;

if pseudo_periodicity == 1
    
    % Define the index to split the geometry 

    Si= 20;   % Split index to cut the geometry
  
    xw_split= x0(:,1:Si);
    yw_split = y0(:,1:Si);
 
    xv_split = xv(1:Si,:);
    yv_split = yv(1:Si,:);
    
    wall_split_fname = 'wall_split.txt';
    vsl_split_fname = 'vsl_split.txt';
    
    write2txt(xw_split,yw_split,wall_split_fname);
    write2txt(xv_split,yv_split,vsl_split_fname);
       
    % 
    x_s = x_(1:Si,:);
    y_s = y_(1:Si,:);
    
    u_s = u_(1:Si,:);
    v_s = v_(1:Si,:);
    nuT_s = muT_/rho;
    nuT_s = nuT_s(1:Si,:);
    
    ux_s = ux_(1:Si,:);
    vx_s = vx_(1:Si,:);
    nuTx_s = muTx_/rho;
    nuTx_s = nuTx_s(1:Si,:);
    
    uy_s = uy_(1:Si,:);
    vy_s = vy_(1:Si,:);
    nuTy_s = muTy_/rho;
    nuTy_s = nuTy_s(1:Si,:);
    
    x_1_s = reshape(x_s,[],1);
    y_1_s = reshape(y_s,[],1);

    u_1_s = reshape(u_s,[],1);
    v_1_s = reshape(v_s,[],1);
    nuT_1_s = reshape(nuT_s,[],1);

    ux_1_s = reshape(ux_s,[],1);
    vx_1_s = reshape(vx_s,[],1);
    nuTx_1_s = reshape(nuTx_s,[],1);

    uy_1_s = reshape(uy_s,[],1);
    vy_1_s = reshape(vy_s,[],1);
    nuTy_1_s = reshape(nuTy_s,[],1);
    
    fu_s = scatteredInterpolant(x_1_s,y_1_s,u_1_s);
    fv_s = scatteredInterpolant(x_1_s,y_1_s,v_1_s);
    fnuT_s = scatteredInterpolant(x_1_s,y_1_s,nuT_1_s);
    
    fux_s = scatteredInterpolant(x_1_s,y_1_s,ux_1_s);
    fvx_s = scatteredInterpolant(x_1_s,y_1_s,vx_1_s);
    fnuTx_s = scatteredInterpolant(x_1_s,y_1_s,nuTx_1_s);
    
    fuy_s = scatteredInterpolant(x_1_s,y_1_s,uy_1_s);
    fvy_s = scatteredInterpolant(x_1_s,y_1_s,vy_1_s);
    fnuTy_s = scatteredInterpolant(x_1_s,y_1_s,nuTy_1_s);
    
    f = fu_s;
    save('u_s.mat','f');

    f = fv_s;
    save('v_s.mat','f');

    f = fnuT_s;
    save('nuT_s.mat','f');

    f = fux_s;
    save('ux_s.mat','f');

    f = fvx_s;
    save('vx_s.mat','f');

    f = fnuTx_s;
    save('nuTx_s.mat','f');

    f = fuy_s;
    save('uy_s.mat','f');

    f = fvy_s;
    save('vy_s.mat','f');

    f = fnuTy_s;
    save('nuTy_s.mat','f');
    
    save output1.mat
    
    figure 
    surf(x_s,y_s,u_s)
    xlabel('Length of mass transfer model');
    ylabel('Viscous sublayer thickness');
    zlabel('u component of velocity');
    saveas(gcf, 'u component of velocity.png');

    figure 
    surf(x_s,y_s,v_s)
    xlabel('Length of mass transfer model');
    ylabel('Viscous sublayer thickness');
    zlabel('v component of velocity');
    saveas(gcf, 'v component of velocity.png');

    figure 
    surf(x_s,y_s,nuT_s)
    xlabel('Length of mass transfer model');
    ylabel('Viscous sublayer thickness');
    zlabel('v component of velocity');
    saveas(gcf, 'Turbulent diffusivity profile.png');

else 

    wall_fname = 'wall.txt';
    vsl_fname = 'vsl.txt';

    write2txt(x0,y0,wall_fname);
    write2txt(xv,yv,vsl_fname);

    x_1 = reshape(x_,[],1);
    y_1 = reshape(y_,[],1);

    u_1 = reshape(u_,[],1);
    v_1 = reshape(v_,[],1);
    nuT_1 = reshape(muT_/rho,[],1);

    ux_1 = reshape(ux_,[],1);
    vx_1 = reshape(vx_,[],1);
    nuTx_1 = reshape(muTx_/rho,[],1);

    uy_1 = reshape(uy_,[],1);
    vy_1 = reshape(vy_,[],1);
    nuTy_1 = reshape(muTy_/rho,[],1);

    fu = scatteredInterpolant(x_1,y_1,u_1);
    fv = scatteredInterpolant(x_1,y_1,v_1);
    fnuT = scatteredInterpolant(x_1,y_1,nuT_1);

    fux = scatteredInterpolant(x_1,y_1,ux_1);
    fvx = scatteredInterpolant(x_1,y_1,vx_1);
    fnuTx = scatteredInterpolant(x_1,y_1,nuTx_1);

    fuy = scatteredInterpolant(x_1,y_1,uy_1);
    fvy = scatteredInterpolant(x_1,y_1,vy_1);
    fnuTy = scatteredInterpolant(x_1,y_1,nuTy_1);

    f = fu;
    save('u.mat','f');

    f = fv;
    save('v.mat','f');

    f = fnuT;
    save('nuT.mat','f');

    f = fux;
    save('ux.mat','f');

    f = fvx;
    save('vx.mat','f');

    f = fnuTx;
    save('nuTx.mat','f');

    f = fuy;
    save('uy.mat','f');

    f = fvy;
    save('vy.mat','f');

    f = fnuTy;
    save('nuTy.mat','f');

    save output1.mat
end

return

%%

load output1.mat

% ScT = 0.85;

% Must run the line below once to open 2nd matlab window 
% An error will occur in the 1st window
% Run the command pwd in the 1st window
% Copy the result and switch to the 2nd window
% Run the command cd <> in the 2nd window
% Where <> is the copied contents from the 1st window
% Then run the line below again from the 1st window

% pipe_tds_fun(ScT,wall_fname,vsl_fname);

%%

% r0 = sqrt(x0.^2+y0.^2);
% 
% figure
% plot(x0,y0,'-k',xv,yv,'-.k')
% 
% figure
% plot3(x_,y_,u_,'ko')
% 
% figure
% plot3(x_,y_,muT_/rho,'ko')
% 
% figure
% plot3(x_,y_,u_,'ko')
% 
% figure
% plot3(x_,y_,v_,'ko')