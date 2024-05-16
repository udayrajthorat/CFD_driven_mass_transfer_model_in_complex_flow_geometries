function [x0,y0,delta,xv,yv,x_,y_,u_,v_,muT_,ux_,vx_,muTx_,uy_,vy_,muTy_,up1,yp1] = EC_cfd_full(rho,mu,U0)

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('C:\Users\mnuyt\OneDrive - University of Leeds\PC050223\Paper 3 - Complex flow geometry\Livelink\Expansion constriction new - 091023');

model.label('EC_cfd_final.mph');

model.param.set('D1', '25.4[mm]', 'Expansion diameter');
model.param.set('L1', 'D1*5', 'Expansion length');
model.param.set('D2', 'D1/2', 'Constriction diameter');
model.param.set('L2', 'D2*10', 'Constriction length');
model.param.set('rho', '1000');
model.param.set('mu', '0.001');
model.param.set('U0', '1');
model.param.set('L3', 'D1/4');
model.param.set('mxe', '1500');
model.param.set('mxec', '300');
model.param.set('mxc', '1500');
model.param.set('mye', '175');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('size', {'L1' 'D1'});
model.component('comp1').geom('geom1').create('ls2', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls2').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('coord1', {'L1' '0'});
model.component('comp1').geom('geom1').feature('ls2').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('coord2', {'L1+2e-2' 'L3'});
model.component('comp1').geom('geom1').create('ls3', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls3').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord1', {'L1' 'D1'});
model.component('comp1').geom('geom1').feature('ls3').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord2', {'L1+2e-2' 'D1-L3'});
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'0.147' 'L3'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'L2' 'D2'});
model.component('comp1').geom('geom1').create('ls4', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls4').selection('vertex1').set('ls3(1)', 1);
model.component('comp1').geom('geom1').feature('ls4').selection('vertex2').set('ls2(1)', 1);
model.component('comp1').geom('geom1').create('ls5', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls5').selection('vertex1').set('r2(1)', 4);
model.component('comp1').geom('geom1').feature('ls5').selection('vertex2').set('r2(1)', 1);
model.component('comp1').geom('geom1').create('csol1', 'ConvertToSolid');
model.component('comp1').geom('geom1').feature('csol1').selection('input').set({'ls2' 'ls3' 'ls4' 'ls5'});
model.component('comp1').geom('geom1').feature('fin').set('repairtoltype', 'absolute');
model.component('comp1').geom('geom1').feature('fin').set('absrepairtol', '1e-9');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('x0', 'x');
model.component('comp1').variable('var1').set('y0', 'y');
model.component('comp1').variable('var1').set('nx0', 'nx');
model.component('comp1').variable('var1').set('ny0', 'ny');
model.component('comp1').variable('var1').selection.geom('geom1', 1);
model.component('comp1').variable('var1').selection.set([2 5 8]);

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Interpolation');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an1', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an2', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an3', 'Analytic');

model.component('comp1').physics.create('spf', 'TurbulentFlowSST', 'geom1');
model.component('comp1').physics('spf').create('inl1', 'InletBoundary', 1);
model.component('comp1').physics('spf').feature('inl1').selection.set([1]);
model.component('comp1').physics('spf').create('out1', 'OutletBoundary', 1);
model.component('comp1').physics('spf').feature('out1').selection.set([10]);

model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('map1').create('dis1', 'Distribution');
model.component('comp1').mesh('mesh1').feature('map1').create('dis2', 'Distribution');
model.component('comp1').mesh('mesh1').feature('map1').create('dis3', 'Distribution');
model.component('comp1').mesh('mesh1').feature('map1').create('dis4', 'Distribution');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis1').selection.set([2 3]);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').selection.set([1 4 7 10]);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis3').selection.set([5 6]);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis4').selection.set([8 9]);

model.component('comp1').view('view1').axis.set('xmin', -0.013191482052206993);
model.component('comp1').view('view1').axis.set('xmax', 0.2796543538570404);
model.component('comp1').view('view1').axis.set('ymin', -0.05547146126627922);
model.component('comp1').view('view1').axis.set('ymax', 0.08162516355514526);

model.component('comp1').material('mat1').label('Water, liquid');
model.component('comp1').material('mat1').set('family', 'water');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'273.15' '413.15' '1.3799566804-0.021224019151*T^1+1.3604562827E-4*T^2-4.6454090319E-7*T^3+8.9042735735E-10*T^4-9.0790692686E-13*T^5+3.8457331488E-16*T^6'; '413.15' '553.75' '0.00401235783-2.10746715E-5*T^1+3.85772275E-8*T^2-2.39730284E-11*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'273.15' '553.75' '12010.1471-80.4072879*T^1+0.309866854*T^2-5.38186884E-4*T^3+3.62536437E-7*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('smooth', 'contd1');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('pieces', {'273.15' '293.15' '0.000063092789034*T^3-0.060367639882855*T^2+18.9229382407066*T-950.704055329848'; '293.15' '373.15' '0.000010335053319*T^3-0.013395065634452*T^2+4.969288832655160*T+432.257114008512'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'273.15' '553.75' '-0.869083936+0.00894880345*T^1-1.58366345E-5*T^2+7.97543259E-9*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('table', {'273' '1403';  ...
'278' '1427';  ...
'283' '1447';  ...
'293' '1481';  ...
'303' '1507';  ...
'313' '1526';  ...
'323' '1541';  ...
'333' '1552';  ...
'343' '1555';  ...
'353' '1555';  ...
'363' '1550';  ...
'373' '1543'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('interp', 'piecewisecubic');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', {'m/s'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(T)*d(rho(T),T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'gamma_w');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '1+(T/Cp(T))*(alpha_p(T)*cs(T))^2');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', '1');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('funcname', 'muB');
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('expr', '2.79*eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an3').set('plotargs', {'T' '273.15' '553.75'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(T)' '0' '0' '0' 'alpha_p(T)' '0' '0' '0' 'alpha_p(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', 'gamma_w(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');

model.component('comp1').physics('spf').prop('TurbulenceModelProperty').set('WallTreatment', 'LowReynoldsNumber');
model.component('comp1').physics('spf').feature('fp1').set('rho_mat', 'userdef');
model.component('comp1').physics('spf').feature('fp1').set('rho', 'rho');
model.component('comp1').physics('spf').feature('fp1').set('mu_mat', 'userdef');
model.component('comp1').physics('spf').feature('fp1').set('mu', 'mu');
model.component('comp1').physics('spf').feature('inl1').set('BoundaryCondition', 'FullyDevelopedFlow');
model.component('comp1').physics('spf').feature('inl1').set('Uavfdf', 'U0');

model.component('comp1').mesh('mesh1').feature('size').set('hauto', 3);
model.component('comp1').mesh('mesh1').feature('size').set('table', 'cfd');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis1').set('numelem', 'mxe');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').set('type', 'predefined');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').set('elemcount', 'mye');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').set('elemratio', 1000);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').set('growthrate', 'exponential');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis2').set('symmetric', true);
model.component('comp1').mesh('mesh1').feature('map1').feature('dis3').set('numelem', 'mxec');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis4').set('numelem', 'mxc');
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('wdi', 'WallDistanceInitialization');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').create('su1', 'StoreSolution');
model.sol('sol1').create('st2', 'StudyStep');
model.sol('sol1').create('v2', 'Variables');
model.sol('sol1').create('s2', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').feature('s2').create('se1', 'Segregated');
model.sol('sol1').feature('s2').create('d1', 'Direct');
model.sol('sol1').feature('s2').create('d2', 'Direct');
model.sol('sol1').feature('s2').create('i1', 'Iterative');
model.sol('sol1').feature('s2').create('i2', 'Iterative');
model.sol('sol1').feature('s2').feature('se1').create('ss1', 'SegregatedStep');
model.sol('sol1').feature('s2').feature('se1').create('ss2', 'SegregatedStep');
model.sol('sol1').feature('s2').feature('se1').create('ll1', 'LowerLimit');
model.sol('sol1').feature('s2').feature('se1').feature.remove('ssDef');
model.sol('sol1').feature('s2').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').create('sc1', 'SCGS');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').create('sc1', 'SCGS');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s2').feature('i2').create('mg1', 'Multigrid');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').create('sl1', 'SORLine');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').create('sl1', 'SORLine');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol1').feature('s2').feature.remove('fcDef');

model.result.dataset.create('edg1', 'Edge2D');
model.result.dataset('edg1').selection.set([2 3 5 6 8 9]);
model.result.dataset.remove('dset2');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result.create('pg3', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg2').create('con1', 'Contour');
model.result('pg2').feature('con1').set('expr', 'p');
model.result('pg3').create('line1', 'Line');
model.result('pg3').feature('line1').set('expr', 'spf.Delta_wPlus');
model.result('pg3').feature('line1').create('hght1', 'Height');
model.result('pg3').feature('line1').feature('hght1').set('expr', 'spf.WRHeightExpr');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Wall Distance Initialization');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').set('stol', 1.0E-6);
model.sol('sol1').feature('s1').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'd1');
model.sol('sol1').feature('s1').feature('fc1').set('initstep', 0.01);
model.sol('sol1').feature('s1').feature('fc1').set('minstep', 1.0E-6);
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', 50);
model.sol('sol1').feature('s1').feature('d1').label('Direct, wall distance (spf)');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s1').feature('i1').label('AMG, wall distance (spf)');
model.sol('sol1').feature('s1').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s1').feature('i1').set('maxlinit', 1000);
model.sol('sol1').feature('s1').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('linemethod', 'uncoupled');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('linemethod', 'uncoupled');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('su1').label('Solution Store 1.1');
model.sol('sol1').feature('st2').label('Compile Equations: Stationary');
model.sol('sol1').feature('st2').set('studystep', 'stat');
model.sol('sol1').feature('v2').label('Dependent Variables 2.1');
model.sol('sol1').feature('v2').set('initmethod', 'sol');
model.sol('sol1').feature('v2').set('initsol', 'sol1');
model.sol('sol1').feature('v2').set('solnum', 'auto');
model.sol('sol1').feature('v2').set('notsolmethod', 'sol');
model.sol('sol1').feature('v2').set('notsol', 'sol1');
model.sol('sol1').feature('v2').set('notsolnum', 'auto');
model.sol('sol1').feature('s2').label('Stationary Solver 2.1');
model.sol('sol1').feature('s2').feature('dDef').label('Direct 3');
model.sol('sol1').feature('s2').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s2').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s2').feature('se1').label('Segregated 1.1');
model.sol('sol1').feature('s2').feature('se1').set('maxsegiter', 300);
model.sol('sol1').feature('s2').feature('se1').set('segstabacc', 'segcflcmp');
model.sol('sol1').feature('s2').feature('se1').set('subinitcfl', 2);
model.sol('sol1').feature('s2').feature('se1').set('subkdpid', 0.03);
model.sol('sol1').feature('s2').feature('se1').set('subcfltol', 0.08);
model.sol('sol1').feature('s2').feature('se1').feature('ss1').label('Velocity u, Pressure p');
model.sol('sol1').feature('s2').feature('se1').feature('ss1').set('segvar', {'comp1_p' 'comp1_u' 'comp1_spf_inl1_Pinlfdf'});
model.sol('sol1').feature('s2').feature('se1').feature('ss1').set('linsolver', 'd1');
model.sol('sol1').feature('s2').feature('se1').feature('ss1').set('subdamp', '0.5');
model.sol('sol1').feature('s2').feature('se1').feature('ss2').label('Turbulence variables');
model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('segvar', {'comp1_k' 'comp1_om'});
model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('linsolver', 'd2');
model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subdamp', '0.35');
model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subtermconst', 'itertol');
model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subiter', 3);
model.sol('sol1').feature('s2').feature('se1').feature('ss2').set('subntolfact', 1);
model.sol('sol1').feature('s2').feature('se1').feature('ll1').label('Lower Limit 1.1');
model.sol('sol1').feature('s2').feature('se1').feature('ll1').set('lowerlimit', 'comp1.k 0 comp1.om 0 ');
model.sol('sol1').feature('s2').feature('d1').label('Direct, fluid flow variables (spf)');
model.sol('sol1').feature('s2').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s2').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s2').feature('d2').label('Direct, turbulence variables (spf)');
model.sol('sol1').feature('s2').feature('d2').set('linsolver', 'pardiso');
model.sol('sol1').feature('s2').feature('d2').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s2').feature('i1').label('AMG, fluid flow variables (spf)');
model.sol('sol1').feature('s2').feature('i1').set('nlinnormuse', true);
model.sol('sol1').feature('s2').feature('i1').set('maxlinit', 1000);
model.sol('sol1').feature('s2').feature('i1').set('rhob', 20);
model.sol('sol1').feature('s2').feature('i1').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('maxcoarsedof', 80000);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('strconn', 0.02);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').label('SCGS 1.1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('iter', 0);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('vankavarsactive', true);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('vankavars', {'comp1_spf_inl1_Pinlfdf'});
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('approxscgs', true);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('pr').feature('sc1').set('scgsdirectmaxsize', 1000);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').label('SCGS 1.1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('iter', 1);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('vankavarsactive', true);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('vankavars', {'comp1_spf_inl1_Pinlfdf'});
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('approxscgs', true);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('po').feature('sc1').set('scgsdirectmaxsize', 1000);
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s2').feature('i1').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').feature('s2').feature('i2').label('AMG, turbulence variables (spf)');
model.sol('sol1').feature('s2').feature('i2').set('nlinnormuse', true);
model.sol('sol1').feature('s2').feature('i2').set('maxlinit', 1000);
model.sol('sol1').feature('s2').feature('i2').set('rhob', 20);
model.sol('sol1').feature('s2').feature('i2').feature('ilDef').label('Incomplete LU 1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').label('Multigrid 1.1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('prefun', 'saamg');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('maxcoarsedof', 50000);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('saamgcompwise', true);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').set('usesmooth', false);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').label('Presmoother 1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('iter', 0);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('linemethod', 'uncoupled');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('pr').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').label('Postsmoother 1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('soDef').label('SOR 1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').label('SOR Line 1.1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('linesweeptype', 'ssor');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('iter', 1);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('linerelax', 0.7);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('linemethod', 'uncoupled');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('po').feature('sl1').set('relax', 0.5);
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').label('Coarse Solver 1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').feature('dDef').label('Direct 2');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').feature('d1').label('Direct 1.1');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').feature('s2').feature('i2').feature('mg1').feature('cs').feature('d1').set('pivotperturb', 1.0E-13);
model.sol('sol1').runAll;

model.save([pwd,'/EC_cfd_solved.mph']);

x0 = mpheval(model,'x','dataonly','on','selection',[2 5 8],'edim','boundary');
y0 = mpheval(model,'y','dataonly','on','selection',[2 5 8],'edim','boundary');
nx0 = mpheval(model,'nx','dataonly','on','selection',[2 5 8],'edim','boundary');
ny0 = mpheval(model,'ny','dataonly','on','selection',[2 5 8],'edim','boundary');
tx0 = mpheval(model,'tx','dataonly','on','selection',[2 5 8],'edim','boundary');
ty0 = mpheval(model,'ty','dataonly','on','selection',[2 5 8],'edim','boundary');
utau0 = mpheval(model,'spf.u_tau','dataonly','on','selection',[2 5 8],'edim','boundary');

r0 = sqrt(x0.^2+y0.^2);
[~,pos] = sort(r0,'ascend');
x0 = x0(pos);
y0 = y0(pos);
nx0 = nx0(pos);
ny0 = ny0(pos);
tx0 = tx0(pos);
ty0 = ty0(pos);

n = numel(x0); % returns the number of elements(n) in array x0
ni = 1000;  % No. of cutlines
n_ = 12;


r1 = zeros(n,ni);
u1 = zeros(n,ni);
v1 = zeros(n,ni);
ut1 = zeros(n,ni);
up1 = zeros(n,ni);
yp1 = zeros(n,ni);

m = zeros(n,ni-1);

delta = zeros(n,1);
xv = zeros(n,1);
yv = zeros(n,1);

x_ = zeros(n,n_);
y_ = zeros(n,n_);
u_ = zeros(n,n_);
v_ = zeros(n,n_);
muT_ = zeros(n,n_);
ux_ = zeros(n,n_);
vx_ = zeros(n,n_);
muTx_ = zeros(n,n_);
uy_ = zeros(n,n_);
vy_ = zeros(n,n_);
muTy_ = zeros(n,n_);

lproj = 1.5e-4;

for i = 1:n
    x1 = linspace(x0(i),x0(i)-lproj*nx0(i),ni);   
    y1 = linspace(y0(i),y0(i)-lproj*ny0(i),ni);
    
    r1(i,:) = sqrt((x1-x0(i)).^2+(y1-y0(i)).^2);  
    
    u1(i,:) = mphinterp(model,'u','coord',[x1;y1]);
    v1(i,:) = mphinterp(model,'v','coord',[x1;y1]);
    
    ut1(i,:) = abs(u1(i,:)*tx0(i) + v1(i,:)*ty0(i));
    up1(i,:) = ut1(i,:)/utau0(i);
    yp1(i,:) = r1(i,:)*utau0(i)/(mu/rho);

    for j = 1:ni-1
        m(i,j) = (up1(i,j+1)-up1(i,j))/(yp1(i,j+1)-yp1(i,j));
    end
    [~,q] = find(abs(m(i,:)-m(i,1))>1e-2);
    q1 = q(1)+1;
    delta(i,1) = r1(i,q1);

    xv(i,1) = x1(q1);
    yv(i,1) = y1(q1);
    
    x_(i,:) = linspace(x0(i),xv(i,1),n_);
    y_(i,:) = linspace(y0(i),yv(i,1),n_);

    u_(i,:) = mphinterp(model,'u','coord',[x_(i,:);y_(i,:)]);
    v_(i,:) = mphinterp(model,'v','coord',[x_(i,:);y_(i,:)]);
    muT_(i,:) = mphinterp(model,'spf.muT','coord',[x_(i,:);y_(i,:)]);

    ux_(i,:) = mphinterp(model,'ux','coord',[x_(i,:);y_(i,:)]);
    vx_(i,:) = mphinterp(model,'vx','coord',[x_(i,:);y_(i,:)]);
    muTx_(i,:) = mphinterp(model,'d(spf.muT,x)','coord',[x_(i,:);y_(i,:)]);

    uy_(i,:) = mphinterp(model,'uy','coord',[x_(i,:);y_(i,:)]);
    vy_(i,:) = mphinterp(model,'vy','coord',[x_(i,:);y_(i,:)]);
    muTy_(i,:) = mphinterp(model,'d(spf.muT,y)','coord',[x_(i,:);y_(i,:)]);
    
end