% To run all scripts and create all figures, tables from Chapter 3

% Displacement formulations
clear all

NX = 80;
NY = 80;

p = 1e-01; pressure='Pa';
run('Displ_R.m');
run('benchmarkD.m');

p = 1e+06; pressure= 'MPa';
run('Displ_R.m');
run('benchmarkD.m');

NX = 30;
NY = 30;

p = 1e-01; pressure='Pa';
run('Displ_R.m');
run('benchmarkD.m');
run('DisplDirichlet_R.m');
run('benchmarkD.m');
run('DisplDirichlet_M.m');
run('benchmarkD.m');

p = 1e+06; pressure= 'MPa';
run('Displ_R.m');
run('benchmarkD.m');
run('DisplDirichlet_R.m');
run('benchmarkD.m');
run('DisplDirichlet_M.m');
run('benchmarkD.m');


% Stress formulations
clear all

NX = 40;
NY = 40;
b = 100;
c = 100;
DOMX = [-c,c];
DOMY = [0,b];
p = 1e-01; pressure='Pa';

letter='A';
run('StressDirichlet_R.m');
run('benchmarkS.m');
run('StressDirichlet_M.m');
run('benchmarkS.m');

letter='B';
run('StressDirichlet_R.m');
run('benchmarkS.m');
run('StressDirichlet_M.m');
run('benchmarkS.m');

NX = 100;
NY = 100;

letter='A';
run('StressDirichlet_R.m');
% Reduce oscillations - smoothdata
MATsigx=smoothdata(MATsigx,1);
MATsigy=smoothdata(MATsigy,1);
MATsh=smoothdata(MATsh,1);
run('benchmarkS.m');
run('StressDirichlet_M.m');
run('benchmarkS.m');

letter='B';
run('StressDirichlet_R.m');
% Reduce oscillations - smoothdata
MATsigx=smoothdata(MATsigx,2);
MATsigy=smoothdata(MATsigy,2);
MATsh=smoothdata(MATsh,2);
run('benchmarkS.m');
run('StressDirichlet_M.m');
run('benchmarkS.m');

NX = 30;
NY = 30;
p = 1e-01; pressure='Pa';
run('StressDirichlet_LS.m');
run('benchmarkS.m');

clear all
% Rectagle sufficiently distant from x=+-a
NX = 40;
NY = 40;
a = 10; % length under the pressure p
b = 100;
c = 100;
DOMX = [a+10,c];
DOMY = [0,b];
letter='A';
p = 1e-01; pressure='Pa';
run('StressDirichlet_R.m');
ScriptName=join(['StressDirichlet_R','_',pressure,'_',string(NX),string(NY),letter,'Domain'],'');
run('benchmarkS.m');