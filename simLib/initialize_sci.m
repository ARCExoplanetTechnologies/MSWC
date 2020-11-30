function sci = initialize_sci(params)
%lambdaRef, flD, samplesPerflD, FOVflD)
lambdaRef = params.lambdaRef;
flD = params.flD;
samplesPerflD = params.samplesPerflD;
FOVflD = params.FOVflD;

% initialize science plane
sci.FOVflD = FOVflD;
sci.flD = flD;
sci.samplesPerflD = samplesPerflD;
sci.pixAngSize = 1/samplesPerflD;
sci.dx = flD / samplesPerflD;
sci.dy = flD / samplesPerflD;
sci.lambda_ref = lambdaRef; % reference lambda for pixel scaling
sci.N = ceil(FOVflD * samplesPerflD);
sci.gridsize = sci.N*sci.dx;
sci.x = -(sci.gridsize - sci.dx)/2 : sci.dx : (sci.gridsize - sci.dx)/2;
sci.xlD = sci.x / sci.flD;
sci.y = -(sci.gridsize - sci.dy)/2 : sci.dy : (sci.gridsize - sci.dy)/2;
sci.ylD = sci.y / sci.flD;
[sci.xx sci.yy] = meshgrid(sci.x, sci.y);
sci.xxlD = sci.xx / sci.flD;
sci.yylD = sci.yy / sci.flD;
sci.rr = sqrt(sci.xx.^2 + sci.yy.^2);
sci.rrlD = sci.rr / sci.flD;

