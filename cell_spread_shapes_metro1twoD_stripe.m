% --------------------------------------------------------------------------------
%                Released under MIT License
% Authors: Hamsini Suresh, Siamak S Shishvan, Andrea Vigliotti 
% Code used to simulate myofibroblasts on adhesive stripes. Results reported in
% Buskermolen et al (2019) "Entropic forces drive cellular contact guidance"
% doi: https://doi.org/10.1101/479071
% --------------------------------------------------------------------------------

% Main function to run MCMC and generate distributions of Gibbs & Helmholtz
% free-energies and cell morphometrics corresponding to the 
% homeostatic temperature beta0
%
function cell_spread_shapes_metro1twoD_stripe(sNote, u0, beta0, coef, Nsteps, swidth, CPucoef, CPvcoef)    

    % ---------------------- Input parameters --------------------------- %
    % sNote  : string to add to output file name
    % u0     : Maximum NURBS control point displacement 
    % beta0  : homeostatic temperature 
    % coef   : double the amplitude of perturbation to control point in MCMC
    % Nsteps : total no. of steps in MCMC
    % swidth : half-width of stripe, in microns
    % CPucoef, CPvcoef : stretch to control points in x- and y-directions
    % respectively, where y-dirn. is parallel to stripe edges    
    % ------------------------------------------------------------------- %
    
    load('INITIALMESH_64.mat','Boun','ldest','Nodes','R0','Tria')
    Rn           = 0.256; % nucleus radius
    r            = swidth*0.018; % convert stripe width in mu-m to simulation units 
    
    sFileName    = sprintf('cell_stripes_%s_beta%.3f.mat', sNote, beta0);
    fprintf(' data saved to %s \n', sFileName);          
               
    %% NURBS paprameters
    p    = 3;
    nSplineNodes = 2;
    knU  = [zeros(1, p), linspace(0, 1, nSplineNodes), ones(1, p)].';
    nCP  = length(knU)-p-1;

    u00 = (linspace(-1, 1, nCP).');
    CPu = u00*ones(1, nCP);
    CPv = ones(nCP, 1)*u00.';               
    ddd = sqrt(CPu.^2+CPv.^2);
    
    %from circle to circle(99%) activate following lines, 
    % otherwise deactivate
    CPu = CPu./ddd;
    CPv = CPv./ddd;
    for i = 1:nCP
        for j = 1:nCP
            if isnan(CPu(i,j))
                CPu(i,j) = 0;
            end
            if isnan(CPv(i,j))
                CPv(i,j) = 0;
            end
        end
    end
    
    CPu = u0*CPucoef*CPu;
    CPv = u0*CPvcoef*CPv;
    CP  = [CPu(:), CPv(:)];    
    
    % Find if an element is within Rn (part of nucleus) or not    
    ctrnuc = 1;
    id_cyto = ones(size(Tria,1),1);    
    for g = 1:1:size(Tria,1)
        idNodes    = Tria(g, :);
        ElNodes    = Nodes(idNodes, :);
        ElCenter   = mean(ElNodes);
        distElCenter = sqrt(ElCenter(1,1)^2 + ElCenter(1,2)^2);
        if distElCenter <= Rn
            id_cyto(g) = 0; % flag = 0 entries form a part of the nucleus
            Nodesnuc1(ctrnuc:ctrnuc+2,:) = ElNodes;
            nucidNodes1(ctrnuc:ctrnuc+2,1) = idNodes';
            ctrnuc = ctrnuc + 3;
        end        
    end        
    Nodesnuc = unique(Nodesnuc1,'rows');        
    nucidNodes = unique(nucidNodes1, 'rows');
    
   
    rrnuc=sqrt((Nodesnuc(:,1)).^2+(Nodesnuc(:,2)).^2);
    k=0;
    sn=size(Nodesnuc,1);
    for i=1:sn
        if(abs(rrnuc(i)-Rn)<0.02) % set tolerance based on boundary points chosen on the nucleus
            k=k+1;
            Bounnuc(k,1)=nucidNodes(i);            
            tet=atan2(Nodesnuc(i,2),Nodesnuc(i,1));
            if (tet < 0)
                tet=tet+2*pi;
            end
            if (abs(tet -2*pi)<1e-10)
                tet=0;
            end
            Bounnuc(k,2)=tet;
            Bounnuc=sortrows(Bounnuc,2);
        end
    end
    Bounnuc=Bounnuc(:,1); % index of boundary nodes of nucleus            
      
    tic;
    rng('shuffle')

    Dispold = zeros(size(Nodes,1),2);
    [rndG, rndGcyto, rndAc, rndAn, rndAsRc, rndAsRn,...
    rndrac, rndran, rndrbc, rndrbn, rndstc, rndstn, ...
    rndphic, rndphin, rndx0c, rndx0n, rndy0c, rndy0n,...
    rndPeric, rndPerin, rndCP, rndtottrac, rndopcsk, acceprate, rnddisp] =  ...
           operations_metro1(sFileName, coef, beta0, Nsteps, CP, knU, p, Nodes, ...
               Tria, Boun, ldest, Dispold, r, id_cyto, Bounnuc);

    toc

    fprintf('done with acceptance rate= %6.2f%% \n', acceprate*100);        
    save (sFileName);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rndG, rndGcyto, rndAc, rndAn, rndAsRc, rndAsRn,...
    rndrac, rndran, rndrbc, rndrbn, rndstc, rndstn, ...
    rndphic, rndphin, rndx0c, rndx0n, rndy0c, rndy0n,...
    rndPeric, rndPerin, rndCP, rndtottrac, rndopcsk, acceprate, rnddisp] = ...
          operations_metro1(sFileName, coef, beta0, Nsteps1, CP, knU, ...
    p, Nodes, Tria, Boun, ldest, Dispold, r, id_cyto, Bounnuc)

     Nsave       = 10000;
     Necho       = 400;
     Burnin       = 1e5;
     idecho      = round(linspace(ceil(Nsteps1/Necho), Nsteps1, Necho));
     Nsteps      = Nsteps1 - Burnin;
    
    nCP = length(knU)-p-1;
    
    rndG     = nan(Nsteps, 1);
    rndGcyto = nan(Nsteps, 1);    
    rndtottrac = nan(Nsteps, 1);
    
    phi6     = (-180:0.25:180)';    
    rndopcsk = zeros(length(phi6),1);
    
    rnddisp  = nan(Nsteps, 128);
    rndAc    = nan(Nsteps, 1);
    rndAn    = nan(Nsteps, 1);
    rndAsRc  = nan(Nsteps, 1);
    rndAsRn  = nan(Nsteps, 1);
    
    rndrac   = nan(Nsteps, 1);
    rndran   = nan(Nsteps, 1);
    rndrbc   = nan(Nsteps, 1);
    rndrbn   = nan(Nsteps, 1);
    rndstc   = nan(Nsteps, 1);
    rndstn   = nan(Nsteps, 1);
    rndphic  = nan(Nsteps, 1);
    rndphin  = nan(Nsteps, 1);
    rndx0c   = nan(Nsteps, 1);
    rndx0n   = nan(Nsteps, 1);
    rndy0c   = nan(Nsteps, 1);
    rndy0n   = nan(Nsteps, 1);
    rndPeric = nan(Nsteps, 1);
    rndPerin = nan(Nsteps, 1);    
    
    Nrjcts  = nan(Nsteps, 1);
    rndCP   = nan(Nsteps, 2*nCP^2);
    Delta   = coef*1;    
    Delta0  = Delta/10;
    
    walking = true;    
    kk      = 1;
    echokk  = 1;
    Nrjct   = 0;        % number of rejected iterations
    
    % For more details on the MCMC calculations, refer to Supplementary S2.3
    % of the above-mentioned paper
    %
    %   build the matrix for the NURBS interpolation
    
    xi  = (Nodes(:,1)-min(Nodes(:,1)))/(max(Nodes(:,1)) - min(Nodes(:,1)));
    eta = (Nodes(:,2)-min(Nodes(:,2)))/(max(Nodes(:,2)) - min(Nodes(:,2)));
    S   = makeNURBSMtx(p, knU, p, knU, [xi, eta]);
    
    % step 1, add some noise to the optimum found previously
    CPkk = CP + (rand(nCP^2, 2)-0.5)*Delta0;  
    Dispkk  = S*CPkk;

    [Gtotsold, Gcytosold, ~, ~, tot_tracold, opcsk_old, Areaold, Areanold, Dispold] ...
                    = get_cell_energy2D_th_twoD_stripe_mex(Dispkk, Nodes, ...
             Tria, Boun, ldest, Dispold, r, id_cyto);         
   
    nodescell = Nodes + Dispold;
    sb = size(Boun,1);
    xc = zeros(sb,1);
    yc = zeros(sb,1);
    for i = 1:sb
       xc(i,1) = nodescell(Boun(i,1),1);
       yc(i,1) = nodescell(Boun(i,1),2);
    end    
    [aoldc, boldc, AsRoldc, phioldc, statusoldc, X0oldc, Y0oldc] = bstelpsmetro(xc,yc);
    perioldc = get_calc_perimeter(xc, yc, sb);
    
    sn = size(Bounnuc,1);
    xn = zeros(sn,1);
    yn = zeros(sn,1); 
    for i = 1:sn
       xn(i,1) = nodescell(Bounnuc(i,1),1);
       yn(i,1) = nodescell(Bounnuc(i,1),2);
    end
    [aoldn, boldn, AsRoldn, phioldn, statusoldn, X0oldn, Y0oldn] = bstelpsmetro(xn,yn);   
    perioldn = get_calc_perimeter(xn, yn, sn);

    while walking
        % step 2, select which control point to modify
        idCP       = round(rand() * (nCP^2-1))+1;
        %step 3, make a trial step
        CPtrial    = CPkk;
        rndstep    = (rand(1,2)-0.5)*Delta;
        CPtrial(idCP,:) =  CPtrial(idCP,:) + rndstep;
        Disptrial  = S*CPtrial;
        % step 4, evalute energy and compute DeltaG
        [Gtotsnew, Gcytosnew, ~, ~, tot_tracnew, opcsk_new, Areanew, Areannew, Dispnew] ...
                    = get_cell_energy2D_th_twoD_stripe_mex(Disptrial, Nodes, ...
             Tria, Boun, ldest, Dispold, r, id_cyto);      
        
        DeltaG = Gtotsnew-Gtotsold;
        
        nodescell = Nodes + Dispnew;
        sb = size(Boun,1);
        xc = zeros(sb,1);
        yc = zeros(sb,1);        
        for i = 1:sb
           xc(i,1) = nodescell(Boun(i,1),1);
           yc(i,1) = nodescell(Boun(i,1),2);
        end
        Dispbordernew = [xc; yc];
        
        sn = size(Bounnuc,1);
        xn = zeros(sn,1);
        yn = zeros(sn,1); 
        for i = 1:sn
           xn(i,1) = nodescell(Bounnuc(i,1),1);
           yn(i,1) = nodescell(Bounnuc(i,1),2);
        end
        
        % step 5 and 6, decide whether accept or reject the step
        % probability of acceptance, 0<pacc<1 uniformly distributed
        pacc = rand();
        % if DeltaG<0 we have passed, otherwise the inequality holds
        cnd1 = exp(-DeltaG*beta0) > pacc;
        % conditions for rejecting iterations
        cnd2 = abs(DeltaG*beta0) < 10;
        
        if cnd1 && cnd2
            CPkk        = CPtrial;
            Gtotsold    = Gtotsnew;
            Gcytosold   = Gcytosnew;      
            tot_tracold = tot_tracnew;
            Areaold     = Areanew;
            Areanold    = Areannew;
            Dispold     = Dispnew; 
            Dispborderold = Dispbordernew;
            opcsk_old   = opcsk_new;
                        
            if kk > Burnin
                kkk = kk - Burnin;
                rndCP(kkk,:) = CPkk(:);
                rnddisp(kkk,:) = Dispborderold;
            end
            
            [anewc, bnewc, AsRnewc, phinewc, statusnewc, X0newc, Y0newc] = bstelpsmetro(xc,yc);
            [anewn, bnewn, AsRnewn, phinewn, statusnewn, X0newn, Y0newn] = bstelpsmetro(xn,yn);
            
            aoldc        = anewc;
            aoldn        = anewn;
            
            boldc        = bnewc;
            boldn        = bnewn;
            
            AsRoldc      = AsRnewc;
            AsRoldn      = AsRnewn;
            
            phioldc      = phinewc;
            phioldn      = phinewn;
            
            statusoldc   = statusnewc;
            statusoldn   = statusnewn;
            
            X0oldc       = X0newc;
            X0oldn       = X0newn;
            
            Y0oldc       = Y0newc;
            Y0oldn       = Y0newn;
            
            perinewc     = get_calc_perimeter(xc, yc, sb);
            perinewn     = get_calc_perimeter(xn, yn, sn);
            
            perioldc     = perinewc;
            perioldn     = perinewn;
        else
            Nrjct       = Nrjct+1;
        end
        
        if kk > Burnin
            kkk = kk - Burnin;
            rndG(kkk)    = Gtotsold;
            rndGcyto(kkk)    = Gcytosold;            
            rndtottrac(kkk) = tot_tracold;
            rndopcsk      = rndopcsk + opcsk_old;

            rndAc(kkk)    = Areaold;
            rndAn(kkk)    = Areanold;

            rndrac(kkk)   = aoldc;
            rndran(kkk)   = aoldn;

            rndrbc(kkk)   = boldc;
            rndrbn(kkk)   = boldn;

            rndAsRc(kkk)  = AsRoldc;
            rndAsRn(kkk)  = AsRoldn;

            rndphic(kkk)  = phioldc;
            rndphin(kkk)  = phioldn;

            rndstc(kkk)   = statusoldc;
            rndstn(kkk)   = statusoldn;

            rndx0c(kkk)   = X0oldc;
            rndx0n(kk)   = X0oldn;

            rndy0c(kkk)   = Y0oldc;
            rndy0n(kkk)   = Y0oldn;

            rndPeric(kkk) = perioldc;
            rndPerin(kkk) = perioldn;                        
            
            Nrjcts(kkk)  = Nrjct;  
        end
                        
        
        if kk >= idecho(echokk) && echokk < Necho && ~mod(kk, 2e3) && kk > Burnin
           kkk = kk - Burnin; 
           meanrndG = mean(rndG(1:kkk));
           meanrndAc = mean(rndAc(1:kkk));
           meanrndAsRc = mean(rndAsRc(1:kkk));
           meanrndAn = mean(rndAn(1:kkk));
           meanrndAsRn = mean(rndAsRn(1:kkk));
        
           fprintf('kk:%8i, done:%6.2f%%, rejected:%6.2f%%, meanG:%6.2f, meanAc:%6.2f, meanAsRc:%6.2f, meanAn:%6.2f, meanAsRn:%6.2f \n', ...
                    kk, 100*kk/Nsteps, 100*Nrjct/kk, meanrndG, meanrndAc, meanrndAsRc, meanrndAn, meanrndAsRn)
            echokk = echokk+1;
        end
        
        acceprate=1-(Nrjct/kk);
        
        if mod(kk,Nsave)==0
            save (sFileName)
        end
        
        kk          = kk+1;
        % check if we are done
        if kk >= Nsteps1+1
            walking = false;
        end
        
    end

end

%%   Cell morphometrics                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [long_axis, short_axis, boacurrent, orient, status, X0_in, Y0_in] = bstelpsmetro(x,y)

status = 1;
% remove bias of the ellipse - to make matrix inversion more accurate. (will be added later on).
ave_size=((max(x)-min(x))+(max(y)-min(y)))/2;
mean_x = mean(x);
mean_y = mean(y);
cent_x=mean_x;
cent_y=mean_y;
x = x-mean_x;
y = y-mean_y;

% the estimation for the conic equation of the ellipse
X = [x.^2, x.*y, y.^2, x, y ];
a = sum(X)/(X'*X);

% extract parameters from the conic equation
[a,b,c,d,e] = deal( a(1),a(2),a(3),a(4),a(5) );

% remove the orientation from the ellipse
    orientation_rad = 1/2 * atan2( b, (c-a) );
%%%%%    orientation_rad = 1/2 * atan( b/(c-a) );
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    [a,b,c,d,e] = deal(...
        a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
        0,...
        a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
        d*cos_phi - e*sin_phi,...
        d*sin_phi + e*cos_phi );
   [mean_x,mean_y] = deal( ...
       cos_phi*mean_x - sin_phi*mean_y,...
       sin_phi*mean_x + cos_phi*mean_y );

% check if conic equation represents an ellipse
test = a*c;

% if we found an ellipse return it's data
if (test>0)

    % make sure coefficients are positive as required
    if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end

    % final ellipse parameters
    X0          = mean_x - d/2/a;
    Y0          = mean_y - e/2/c;
    F           = 1 + (d^2)/(4*a) + (e^2)/(4*c);
    [a,b]       = deal( sqrt( F/a ),sqrt( F/c ) );
    long_axis   = 2*max(a,b);
    short_axis  = 2*min(a,b);
    boacurrent=long_axis/short_axis;
    
    orientation_rad = -orientation_rad;
    if orientation_rad < 0
        orient = orientation_rad +pi;
    else
        orient = orientation_rad;
    end
    % rotate the axes backwards to find the center point of the original TILTED ellipse
    R           = [ cos_phi sin_phi; -sin_phi cos_phi ];
    P_in        = R * [X0;Y0];
    X0_in       = P_in(1);
    Y0_in       = P_in(2);

    diff=sqrt( (X0_in-cent_x)^2+(Y0_in-cent_y)^2 );
    if diff > ave_size/2
        [long_axis, short_axis, boacurrent, orient, X0_in, Y0_in] = bstelpsmetro_alternative(x,y);
        status = 2;
%         disp( 'wrong fit!! use alternative function' );
    end

elseif (test==0)
    [long_axis, short_axis, boacurrent, orient, X0_in, Y0_in] = bstelpsmetro_alternative(x,y);
    status = 3;
%     disp( 'bstelps=0: Did not locate an ellipse; use alternative function' );
else
    [long_axis, short_axis, boacurrent, orient, X0_in, Y0_in] = bstelpsmetro_alternative(x,y);
    status = 4;
%     disp( 'bstelps<0: Did not locate an ellipse; use alternative function' );
end

end
function [long_axis, short_axis, boacurrent, orient, uCentre, vCentre] = bstelpsmetro_alternative(X,Y)  

    % normalize data
    mx = mean(X);
    my = mean(Y);
    sx = (max(X)-min(X))/2;
    sy = (max(Y)-min(Y))/2;

    x = (X-mx)/sx;
    y = (Y-my)/sy;

    % Force to column vectors
    x = x(:);
    y = y(:);

    % Build design matrix
    D = [ x.^2  x.*y  y.^2  x  y  ones(size(x)) ];

    % Build scatter matrix
    S = D'*D;

    % Build 6x6 constraint matrix
    % C(6,6) = 0;
    C = zeros(6);
    C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

    % Solve eigensystem
    % New way, numerically stabler in C [gevec, geval] = eig(S,C);

    % Break into blocks
    tmpA = S(1:3,1:3);
    tmpB = S(1:3,4:6);
    tmpC = S(4:6,4:6);
    tmpD = C(1:3,1:3);
    % tmpE = inv(tmpC)*tmpB';
    tmpE = tmpC\tmpB.';
    % [evec_x, eval_x] = eig(inv(tmpD) * (tmpA - tmpB*tmpE));
    [evec_x, eval_x] = eig(tmpD\(tmpA - tmpB*tmpE));

    % Find the positive (as det(tmpD) < 0) eigenvalue
    I = real(diag(eval_x)) < 1e-8 & ~isinf(diag(eval_x));

    % Extract eigenvector corresponding to negative eigenvalue
    A = real(evec_x(:,I));

    % Recover the bottom half...
    evec_y = -tmpE * A;
    A = [A; evec_y];

    % unnormalize
    par = [
        A(1)*sy*sy,   ...
        A(2)*sx*sy,   ...
        A(3)*sx*sx,   ...
        -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
        -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
        A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
        - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
        + A(6)*sx*sx*sy*sy   ...
        ]';

    % Convert to geometric radii, and centers

    thetarad = 0.5*atan2(par(2),par(1) - par(3));
    cost = cos(thetarad);
    sint = sin(thetarad);
    sin_squared = sint.*sint;
    cos_squared = cost.*cost;
    cos_sin = sint .* cost;

    Ao = par(6);
    Au =   par(4) .* cost + par(5) .* sint;
    Av = - par(4) .* sint + par(5) .* cost;
    Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
    Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;

    % ROTATED = [Ao Au Av Auu Avv]

    tuCentre = - Au./(2.*Auu);
    tvCentre = - Av./(2.*Avv);
    wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;

    uCentre = tuCentre .* cost - tvCentre .* sint;
    vCentre = tuCentre .* sint + tvCentre .* cost;

    Ru = -wCentre./Auu;
    Rv = -wCentre./Avv;

    Ru = sqrt(abs(Ru)).*sign(Ru);
    Rv = sqrt(abs(Rv)).*sign(Rv);

   if (Ru > Rv)
        long_axis = Ru;
        short_axis = Rv;
        boacurrent = Ru / Rv;
        orientation_rad = thetarad;
   else
        long_axis = Rv;
        short_axis = Ru;
        boacurrent = Rv / Ru;
        orientation_rad = thetarad-pi/2;
   end

    if orientation_rad < 0
        orient = orientation_rad +pi;
    else
        orient = orientation_rad;
    end

end
function [ldist] = get_calc_perimeter(x, y, segi)

    ldist = 0;    
    for k=1:1:segi
        start = k;

        if k==segi
            stop = 1;        
        else
            stop = k + 1;
        end

        ldist = ldist + sqrt( (x(stop) - x(start))^2 + (y(stop) - y(start))^2 );
    end      
    
end
%%  NURBS functions                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S              = makeNURBSMtx(p, knxi, q, kneta, xieta)

    n     = length(knxi) -p-1;        % number of functions in xi
    m     = length(kneta)-q-1;        % number of functions in eta
    nPts  = length(xieta);

    S     = zeros(nPts, m*n);

    for ii=1:nPts

        [Nxi, idxi]     = NURBS1D(n-1, p, knxi,  xieta(ii, 1));
        [Neta, ideta]   = NURBS1D(m-1, q, kneta, xieta(ii, 2));
        %
        % this is for the matrix type arrangement of the control points
        %
        S0      = Nxi.'*Neta;
        nidxi   =  numel(idxi);
        nideta  =  numel(ideta);
        idx0    = idxi.'*ones(1, nideta) + ...
            ones(nidxi, 1)*((ideta-1)*n);
        
        S(ii, idx0(:)) = S0(:);
    end
end
function [Nxi, idxi] = NURBS1D(n, p, knxi, xi)
    
    uspan   = FindSpan(n,p,xi,knxi);
    Nxi     = Der1BasisFun(uspan, xi,p,knxi);
    idxi    = (uspan-p+1) + (0:p);
end
function knotSpanIndex = FindSpan(n,p,u,U)
    %--------------------------------------------------------------
    %function knotSpanIndex = FindSpan(n,p,u,U)
    % NURBS-Book (algorithm A2.1)
    % find the knot span index for one variable u
    %INPUT:
    % n          : number of basis function -1
    %        NURBS-Book: np+1 # basis, np max index (startindex 0)
    %        here        np   # basis and max index (startindex 1)
    % p          : degree of the basis functions
    % u          : evaluation point
    % U          : knot vector (row vector)
    %OUTPUT:
    % knotSpanIndex : index of knot span
    %--------------------------------------------------------------
    if (u == U(n+2))
        knotSpanIndex= n;
        return
    end
    low = p;
    high = n+1;
    mid = floor((low + high)/2);
    while (u <U(mid+1) || u >= U(mid+2) )
        if( u < U(mid+1))
            high = mid;
        else
            low = mid;
        end
        mid = floor((low+high)/2);
    end
    knotSpanIndex = mid;
end
function ders = Der1BasisFun(i,u,p,U)
    %--------------------------------------------------------------
    %function ders = Der1BasisFun(i,u,p,U)
    % NURBS-Book modified (algorithm A2.3)
    % evalute nonzero basis functions and first derivative
    %INPUT:
    % i          : current knotspan
    % u          : evaluation point
    % p          : degree of the basis functions
    % U          : knot vector
    %OUTPUT:
    % ders       : matrix (2 , p+1)
    %              1. row vector (dim p+1)
    %              values of the basis function N_(i-p) ... N_(i)
    %              at the evaluation point u
    %              2. first derivation
    %--------------------------------------------------------------
    ders = zeros(1,p+1);
    N=zeros(p+1,p+1);
    N(1,1)=1;
    left=zeros(1,p+1);
    right=zeros(1,p+1);
    for j=1:p
        left(j+1) = u-U(i+1-j+1);
        right(j+1) = U(i+j+1)-u;
        saved = 0;
        for r=0:j-1
            N(j+1,r+1) = right(r+2) + left(j-r+1);
            temp = N(r+1,j)/N(j+1,r+1);
            N(r+1,j+1) = saved + right(r+2)*temp;
            saved = left(j-r+1)*temp;
        end
        N(j+1,j+1) = saved;
    end
    for j=0:p
        ders(1,j+1) = N(j+1,p+1);
    end
end
