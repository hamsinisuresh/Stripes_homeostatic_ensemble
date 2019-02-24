% --------------------------------------------------------------------------------
%                Released under MIT License
% Authors: Hamsini Suresh, Siamak S Shishvan, Andrea Vigliotti 
% Numerical simulation of myofibroblasts on adhesive stripes. Results reported in
% Buskermolen et al (2019) "Entropic forces drive cellular contact guidance"
% doi: https://doi.org/10.1101/479071
% --------------------------------------------------------------------------------

% Code to calculate strains, stresses within the cell and nucleus,
% Helmholtz and Gibbs free-energies and cellular tractions under
% specified displacement.
% 
%{
    % to produce the mex version run the following commands at prompt:

    double1Da = coder.typeof(double(1), [inf, 1], [false, false]);
    double2Da = coder.typeof(double(1), [inf, 2], [false, false]);
    double2Db = coder.typeof(double(1), [inf, 3], [false, false]);
    double2Dc = coder.typeof(double(1), [inf, 9], [false, false]);
    sArgs = {double2Da, ...
        double2Da, double2Db, double1Da, ...
        double2Dc, double2Da, double(1), double1Da};
    codegen get_cell_energy2D_th_twoD_stripe -args sArgs

%}
function [Gtotstar, Gcytostar, ...
        Gelasstar, GelasEl, tot_trac, opcsk_dist, ...
        Area, Arean, Disp, epsn, ...
        lmdn, Cel, Fel] = get_cell_energy2D_th_twoD_stripe(Disptrial, ...
        Nodes, Tria, Boun, ...
        ldest, Dispold, r, id_cyto)
    
    coder.extrinsic('fwrite');        
    %% model parameters
    
    kt          = 4.279e-6;     % boltzmann constant nN*mum
    muT         = 20*kt;        % transition potential
    muU0        = 8*kt;         % unbound potential
    muB0        = 9*kt;         % bound potential
    r_Deltamu   = 0;            % \Delta\mu/(muT-muB0)
    sigmax      = 240;          % max fibre tension kPa
    eps_s       = 0.3;          % isometric tension plateau parameter
    eps_p       = 0.6;          % isometric tension slope parameter
    edot0       = 0.8/6;        % hill relation
    Omega       = 10^(-7.101);  %
    beta_       = 1.2;
    Cmax        = 1.0;
    n0          = 3e6;          % total number of proteins within the RVE    
    f0          = 0.1/pi;
    b0          = 0.2*(1/sqrt(pi)); % initial thickness =R0/5; ALSO V0=pi*R0^2*b0=b0
    K2          = 10^6;         % co-efficient for elastic penalty term
    Jc          = 0.6;          % threshold used in elastic penalty term
    
    %% HVSC parameters
    
    Ecell       = 5;     % Young's modulus of cell, in kPa
    Enuc        = 10;    % Young's modulus of nucleaus, in kPa   
    alpha       = 5;     % Material constant for cell in Ogden-type 2D model 
    Kappa       = 35;    % In-plane bulk modulus for cell
    etamax      = 0.75;  % Limit on SF formation at each angle in an RVE
    alpha_nuc   = 20;    % Material constant for nucleus in Ogden-type 2D model
    Kappa_nuc   = 35;    % In-plane bluk modulus for nucleus
        
    
    %% phi and RR
    
    nphi        = 24;           % number of directions    
    phi         = ((0:nphi-1)/nphi-1/2)*pi;
    phi         = phi(:);
    cosphi2     = cos(phi).^2;
    sinphi2     = sin(phi).^2;
    sin2phi     = sin(2*phi);
    RR          = [cosphi2,  sinphi2, sin2phi];    
    phi6        = -180:0.25:180;   % finer discretisation to get smooth CSK order parameter distribution 
    opcsk_dist  = zeros(length(phi6),1);
          
    Disp = interx(Nodes, Boun, Dispold, Disptrial, r);
    
    % "nodescell" is the current shape of the cell    
    nodescell = Nodes + Disp;
       
    %% loop over the elements    
    nTria    = size(Tria, 1);
    GelasEl  = zeros(nTria,1);
    Area0El  = zeros(nTria,1);            
    Nb_Nu_pi = nan(nTria, 1);
    n_f      = nan(nphi, nTria);    
    eterm    = nan(nphi, nTria);
    etahat   = nan(nphi, nTria);
    netterm  = nan(nphi, nTria);
    
    lmdn     = nan(nphi, nTria);
    epsn     = nan(3, nTria);
    Cel      = nan(3, nTria);
    Fel      = nan(4, nTria);
    thet     = nan(1,nTria);
    detFEl   = nan(1,nTria);
    thickness= zeros(1,nTria);
    sigcpel  = zeros(3,nTria);
    lmd12    = zeros(2,nTria);
    nu1      = nan(2,nTria);
    nu2      = nan(2,nTria);
      
    Area0 = 0;
    Area  = 0;    
    Areac = 0;
            
    for kk = 1:nTria
        idNodes    = Tria(kk, :);
        ElNodes    = Nodes(idNodes, :);
        ElDisp     = Disp(idNodes, :);
        
        % Calculate principal stretches in each RVE and areas of cell and nucleus
        [Area0El(kk), Cel(:,kk), ...
            epsn(:,kk), lmd12(:,kk), nu1(:,kk), nu2(:,kk), detFEl(kk), Fel(:,kk), thet(kk)] ...
            = calcA0andC(ElNodes, ElDisp);  
        
        Area0      = Area0 + Area0El(kk);
        Area       = Area  + Area0El(kk)*detFEl(kk);
        Areac      = Areac + id_cyto(kk)*Area0El(kk)*detFEl(kk);
        
        lmdn(:,kk) = RR*(epsn(:,kk)+[1;1;0]);                
                   
        % Calculate the distribution of bound SF units in the cell
        % 
        if(id_cyto(kk))             % part of the cytosol        
            [Nb_Nu_pi(kk), n_f(:,kk), eterm(:,kk)] = get_Nb_Nu_pi(lmdn(:,kk), beta_, ...
            edot0, sigmax, eps_s, eps_p, ...
            kt, muT, muU0, muB0, ...
            r_Deltamu, Omega, Cmax);      
        else                        % part of the nucleus
            Nb_Nu_pi(kk) = 0;
            n_f(:,kk) = nan(nphi,1);
            eterm(:,kk) = nan(nphi,1);
        end                
    end
    Arean = Area - Areac;
    
 flipping = min(detFEl);
 
 % If elements in mesh have NOT flipped
 if flipping > (Jc-1)

    % Ogden-type 2D model
    for kk = 1:nTria
        G0 = Ecell/3;
        Gnuc = Enuc/3;
        lamthick = 1/(lmd12(1,kk)*lmd12(2,kk));
        thickness(kk) = b0*lamthick;
        lmdb1 = sqrt(lmd12(1,kk)/lmd12(2,kk));
        lmdb2 = 1/lmdb1;
        J = lmd12(1,kk)*lmd12(2,kk);
        hj = Jc>J;                        
        
             if(id_cyto(kk)==0)         % part of nucleus
                 simple = lmdb1^alpha_nuc + lmdb2^alpha_nuc - 2;
                 GelasEl(kk) = 2*Gnuc/alpha_nuc^2*simple + Kappa_nuc/2*(abs(J-1))^2 - K2*log(J+1-Jc)*hj;
                 term1 = Gnuc/alpha_nuc*(lmdb1^alpha_nuc - lmdb2^alpha_nuc);
                 term2 = Kappa_nuc*J*(J-1) - K2*(J/(J+1-Jc))*hj;
             else                       % part of the cytosol
                 simple = lmdb1^alpha + lmdb2^alpha - 2;
                 GelasEl(kk) = 2*G0/alpha^2*simple + Kappa/2*(abs(J-1))^2 - K2*log(J+1-Jc)*hj;
                 term1 = G0/alpha*(lmdb1^alpha - lmdb2^alpha);
                 term2 = Kappa*J*(J-1) - K2*(J/(J+1-Jc))*hj;
             end
                        
        sig1 = term1 + term2;
        sig2 =-term1 + term2;
        nnu1 = nu1(:,kk);
        nnu2 = nu2(:,kk);
        sigmacart = sig1*(nnu1*nnu1') + sig2*(nnu2*nnu2');
        
        % stresses within the cell from passive hyperelastic components
        sigcpel(1,kk) = sigmacart(1,1);
        sigcpel(2,kk) = sigmacart(2,2);
        sigcpel(3,kk) = sigmacart(1,2);
    end            
       
    Nuhat     = 1/( 1 + sum(Nb_Nu_pi.*Area0El)/sum(id_cyto.*Area0El) );
        
    delphi = pi/nphi;    
    % Calculate concentration of unbound proteins in all RVEs in cell
    Nuhat_Actual = bisection(Nuhat, 1, etamax, eterm, n_f, Area0El, nTria, delphi, id_cyto); % lower bound = old Nuhat, upper bound = 1

    %% done
    % cytoskeleton free energy
    Gcytostar = (n0*sum(id_cyto.*Area0El)*kt*log(Nuhat_Actual))/Area0;
    
    % elastic free-energy
    Gelasstar = sum(GelasEl.*Area0El)/Area0;    
    

     for kk = 1:nTria
        if(id_cyto(kk))
            phistar     = phi+thet(kk);
            cosphistar2 = cos(phistar).^2;
            sinphistar2 = sin(phistar).^2;
            sin2phistar = sin(2*phistar);

            etahat(:,kk) = (Nuhat_Actual*eterm(:,kk)/pi)./(n_f(:,kk)+Nuhat_Actual*eterm(:,kk)./(pi*etamax));
            netterm(:,kk) = etahat(:,kk).*n_f(:,kk);           
          % Nthat_Actual(kk) = Nuhat_Actual + pi*(sum(netterm(:,kk)))/nphi;

            etahatss = (Nuhat_Actual*eterm(:,kk)/pi)./(n_f(:,kk)+Nuhat_Actual*eterm(:,kk)./(pi*etamax));
            % note that sigcpel:(sig11,sig22,sig12)
            here1= f0*sigmax/detFEl(kk);
            here=etahatss.*lmdn(:,kk);
            
            % stresses in cell from active SF cytoskeleton
            sigcpel(1,kk) = sigcpel(1,kk) + here1 * mean(here.*cosphistar2)*pi;
            sigcpel(2,kk) = sigcpel(2,kk) + here1 * mean(here.*sinphistar2)*pi;
            sigcpel(3,kk) = sigcpel(3,kk) + here1 * mean(here.*sin2phistar)*pi*0.5;
        end
     end
    
     % Calculate distribution of rotationally-invariant SF orientations in 
     % a given configuration
     nettermx = nan(nphi, nTria);
     nettermy = nan(nphi, nTria);
     for nel = 1:nTria
            if id_cyto(nel)         
                  nettermx(:,nel) = netterm(:,nel).*cos(phi+thet(nel));
                  nettermy(:,nel) = netterm(:,nel).*sin(phi+thet(nel));
            end        
     end


    maxphi2 = nan(nTria,1);
    
    for nel = 1:nTria
       total_dirn = nan(nphi,2);       

       for gg = 1:1:nphi
        if ~isnan(nettermx(gg, nel)) && ~isnan(nettermy(gg, nel))
            mag = norm(nettermx(gg, nel)+1i*nettermy(gg,nel));
            dirn = atand(nettermy(gg,nel)/nettermx(gg,nel));

            if dirn < 0
                dirn = dirn + 180;
            end
            total_dirn(gg,1) = mag*Area0El(nel)/Area0;
            total_dirn(gg,2) = dirn;                 
        end
       end   

       % Dirn. of maximum concentration of SF in RVE
       if id_cyto(nel)   
           [~,maxi2] = max(total_dirn(:,1));           
           maxphi2(nel) = total_dirn(maxi2,2);           
       end
    end

    maxphi2csk = maxphi2(~isnan(maxphi2));   

    % average csk orientation
    csk_avg = mean(maxphi2csk);
    maxphi2csk = maxphi2csk - csk_avg;       
    [maxphi21csk, ~] = sort(maxphi2csk);  

    % distribution of SF orientations in given cell configuration
    for pp = 1:length(maxphi21csk)
        [~, ll] = min(abs(maxphi21csk(pp) - phi6));    
        opcsk_dist(ll) = opcsk_dist(ll) + 1;        
    end   

    % calculate cellular tractions (in kPa)
    [~,~,~,~,tot_trac] = calc_traction(sigcpel, thickness, nodescell, Tria, ldest);   

    % calculate total Gibbs free-energy of the system
    Gtotstar  = Gcytostar + Gelasstar;  
    
 else
    Gtotstar  = inf;
    Gcytostar = inf;
    Gelasstar = inf;
    GelasEl   = inf;    
    tot_trac = inf;
    maxphi = nan;
 end
%    fprintf('Gtotstar : %.8f, Gcytostar : %.8f, Gelasstar : %.8f, Area : %.4f\n', Gtotstar, Gcytostar, Gelasstar, Area);
    
end
%
function [A, Nx, Ny] = calcNvectors(p0)
    %
    x1 = p0(1,1); x2 = p0(2,1); x3 = p0(3,1);
    y1 = p0(1,2); y2 = p0(2,2); y3 = p0(3,2);
    %
    Delta = x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2;
    Nx    = [y2 - y3; y3 - y1; y1 - y2]/Delta;
    Ny    = [x3 - x2; x1 - x3; x2 - x1]/Delta; 
    A     = abs(Delta)/2;
end
%
function [A0, C, epsn, lmd, nu1, nu2, detF, F, thet] = calcA0andC(p0, Disp)

    idxu = (1:2:6); idxv = (2:2:6);
    [A0, Nx, Ny] = calcNvectors(p0);
    
    Nux = zeros(1,6); Nuy = zeros(1,6);
    Nvx = zeros(1,6); Nvy = zeros(1,6);
    
    Nux(idxu) = Nx; Nuy(idxu) = Ny;
    Nvx(idxv) = Nx; Nvy(idxv) = Ny;
    
    % G = [dudX, dudY; dvdX, dvdY];
    u0       = zeros(6,1);
    u0(idxu) = Disp(:,1);
    u0(idxv) = Disp(:,2);
    F        = zeros(2);
    F(:)     = [1; 0; 0; 1] + [Nux; Nvx; Nuy; Nvy]*u0;    
    C        = F.'*F;
 
    if det(F) > 1e-14
        
        [u, lmdn2] = eig2x2sym(C);
        
        if min(lmdn2) > 1e-15
            lmd      = sqrt(lmdn2);
            mtstrch  = u*diag(lmd)*u.';
            rmat=F/mtstrch;

            if abs(rmat(2,1))>1
                disp(rmat)
                rmat(2,1)=sign(rmat(2,1));
            end

            thet     = asin(rmat(2,1));

            nu1=rmat*u(:,1);
            nu2=rmat*u(:,2);

            epsn     = u*diag(lmd-1)*u.';
            epsn     = epsn([1; 4; 2]);
            C        = C([1; 4; 2]);
            detF     = F(1)*F(4) - F(2)*F(3);
            F        = F(:);
        else
%         disp('lmdn2 is negative'); %this configuration has to be rejected.
            C = C([1; 4; 2]); epsn = [1; 1; 1]; lmd =[1; 1]; 
            nu1 = [1; 1]; nu2 = [1; 1]; detF = -1; F = F(:); thet =1;
        end
    else
%         disp('det(F) is zero or negative'); %this configuration has to be rejected.
        C = C([1; 4; 2]); epsn = [1; 1; 1]; lmd =[1; 1]; 
        nu1 = [1; 1]; nu2 = [1; 1]; detF = -1; F = F(:); thet =1;
    end
end
%%
function [u, lmd] = eig2x2sym(C)

m11s=C(1,1);
m21s=C(2,1);
m22s=C(2,2);

order=true;

tau = 0.5*(m22s - m11s)./m21s;
sign_tau=sign(tau);
if tau==0
    sign_tau=1;
end
t = sign_tau./(abs(tau) + sqrt(1 + tau.^2));
c = 1./sqrt(1 + t.^2);
s = c.*t;
l1s = m11s - t.*m21s;
l2s = m22s + t.*m21s;
e11s = c;
e21s = -s;

swaps = abs(l2s) > abs(l1s);
e11s(swaps) = s(swaps);
e21s(swaps) = c(swaps);
temp = l1s(swaps);
l1s(swaps) = l2s(swaps);
l2s(swaps) = temp;

divzeros = m21s == 0;
order = m11s(divzeros) >= m22s(divzeros);
e11s(divzeros) = order;
e21s(divzeros) = 1 - order;
l1s(divzeros) = m11s(divzeros).*order + m22s(divzeros).*~order;
l2s(divzeros) = m22s(divzeros).*order + m11s(divzeros).*~order;
%eigenvectors are orthogonal
if abs(e11s)>1e-15
    temp1=e21s/e11s;
    e22s=1/sqrt(1+temp1^2);
    e12s=-(temp1)*e22s;
else
    temp1=e11s/e21s;
    e12s=1/sqrt(1+temp1^2);
    e22s=-(temp1)*e12s;
end
u=[e11s e12s; e21s e22s];
lmd=[l1s;l2s];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions for the cell energy
function [Nb_Nu_pi, n_f, eterm] = get_Nb_Nu_pi(lmdn, beta_, ...
        edot0, sigmax, eps_s, eps_p, ...
        kt, muT, muU0, muB0, ...
        r_Deltamu, Omega, Cmax)
    %
    % 
    %
    nphi        = length(lmdn(:));
    epstldnss   = ones(nphi, 1) * sqrt(1+1/beta_)-1;
    n_f         = lmdn./(1+epstldnss);
    
    %
    tau_00      = calc_tau(epstldnss, zeros(nphi, 1), ...
        [], [], [], [], [], ...
        edot0, sigmax, eps_s, eps_p);
    Psi         = calcPsi(epstldnss, 0, beta_, ...
        muT, muB0, r_Deltamu);
   
    muu         = Cmax*muU0;
    mub         = Psi - Omega*tau_00.*(1+epstldnss);
  
    eterm           = exp(n_f.*(muu-mub)/kt);
    Nb_Nu_pi    = mean(eterm);

end
function [tau_f, dtau_fdlmdn, ...
        dtau_fdn_f, dtau_fdeta] = calc_tau(epstldn, epstlddot, ...
        depstldndlmdndot, depstlddotdlmdn, ...
        depstldndn_f, depstlddotdn_f, depstlddotdeta, ...
        edot0, sigmax, eps_s, eps_p)
    
    nphi  = size(epstldn,1);
    tau0  = zeros(nphi,1);
    
    % this section finds the Hill's forces
    foo              = abs(epstldn);
    dfoodepstldn     = sign(epstldn);
    id1              = foo <= eps_p;
    id2              = ~id1;
    tau0(id1)        = sigmax;
    tau0(id2)        = sigmax * exp(-((foo(id2)-eps_p)/eps_s).^2);
    
    dtau0depstldn      = zeros(nphi,1);
    dtau0depstldn(id2) = tau0(id2) .* (-2/eps_s^2) .* ...
        (foo(id2)-eps_p) .* dfoodepstldn(id2);
    
    % myerf is normalized to be between 0 an 1
    [func, dfunc]    = myerf(epstlddot/edot0+2);
    dfuncdepstlddot  = dfunc/edot0;
    
    tau_f            = tau0.*func;           % tension in the single fiber
    
    if nargout ~= 1
        dtau_fdepstldn   = dtau0depstldn.*func;
        dtau_fdepstlddot = tau0.*dfuncdepstlddot;
        
        dtau_fdlmdn = dtau_fdepstldn.*depstldndlmdndot + ...
            dtau_fdepstlddot.*depstlddotdlmdn;
        dtau_fdn_f   = diag(dtau_fdepstldn.*depstldndn_f) + ...
            diag(dtau_fdepstlddot)*depstlddotdn_f;
        dtau_fdeta   = diag(dtau_fdepstlddot)*depstlddotdeta;
    end
end
function [erf, derfdx] = myerf(x)
    
    p   = 0.47047;
    a_1 = 0.3480242;
    a_2 = -0.0958798;
    a_3 = 0.7478556;
    
    signx = sign(x);
    x = abs(x);
    expx2 = exp(-x.^2);
    
    t = 1./(1+p.*x);
    
    func1    = a_1.*t + a_2.*t.^2 + a_3.*t.^3;
    dfunc1dt = a_1 + 2.*a_2.*t + 3.*a_3.*t.^2;
    dtdx     = -p./(1+p.*x).^2;
    
    erf      = signx.*(1 - func1 .* expx2);
    
    derfdx   = -expx2.*dfunc1dt.*dtdx + ...
        2.*x.*func1.*expx2;
    
    % normalise between 0 and 1
    erf     = erf/2+1/2;
    derfdx  = derfdx/2;
end
function [Psi, dPsidepstldn, ...
        d2Psidepstldn2] = calcPsi(epstldn, eta, beta_, ...
        muT, muB0, r_Deltamu)
    
    nphi           = size(epstldn, 1);
    Psi            = zeros(nphi, 1);
    dPsidepstldn   = zeros(nphi, 1);
    d2Psidepstldn2 = zeros(nphi, 1);
    
    etamax = max(eta);
    etamin = min(eta);
    if etamax == 0
        kappa = 0;
    else
        kappa = (etamax-etamin)/etamax;
    end
    DeltamuB0      = (muT-muB0)*r_Deltamu;
    treshold       = muT - DeltamuB0;
    foo            = muB0*(1+beta_*epstldn.^2) + (1-kappa)*DeltamuB0;
    id1            = foo < treshold;
    id2            = ~id1;
    
    Psi(id1)            = foo(id1);
    Psi(id2)            = treshold;
    
    dPsidepstldn(id1)   = 2*muB0*beta_*epstldn(id1);
    dPsidepstldn(id2)   = 0;
    
    d2Psidepstldn2(id1) = 2*muB0*beta_;
    d2Psidepstldn2(id2) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Nuhat_p = bisection(Nuhat_a, Nuhat_b, etamax, eterml, nhatl, ...
    Area0Ell, numEl, delphi, id_cyto)

    %-------------------------------------------------------------------
    %   Search for Nuhat_p in the interval [Nuhat_a, Nuhat_b) using the
    %   method of bisection.
    %
    %   Modified the code published by Brato Chakrabarti on 12 Nov 2011
    %-------------------------------------------------------------------

    Nuhat_p = Nuhat_a;
    aterm = zeros(numEl,1);
    eterml_pi=eterml/pi;
    nhatlmetamax=etamax*nhatl;

    for jjj=1:1:numEl       
        if (~max(isnan(eterml_pi(:,jjj)))) 
            
            % integrate over phi for one RVE 
            aterm(jjj,1) = delphi*sum((Nuhat_a*eterml_pi(:,jjj))./ ...
                (1+Nuhat_a*eterml_pi(:,jjj)./nhatlmetamax(:,jjj)));    
        end
    end
    
    % integrate over all RVEs in a given walk
    f_a = Nuhat_a + sum(aterm.*Area0Ell)/sum(id_cyto.*Area0Ell) - 1; 


    for jjj=1:1:numEl      
        if (~max(isnan(eterml_pi(:,jjj))))    
            aterm(jjj,1) = delphi*sum((Nuhat_b*eterml_pi(:,jjj))./ ...
                (1+Nuhat_b*eterml_pi(:,jjj)./nhatlmetamax(:,jjj))); 
        end
    end
    f_b = Nuhat_b + sum(aterm.*Area0Ell)/sum(id_cyto.*Area0Ell) - 1;    


    if f_a*f_b>0 
        disp('Try different initial guesses\n')
    else
        Nuhat_p = (Nuhat_a + Nuhat_b)/2;
        for jjj=1:1:numEl 
            if (~max(isnan(eterml_pi(:,jjj))))    
                aterm(jjj,1) = delphi*sum((Nuhat_p*eterml_pi(:,jjj))./ ...
                    (1+Nuhat_p*eterml_pi(:,jjj)./nhatlmetamax(:,jjj))); 
            end
        end
        f_p = Nuhat_p + sum(aterm.*Area0Ell)/sum(id_cyto.*Area0Ell) - 1;   

        err = abs(f_p);
        while err > 1e-6
           if f_a*f_p<0 
               Nuhat_b = Nuhat_p;
           else
               Nuhat_a = Nuhat_p;          
           end
            Nuhat_p = (Nuhat_a + Nuhat_b)/2; 
            for jjj=1:1:numEl
                if (~max(isnan(eterml_pi(:,jjj))))    
                    aterm(jjj,1) = delphi*sum((Nuhat_p*eterml_pi(:,jjj))./ ...
                        (1+Nuhat_p*eterml_pi(:,jjj)./nhatlmetamax(:,jjj))); 
                end
            end
            f_p = Nuhat_p + sum(aterm.*Area0Ell)/sum(id_cyto.*Area0Ell) - 1;   
            err = abs(f_p);
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cell_tips, cent_Tria, AreaNodal, trcR, tot_trac] = ...
    calc_traction(sigcpel, thickness, nodescell, Tria, ldest)

cell_tips=zeros(1,4);
cell_tips(1) = min(nodescell(:,1));
cell_tips(2) = max(nodescell(:,1));
cell_tips(3) = min(nodescell(:,2));
cell_tips(4) = max(nodescell(:,2));

nnodcell=size(nodescell,1);
ndofscell=3*nnodcell;

nTria=size(Tria,1);

trcR=zeros(nnodcell,1);
AreaNodal=zeros(nnodcell,1);

cent_Tria=zeros(nTria,2);
fn=zeros(ndofscell,1);
anod=zeros(ndofscell,1);
x=zeros(1,3);
y=zeros(1,3);
totf=0;
totA=0;

for nel=1:nTria
    Elcor=nodescell(Tria(nel,:),:);
    x(1)=Elcor(1,1);
    y(1)=Elcor(1,2);
    x(2)=Elcor(2,1);
    y(2)=Elcor(2,2);
    x(3)=Elcor(3,1);
    y(3)=Elcor(3,2);
    cent_Tria(nel,1)=(x(1)+x(2)+x(3))/3;
    cent_Tria(nel,2)=(y(1)+y(2)+y(3))/3;
    A=abs(x(1)*(y(2)-y(3))+x(2)*(y(3)-y(1))+x(3)*(y(1)-y(2)))/2.0;
    BmatT=zeros(6,3);
    BmatT(1,1)=y(2)-y(3);
    BmatT(3,1)=y(3)-y(1);
    BmatT(5,1)=y(1)-y(2);
    BmatT(2,2)=x(3)-x(2);
    BmatT(4,2)=x(1)-x(3);
    BmatT(6,2)=x(2)-x(1);
    BmatT(1,3)=BmatT(2,2);
    BmatT(2,3)=BmatT(1,1);
    BmatT(3,3)=BmatT(4,2);
    BmatT(4,3)=BmatT(3,1);
    BmatT(5,3)=BmatT(6,2);
    BmatT(6,3)=BmatT(5,1);
    BmatT=BmatT/(2*A);
    ftemp=thickness(nel)*(BmatT*sigcpel(:,nel))*A;
    ldest1=ldest(nel,1);
    ldest2=ldest(nel,2);
    ldest4=ldest(nel,4);
    ldest5=ldest(nel,5);
    ldest7=ldest(nel,7);
    ldest8=ldest(nel,8);
    
    fn(ldest1)=fn(ldest1)+ftemp(1,1);
    fn(ldest2)=fn(ldest2)+ftemp(2,1);
    fn(ldest4)=fn(ldest4)+ftemp(3,1);
    fn(ldest5)=fn(ldest5)+ftemp(4,1);
    fn(ldest7)=fn(ldest7)+ftemp(5,1);
    fn(ldest8)=fn(ldest8)+ftemp(6,1);
    Aover3=A/3;
    anod(ldest1)=anod(ldest1)+Aover3;
    anod(ldest2)=anod(ldest2)+Aover3;
    anod(ldest4)=anod(ldest4)+Aover3;
    anod(ldest5)=anod(ldest5)+Aover3;
    anod(ldest7)=anod(ldest7)+Aover3;
    anod(ldest8)=anod(ldest8)+Aover3;
    if thickness(nel)>0
        totA= totA+A;
    end
end

for nn=1:nnodcell
    nn1=3*(nn-1)+1;
    nn2=3*(nn-1)+2;
    temp=sqrt((fn(nn1))^2+(fn(nn2))^2);
    totf=totf+temp;
    trcR(nn)=temp/anod(nn1);
    AreaNodal(nn)=anod(nn1);
end

tot_trac = totf / totA;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Disp] = interx(Nodes, Boun, Disp_old, Disp, r)
    %-------------------------------------------------------------------
    %   INTERX Intersection of curves
    %   Author : NS
    %   Version: 3.0, 21 Sept. 2010
    %   Modified to find intersection of cell outline with stripe edges
    %
    %   Nodes outside the stripe are pushed back into the stripe along a
    %   line joining the current nodal disp. with nodal disp. from the
    %   previous accepted configuration
    %-------------------------------------------------------------------

    outcounter = 0; 
    inout = zeros(2,2);  
    Boun_out = zeros(size(Boun,1),1);
    P = zeros(size(Boun,1),2);
    for i=1:size(Boun,1)        
        if(abs(Nodes(Boun(i),1)+Disp(Boun(i),1))>r)
            
            % point outside stripe
            outcounter = outcounter + 1;
            Boun_out(outcounter) = Boun(i);
            
            inout(1,:) = Nodes(Boun(i),:)+Disp_old(Boun(i),:);
            inout(2,:) = Nodes(Boun(i),:)+Disp(Boun(i),:);
            
            slope_m = (inout(2,2)-inout(1,2))/(inout(2,1)-inout(1,1));            
             
            P(outcounter,1) = sign(Nodes(Boun(i),1)+Disp(Boun(i),1))*r;
            P(outcounter,2) = slope_m*P(outcounter,1) + inout(1,2) - slope_m*inout(1,1);    
        end
        
    end
    
    if(outcounter>0)  % if atleast one node is outside the stripe
        for i=1:outcounter
            Disp(Boun_out(i),:) = P(i,:) - Nodes(Boun_out(i),:); 
        end
    end
end
