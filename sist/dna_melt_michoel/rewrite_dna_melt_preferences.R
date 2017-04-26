###function dna_melt=dna_melt_preferences()

###global dna_melt
###% dna_MELT_PREFERENCES - Set preferences for dna_melt programs
###%   set model and integration parameters at will, then run this function.   

###% MODEL PARAMETERS
kB = 8.617385*10^(-5);# % Boltzmann constant [eV/K]
Ds = 0.18; #% Morse potential depth C-G [eV]
as = 6.9; #% Morse potential inv. width C-G [1/Angstrom]
Dw = 0.12; #% Morse potential depth A-T [eV]
aw = 4.2; #% Morse potential inv. width A-T [1/Angstrom]
###D =[Dw Ds Ds Dw]';
D<-t(c(Dw, Ds, Ds, Dw))
###a=[aw as as aw]';
a<-t(c(aw, as, as, aw))

r0 = 10.0;# % closed bp radius
###sigmaslide = [0.28 0.71 0.82 0.48;
   ###           1.23 1.17 1.02 0.82;
      ###        0.69 0.86 1.17 0.71;
         ###     1.09 0.69 1.23 0.28]; % slide variance
sigmaslide<-array(c(0.28, 0.71, 0.82, 0.48,
        1.23, 1.17, 1.02, 0.82,
        0.69, 0.86, 1.17, 0.71,
        1.09, 0.69, 1.23, 0.28), dim=c(4,4))


###K = 0.1*sigmaslide.^(-1); % anharmonic stacking strength [eV]
K = 0.1*sigmaslide^(-1);#% anharmonic stacking strength [eV]
###%K = (sum(sum(K))/16)*ones(4);
###%K = 0.65*ones(4);

###alpha = repmat(0.5, [4 4]); % anharmonic stacking inv. width [1/Angstrom]
alpha<-matrix(0.5,4,4)

###sigmatwist = [3.3 3.8 4.8 2.8;
   ###           9.5 3.7 5.3 4.8;
      ###        3.8 4.0 3.7 3.8;
         ###     6.7 3.8 9.5 3.3]; % twist variance

sigmatwist<-array(c(3.3, 3.8, 4.8, 2.8,
                     9.5, 3.7, 5.3, 4.8,
                     3.8, 4.0, 3.7, 3.8,
                     6.7, 3.8, 9.5, 3.3), dim=c(4,4))




###E =0.4*sigmatwist.^(-1); % twist energy strength [eV]
E =0.4*sigmatwist^(-1); #% twist energy strength [eV]


###%E = (sum(sum(E))/16)*ones(4);
###%E = 0.04*ones(4);

h = 3.4; #% vertical distance between bps [Angstrom]
###theta0 = (pi/180).*[35.9 32.9 34.8 32.4;
   ###                 37.4 31.9 35.1 34.8;
      ###              37.8 37.4 31.9 32.9;
         ###           30.6 37.8 37.4 35.9]; % average twist angles
theta0 = (pi/180)*array(c(35.9, 32.9, 34.8, 32.4,
                    37.4, 31.9, 35.1, 34.8,
                    37.8, 37.4, 31.9, 32.9,
                    30.6, 37.8, 37.4, 35.9), dim=c(4,4));# % average twist angles


###%theta0 = (sum(sum(theta0))/16)*ones(4);
###%theta0 = 0.60707*ones(4);


###l0 = sqrt(h^2*ones(4) + 4*r0^2*(sin(0.5*theta0)).^2); 
l0 = sqrt(h^2*matrix(1,4,4) + 4*r0^2*(sin(0.5*theta0))^2); 


###% ground state length between nucleotides on same strand [Angstrom]

###% NUMERICAL INTEGRATION PARAMETERS
ML = 36; ###% size of Gauss-Legendre grid
###% if you change ML, replace legz by "legz=legpolzeros(ML)"
###legz = [-0.997830462484086;
   ###   -0.972027691049698;
      ##  -0.948272984399508;
        #-0.917497774515659;
        #-0.879929800890397;
#        -0.835847166992475;
        #-0.785576230132207;
        #-0.729489171593557;
#        -0.668001236585521;
 #       -0.601567658135981;
  #      -0.530680285926245;
   #     -0.45586394443342;
    #    -0.377672547119689;
     #   -0.296684995344028;
      #  -0.213500892316866;
       # -0.128736103809385;
#        -0.0430181984737086;
 #       0.0430181984737086;
  ##     0.213500892316866;
     #   0.296684995344028;
    #    0.377672547119689;
      #  0.45586394443342;
       # 0.530680285926245;
      #0.601567658135981;
       # 0.668001236585521;
      #0.729489171593557;
       # 0.785576230132207;
      #  0.835847166992475;
       # 0.879929800890397;
        #0.917497774515659;
#        0.948272984399508;
 #       0.972027691049698;
  #      0.988586478902212;
   #     0.997830462484086]; % roots of Legendre polynomial P_ML

legz = rbind(-0.997830462484086,
        -0.988586478902212,
        -0.972027691049698,
        -0.948272984399508,
        -0.917497774515659,
        -0.879929800890397,
        -0.835847166992475,
        -0.785576230132207,
        -0.729489171593557,
        -0.668001236585521,
        -0.601567658135981,
        -0.530680285926245,
        -0.45586394443342,
        -0.377672547119689,
        -0.296684995344028,
        -0.213500892316866,
        -0.128736103809385,
        -0.0430181984737086,
        0.0430181984737086,
        0.128736103809385,
        0.213500892316866,
        0.296684995344028,
        0.377672547119689,
        0.45586394443342,
        0.530680285926245,
        0.601567658135981,
        0.668001236585521,
        0.729489171593557,
        0.785576230132207,
        0.835847166992475,
        0.879929800890397,
        0.917497774515659,
        0.948272984399508,
        0.972027691049698,
        0.988586478902212,
        0.997830462484086) # % roots of Legendre polynomial P_ML

######% $$$ ML = 70;
######% $$$ legz = [ -0.999418285973576
######% $$$          -0.99693625196168
######% $$$          -0.99247605521169
######% $$$         -0.986045558070399
######% $$$         -0.977657405957592
######% $$$         -0.967328223664986
######% $$$         -0.955078509114293
######% $$$         -0.940932579003815
######% $$$         -0.924918516897934
######% $$$         -0.907068116260923
######% $$$         -0.887416816863348
######% $$$         -0.866003634213859
######% $$$          -0.84287108199898
######% $$$         -0.818065087625441
######% $$$         -0.791634901007893
######% $$$           -0.7636329967719
######% $$$         -0.734114970060943
######% $$$         -0.703139426151529
######% $$$         -0.670767864094077
######% $$$         -0.637064554609778
######% $$$         -0.602096412485356
######% $$$         -0.565932863718808
######% $$$         -0.528645707679711
######% $$$         -0.490308974557637
######% $$$         -0.450998778381648
######% $$$         -0.410793165902631
######% $$$         -0.369771961638462
######% $$$         -0.328016609389643
######% $$$         -0.285610010540038
######% $$$         -0.242636359463741
######% $$$         -0.199180976364858
######% $$$          -0.15533013788207
######% $$$         -0.111170905794299
######% $$$        -0.0667909541675513
######% $$$        -0.0222783952861403
######% $$$         0.0222783952861403
######% $$$         0.0667909541675513
######% $$$          0.111170905794299
######% $$$           0.15533013788207
######% $$$          0.199180976364858
######% $$$          0.242636359463741
######% $$$          0.285610010540038
######% $$$          0.328016609389643
######% $$$          0.369771961638462
######% $$$          0.410793165902631
######% $$$          0.450998778381648
######% $$$          0.490308974557637
######% $$$          0.528645707679711
######% $$$          0.565932863718808
######% $$$          0.602096412485356
######% $$$          0.637064554609778
######% $$$          0.670767864094077
######% $$$          0.703139426151529
######% $$$          0.734114970060943
######% $$$            0.7636329967719
######% $$$          0.791634901007893
######% $$$          0.818065087625441
######% $$$           0.84287108199898
######% $$$          0.866003634213859
###% $$$          0.887416816863348
###% $$$          0.907068116260923
###% $$$          0.924918516897934
###% $$$          0.940932579003815
###% $$$          0.955078509114293
###% $$$          0.967328223664986
###% $$$          0.977657405957592
###% $$$          0.986045558070399
###% $$$           0.99247605521169
###% $$$           0.99693625196168
###% $$$          0.999418285973576];

legw <- rbind( 
                c(             2*(1-t(legz)^2)/((ML+1)*t(legz)*legpol(legz,ML))            ),
                     - c(             (ML+1)*legpol(legz,ML+1)^2         )  
               )###; % weights


###legw = 2.*(1-legz.^2)./((ML+1).*legz.*legpol(legz,ML) ...
   ###                     - (ML+1).*legpol(legz,ML+1)).^2; % weights


inta <-  9.3
intb <- 40####;  % integration interval [inta,intb]
xi <- (intb-inta)/2*legz + (intb+inta)/2##;   % legendre nodes in [inta,intb]


###inta = 9.3;
###intb = 40;  % integration interval [inta,intb]
###xi = (intb-inta)/2.*legz + (intb+inta)/2;   % legendre nodes in [inta,intb]

MC <- 24#; % size of Gauss-Chebyshev grid
mc <- 1:MC#;
chebz = cos(pi.*(mc-0.5)./MC)'; % roots of Chebyshev polynomial T_MC
chebw = pi/MC; % weights 

###MC = 24; % size of Gauss-Chebyshev grid
###mc=1:MC;
###chebz = cos(pi.*(mc-0.5)./MC)'; % roots of Chebyshev polynomial T_MC
###chebw = pi/MC; % weights 


% DO NOT CHANGE BELOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setpref('dna_melt',{'kB','D','a','r0','K','alpha','E','h','theta0','l0',...
%                    'ML','legw','xi','MC','chebz','chebw'},...
%                   {kB,D,a,r0,K,alpha,E,h,theta0,l0,ML,legw,xi,MC,chebz,chebw});
%global kB,D,a,r0,K,alpha,E,h,theta0,l0,ML,legw,xi,MC,chebz,chebw
dna_melt=struct ("kB",kB, "D",D, "a",a, "r0",r0, "K",K, "alpha",alpha, "E",E, "h",h,...
"theta0",theta0, "l0",l0, "ML",ML, "legw",legw, "xi",xi, "MC",MC, "chebz",chebz, "chebw",chebw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function z = legpolzeros(n)
% LEGPOLZEROS - find the zeros of the Legendre polynomial of order n
%   


#fzero from pracma R library

fzero <- function(f, x, ..., maxiter = 100, tol = .Machine$double.eps^(1/2)) {
if (!is.numeric(x) || length(x) > 2)
stop("Argument 'x' must be a scalar or a vector of length 2.")

err <- try(fun <- match.fun(f), silent = TRUE)
if (class(err) == "try-error") {
stop("Argument function 'f' not known in parent environment.")
} else {
f <- function(x) fun(x, ...)
}

zin <- NULL
if (length(x) == 2) {
if (x[1] <= x[2]) {
a <- x[1]; b <- x[2]
} else {
warning("Left endpoint bigger than right one: exchanged points.")
a <- x[2]; b <- x[1]
}
zin <- .zeroin(f, a, b, maxiter = maxiter, tol = tol)

} else {  # try to get b
a <- x; fa <- f(a)
if (fa == 0) return(list(x = a, fval = fa))
if (a == 0) {
aa <- 1
} else {
aa <- a
}
bb <- c(0.9*aa, 1.1*aa, aa-1, aa+1, 0.5*aa, 1.5*aa,
-aa, 2*aa, -10*aa, 10*aa)
for (b in bb) {
fb <- f(b)
if (fb == 0) return(list(x = b, fval = fb))
if (sign(fa) * sign(fb) < 0) {
zin <- .zeroin(f, a, b, maxiter = maxiter, tol = tol)
break
}
}
}

if (is.null(zin)) {
warning("No interval w/ function 'f' changing sign was found.")
return(list(x = NA, fval = NA))
} else {
x1 <- zin$bra[1]; x2 <- zin$bra[2]
f1 <- zin$ket[1]; f2 <- zin$ket[2]
x0 <- sum(zin$bra)/2; f0 <- f(x0)
if (f0 < f1 && f0 < f2) {
return(list(x = x0, fval = f0))
} else if (f1 <= f2) {
return(list(x = x1, fval = f1))
} else {
return(list(x = x2, fval = f2))
}
}
}


legpolzeros<-function(n) {
#global dna_melt
  ###ZR = zeros(ceil(n/2)+1,n);

ZR<-matrix(0, ceiling(n/2)+1, n)
  ###leg = inline('legpol(x)');
leg<-function(x) legpol(x)#???
  ###ZR(2,1)=1;
  ZR(2,1)<-1
  ###for l = 2:n
    ###ZR(ceil(l/2)+1,l)=1;
    ###ord = l;
    ###for k = 1+mod(l,2):ceil(l/2)
      ###ZR(k,l) = fzero(leg, [ZR(k-mod(l,2),l-1) ZR(k+1-mod(l,2),l-1)]);
    ###end
  ###end

for (l in 2:n) {
    ZR(ceiling(l/2)+1,l)<-1
    ord<-l
      for (k in (1 + l %% 2):(ceiling(l/2))) {
        ZR(k,l) <- fzero(leg, c(ZR[(k-l%%2),(l-1)], ZR[(k+1-l%%2),(l-1)])) 
                              ###[ZR(k-mod(l,2),l-1) 
                              ###                    ZR(k+1-mod(l,2),l-1)]);
      }
}
   
  if (n%%2 ==0) {
    z <- sort(t(c(-t(ZR[1:ceiling(n/2),n]), t(ZR[1:ceiling(n/2),n]))))
###  if mod(n,2) == 0
   ### z = sort([-ZR(1:ceil(n/2),n).' ZR(1:ceil(n/2),n).'].');
  } else {
    z <- sort(t(c(-t(ZR[2:ceiling(n/2),n]), 0, t(ZR[2:ceiling(n/2),n]))))
  }


###else
  ###z = sort([-ZR(2:ceil(n/2),n).' 0  ZR(2:ceil(n/2),n).'].');
  ###end


global dna_melt
  ZR = zeros(ceil(n/2)+1,n);
  leg = inline('legpol(x)');
  ZR(2,1)=1;
  for l = 2:n
    ZR(ceil(l/2)+1,l)=1;
    ord = l;
    for k = 1+mod(l,2):ceil(l/2)
      ZR(k,l) = fzero(leg, [ZR(k-mod(l,2),l-1) ZR(k+1-mod(l,2),l-1)]);
    end
  end
   
  if mod(n,2) == 0
    z = sort([-ZR(1:ceil(n/2),n).' ZR(1:ceil(n/2),n).'].');
  else
    z = sort([-ZR(2:ceil(n/2),n).' 0  ZR(2:ceil(n/2),n).'].');
  end
  
    
    
    #new function in R for nargin return
    
    
    
    nargin <- function(...){
      
      print(match.call())
      nargin <- length(as.list(match.call())) 
      print(nargin)
    }
    
    library(pracma)
    
function l = legpol(x,n)
% LEGPOL - legendre polynomial of order n evaluated at
%   x is either real value in [-1,1] or a (column) vector of such values

    
    
legpol<-function(x, n) {
  ord <- n
  if (nargin(legpol)==1) {
  ###if nargin == 1
    tmp<-legendre(ord, t(x))
  ###tmp = legendre(ord,x');
    s<-dim(x)#??
    ###s = size(x);
    if (s[1,1] ==1){
      l <-tmp[1]
    } else {
      l<-t(tmp[1,])
    }
  } else {
    tmp <- legendre(n,t(x))
    s <- size(x)
    if (s[1,1] == 1) {
      l = tmp[1]
      } else { 
      l = t(tmp[1,]) #?? non -conjugated
      }
  }
}

#alternative legpol

legpol<-function(x, n) {
  ord <- n
  if (nargs()==1) {
    ###if nargin == 1
    tmp<-legendre(ord, t(x))
    ###tmp = legendre(ord,x');
    s<-dim(x)#??
    ###s = size(x);
    if (s[1] ==1){ #??
      l <-tmp[1]
    } else {
      l<-t(tmp[1,])
    }
  } else {
    tmp <- legendre(n,t(x))
    s <- size(x)
    if (s[1] == 1) {
      l = tmp[1]
    } else { 
      l = t(tmp[1,]) #?? non -conjugated
    }
  }
}


    ###if s(1,1) == 1
      ###l = tmp(1);
    ###else
      ###l = tmp(1,:).';
    ###end
       ###else
     ###tmp = legendre(n,x');
    ###s = size(x);
    ###if s(1,1) == 1
      ###l = tmp(1);
    ###else
      ###l = tmp(1,:).';
         ###                         end
            ###                    end
                                  
                                  
                                  
  
}

  
  if nargin == 1
    tmp = legendre(ord,x');
    s = size(x);
    if s(1,1) == 1
      l = tmp(1);
    else
      l = tmp(1,:).';
    end
  else
    tmp = legendre(n,x');
    s = size(x);
    if s(1,1) == 1
      l = tmp(1);
    else
      l = tmp(1,:).';
    end
  end
  















Legendre <- function(x, n, normalized=TRUE, intercept=FALSE, rescale=TRUE)
{
  ## Create a design matrix for Legendre polynomials
  ## x - numeric
  ## n - see orthopolynom
  ## normalized - logical, see orthopolynom
  ## intercept - logical, add intercept
  tmp <- legendre.polynomials(n=n, normalized=normalized)
  if(!intercept) tmp <- tmp[2:length(tmp)]
  polynomial.values(polynomials=tmp, x=x, matrix=TRUE)
}