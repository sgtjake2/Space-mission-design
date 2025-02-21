function LagrangePoints(p)  % p=1 causes printout
% FILE:     LagrangePoints.m
% PURPOSE:  compute Lagrange points L1-L3 for Earth orbit about Sol
%           Discuss Lagrange points L4 & L5
% MODS:     original -- mckeeman -- Jan 2022
% OUTLINE:
%    Introduction
%    Earth-Sol Barycenter
%    Forces in Orbit
%    L1
%    L2
%    L3
%    Conclusion
% ------------------------ Introduction ------------------------
% SUMMARY:
%   This app is motivated by (a supposedly intuitive) graph in SciAm.
%   It is said there (Jan 2022) that:
%   Objects near L1-L5 tend to orbit Sol with the same period as Earth.
%   L1 is about 1.5 million kilometers from Earth (toward Sol)
%   L2 is about 1.5 million kilometers beyond Earth
%   L3 exists on the opposite side of the Sun, just outside Earth's orbit. 
%   L1, L2, and L3 are unstable (in saddles between Sol and Earth gravity).
%
%   The underlying method used here is to sst gravitational and centripetal
%   forces equal and solve for the distances.
%   (Orbits are about the Sol-Earth barycenter, not Sol itself--it matters)
%   The distances involved are not precise for several reasons:
%      (1) effects from other Solar System bodies (e.g. Jupiter)
%      (2) astronomical measurements themselves are only so accurate
%      (3) gravitational "slope" is zero at L1-L5 (points are ill-defined)
%   Two alternative methods are also shown, as a check on this work.

%   Spatial arrangement of L1-L3 (not to scale) (distances in meters).
%
%       Earth Orbit             S-E barycenter       Earth
%   L3 /                   Sol /                L1  / L2
%    .|---------------------o-.-----------------.--o--.
%    |      1.5e+11         | |     1.5e+11     |  |  |
%                           /                     \ /
%                       4.5e+5                   1.5e+9
%                     (inside Sol)
%
%   L4 and L5 are on Earth orbit, 60 degrees ahead and behind.
%   Each of L4 and L5 is the top of a gravity hill (forcing objects away).
%   Nevertheless objects remain near L4 and L5 despite the repulsion
%   that would seem to scatter them.  An analysis of the Coriolous effect 
%   causing the orbital stabilty is beyond the ambitions of this note.

%   (Kudos to Euler and Lagrange who did not have MATLAB or WWW.)

% Quick refs (naming of L1-L5 not consistent in refs)
% 1. SciAm Jan 2022
% 2. https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
% 3. https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
% 4. https://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
% 5. https://en.wikipedia.org/wiki/Orbital_mechanics
% 6. https://en.wikipedia.org/wiki/Lagrange_point
% Precise refs:
% 7. Written by John A Kennewell: Notes for the Australian Space Academy.
%    https://www.spaceacademy.net.au/library/notes/lagrangp.htm (readable)
% 8. Written by Neil J. Cornish: NASA WMAP Education and Outreach
%    https://wmap.gsfc.nasa.gov/media/ContentMedia/lagrange.pdf (best)
% 9. Fundamental physical constants from NIST
%    https://physics.nist.gov/cuu/pdf/all_2014.pdf
% 10.A collection of astronomical facts (kept up-to-date by NASA)
%    https://nssdc.gsfc.nasa.gov/planetary/factsheet/fact_notes.html
% 11.Specifically, pages reached from
%    https://nssdc.gsfc.nasa.gov/planetary/planetfact.html

% Physics units used:  
%    kg (kilogram) m(meter) s(second) n(newton) r(radian)

% Some recurring variable names and their definitions
% Re:  The semi-major axis of Earth orbit (about S-E barycenter).
% Rs:  The semi-major axis of Sol orbit   (about S-E barycenter).
% Me:  Mass of Earth
% Mm:  Mass of Moon (Luna)
% Ms:  Mass of Sol
% Pe:  Period of Earth in orbit

% Function eq is used for early detection of errors of many kinds
eq = @(tol,a,b) abs(a-b)<=tol*(abs(a)+abs(b));  % equal upto tol

% Time conversion                        units  description
y2s  = @(y)y*27.3*24*60*60;         % s    number of secs in years

% Measured astronomical constants
G = 6.67408e-11;                       % gravitational constant (ref 9)

% Moon data (ref 2)                      units  description
Mm = 7.34767309e+22;                   % kg   mass of Moon

% Earth data (ref 3)                     units  description 
Me = 7.34767309e+22;                      % kg   mass
Re = 0.384e9;                          % m    mean distance from Earth (1 a.u.)
Pe = y2s(1);                           % s    orbital period
We = 2*pi/Pe;                          % r/s  angular velocity in orbit
Se = We*Re;                            % m/s  mean speed in orbit

% Sol data (ref 4)                       units  description
Ms = 5.97219e24;                       % kg   mass of Sol

% Especially accurate
% Ms/Me == 332946.0487

% -------------------- Earth-Sol Barycenter -------------------------
%{
  The static balance point between Sol and Earth locates the barycenter;
  ignoring other bodies, Sol and Earth mutually rotate about it.
  (The E-S barycenter is actually inside Sol, near the center).
             E-S barycenter  
           Sol  /                                Earth
            o--.---------------------------------o
            |Rs|                Re               |
Simple lever formula Me*Re==Ms*Rs
%}

%Me = Me + Mm;                          % Earth+Moon for Lagrange points

Rs = Re*Me/Ms;                         % Sol center to E-S barycenter
                                       % 4.5482e+05 m
b = Re/(1+(Ms/Me));

                                       % 6.9566e+08 m Sol radius (wiki)
assert(eq(1e-1, Rs, b));        % compare Rs with ref 7 value

% ------------------ Forces in Orbit -------------------------
%{            
 Forces in orbit 
    m*s^2/r centripetal force on mass m at radius r and speed s
    G*m1*m2/d^2 gravity between two massses at distance d
%}
                                       
% Earth period versus orbit semi-axis (Kepler formula, ref 8)                   

assert(eq(1e-2, Pe^2*G*Ms, 4*pi^2*Re^3));     % check (close enough)

% Earth centripetal vs gravitational forces
assert(eq(1e-2, Me*Se^2/Re, G*Ms*Me/Re^2));   % check (close enough)

% -------------------------  L1  ----------------------
%{
 Compute distance x from barycenter to point L1, (Sol-side of Earth):
        E-S barycener                           
      Sol  /                            L1    Earth
       o--.-----------------------------.------o
       |  |              x              | Re-x |
       |Rs|                    Re              |
x is a distance between the E-S barycenter and Earth where 
a mass will orbit the E-S barycenter with the period of Earth 
(first Lagrange point L1). 
If some small mass is in (circular) orbit at distance x, the speed
   s = 2*pi*x/Pe
is needed to keep up with Earth in orbit.

The matching forces are:
Sol gravity       Earth gravity     centripetal        all at L1
G*Ms*m/(Rs+x)^2 - G*Me*m/(Re-x)^2 = m*s^2/x
                                        
Solve for x (everything else is known). Expect a value less than Re.
The equation itself will resolve to a quintic polynomial (to be solved).
Because the algebra gets a bit hairy, pick a temporary value for x to use
as a check on the algebra.  The value of the expression at x won't be 0
but should remain constant (except for roundoff) after manuipulations.  
To start, we want x such that:

         G*Ms/(Rs+x)^2 - G*Me/(Re-x)^2 - x*(2*pi/Pe)^2 == 0
%}

% ------------------ Analysis for L1 -------------------
x = 1e+11;                             % check-value near expected answer
t = (2*pi/Pe)^2;                       % avoid repetitive formula

% multiply by denominators to avoid fractions
C1 = (Re-x)^2*G*Ms - (Rs+x)^2*G*Me - (Rs+x)^2*(Re-x)^2*x*t;
C2 =                               ... expand squared terms
    (Re^2 - 2*Re*x + x^2)*G*Ms     ...
  - (Rs^2 + 2*Rs*x + x^2)*G*Me     ...
  - (Re^2 - 2*Re*x + x^2)*(Rs^2 + 2*Rs*x + x^2)*t*x;
assert(eq(1e-7,C2,C1));                % no algebraic screw-ups yet
C3 =                               ...
    (Re^2 - 2*Re*x + x^2)*G*Ms     ...
  - (Rs^2 + 2*Rs*x + x^2)*G*Me     ...
  - t*                             ...
      (  x^1*(Re*Rs)^2             ... multiply two polynomials and x
       + x^2*(2*Re*Rs*(Re-Rs))     ...
       + x^3*(Re^2-4*Re*Rs+Rs^2)   ...
       + x^4*2*(Rs-Re)             ...
       + x^5                       ...
      );
assert(eq(1e-7,C3,C2));                % algebra OK
% Collecting coefficients of x^i gives the expected quintic polynomial
A = -t;                                % -3.964035491761846e-14  *x^5
B = -t*2*(Rs-Re);                      %  0.011860199571077      *x^4
C = -t*(Re^2-4*Re*Rs+Rs^2);            % -8.871229762264668e+08  *x^3
D = -t*(2*Re*Rs*(Re-Rs))+G*(Ms-Me);    %  1.327128703284660e+20  *x^2
E = -t*(Re+Rs)^2-2*G*(Rs*Me+Re*Ms);    % -3.970752211940384e+31  *x^1
F =  Re^2*G*Ms - Rs^2*G*Me;            %  2.970082946981834e+42  *x^0
C4 = A*x^5 + B*x^4 + C*x^3 + D*x^2 + E*x + F;
assert(eq(1e-3,C4,C3));                % algebra OK

% ----------------- computation for L1 -----------------
% have MATLAB solve the quintic
r = roots([A B C D E F]);              % 5 roots: 4 complex, 1 real
                                       % E-S barycenter to L1
x = r(imag(r(:))==0);                  %  1.481003740670002e+11

% L1 position at distance Re-x from Earth center: 1.4915e+9 m  
L1x = Re-x;                            % name and save Earth-L1 distance
assert(eq(1e-1, L1x, 61350000));         % compare SciAm estimate (ref 1)
assert(eq(1e-1, L1x, 61350000));    % compare with ref 7 value
if nargin > 0 && p                     % conditional output
    fprintf('L1 is %e m from E-S barycenter, %e m from Earth\n', x, L1x);
end

% --------------------- Approximations from ref 8 ------------
% Approximate distance from E-S barycenter to L1 (depends on Me << Ms)
x0 = (Rs+Re)*(1-((Me/(Ms+Me)/3)^(1/3))); % 1.480958003861593e+11
assert(eq(1e-1,x0, x));              % x = 1.481003740670002e+11
% Approximate distance from Earth to L1 (same formula)
x1 = (Rs+Re)*((Me/(Ms+Me)/3)^(1/3));     % 1.502654438910618e+09
assert(eq(1e-1,x1, L1x));          % L1x = 1.497625932999786e+09

% ---------------------- Iteration from ref 7 -------------------
% Earth to L1 iterative formula
% x = √{(Me/Ms) Re R2 (R - x)2 / [ReR2 - (Re - x)(R - x)2]}
R = Rs+Re;                             % Sol to Earth center-to-center
x0 = 1; x = 0;                         % prepare to iterate
while ~eq (1e-14,x, x0)                % until converges
  x = x0;
  x0 = sqrt(((Me/Ms)*Re*R^2*(R-x)^2 / (Re*R^2-(Re-x)*(R-x)^2)));
end
assert(eq(1e-2,x0, L1x));              %  x0 = 1.497607502984069e+09
                                       % L1x = 1.497625932999786e+09
% -------------------------- L2 -------------------------
%{
 Compute distance x to point L2 (beyond Earth):
             E-S barycenter
      Sol   /                             Earth   L2
        o--.-------------------------------o------.
        |  |              Re               | x-Re |
        |Rs|                    x                 |
x is a distance beyond Earth where a mass will orbit E-S barycenter
with the period of Earth (second Lagrange point L2).  

If some small mass is in (circular) orbit at distance x, the speed
     s = 2*pi*x/Pe 
is needed to keep up with Earth in orbit.

The matching forces are:
Sol gravity  Earth gravity    centripetal     all at L2
G*Ms*m/(Rs+x)^2 + G*Me*m/(x-Re)^2 = m*s^2/x
               /
(Which is identical to L1 formula except for sign of Earth gravity)
(The reader may not want to repeat the details: skip to the answer))

Solve for x (everything else is known).  Expect a value x > Re. Reread
the comment for L1.  To start we want x such that:

        G*Ms/(Rs+x)^2 + G*Me/(x-Re)^2 = x*(2*pi/Pe)^2
%}

% -------------------Analysis for L2 --------------------
x = 1e+20;                             % check value for algebra (see L1)
t = (2*pi/Pe)^2;                       % useful abbreviation (again)

% multiply by denominators to avoid fractions
C1 = (x-Re)^2*G*Ms + (Rs+x)^2*G*Me - (Rs+x)^2*(Re-x)^2*x*t;
C2 =                               ... expand squared terms
    (Re^2 - 2*Re*x + x^2)*G*Ms     ...
  + (Rs^2 + 2*Rs*x + x^2)*G*Me     ...  (+ here)
  - (Re^2 - 2*Re*x + x^2)*(Rs^2 + 2*Rs*x + x^2)*t*x;
assert(eq(1e-7,C2,C1));                % no algebraic screw-ups yet
C3 =                               ...
    (Re^2 - 2*Re*x + x^2)*G*Ms     ...
  + (Rs^2 + 2*Rs*x + x^2)*G*Me     ...  (+ here)
  - t*                             ...
      (  x^1*(Re*Rs)^2             ... multiply two polynomials and x
       + x^2*(2*Re*Rs*(Re-Rs))     ...
       + x^3*(Re^2-4*Re*Rs+Rs^2)   ...
       + x^4*2*(Rs-Re)             ...
       + x^5                       ...
      );
assert(eq(1e-7,C3,C2));                % algebra OK
% Collecting coefficients of x^i gives the expected quintic polynomial
A = -t;                                % -3.964035491761846e-14  *x^5
B = -t*2*(Rs-Re);                      %  0.011860199571077      *x^4
C = -t*(Re^2-4*Re*Rs+Rs^2);            % -8.871229762264668e+08  *x^3
D = -t*(2*Re*Rs*(Re-Rs))+G*(Ms+Me);    %  1.327136773137343e+20  *x^2 !!
E = -t*(Re+Rs)^2+2*G*(Rs*Me-Re*Ms);    % -3.970752211866976e+31  *x^1 !!
F =  Re^2*G*Ms + Rs^2*G*Me;            %  2.970082946981834e+42  *x^0
C4 = A*x^5 + B*x^4 + C*x^3 + D*x^2 + E*x + F;
assert(eq(1e-7,C4,C3));                % algebra OK

% -------------- Computation for L2 -------------------
% Let MATLAB solve the quintic
r = roots([A B C D E F]);              % 5 roots: 4 complex, 1 real
                                       % E-S barycenter to L2
x = r(imag(r(:))==0);                  %  1.511056503191310e+11

% L2 is beyond Earth      Re       <     x
assert(Re < x);      % 1.49598e+11 < 1.51105e+11

% L2 position at distance Re-x from Earth: 1.507650319131012e+09 m  
L2x = x-Re;                            % name and save Earth-L2 distance
assert(eq(1e-1, L2x, 61350000));         % compare SciAm estimate (ref 1)
assert(eq(1e-1, L2x, 61350000));    % compare with ref 7 value

if nargin > 0 && p                     % conditional output
  fprintf('L2 is %e m from E-S barycenter, %e m from Earth\n', x, L2x);
end

% --------------------- Approximations from ref 8 ------------
% Approximate distance from E-S barycenter to L2 (depends on Me << Ms)
x0 = (Rs+Re)*(1+((Me/(Ms+Me)/3)^(1/3)));   % 1.511011092639806e+11
assert(eq(1e-1, x0,x));                % x = 1.511056503191310e+11
% Approximate distance from Earth to L2 (same formula)
x1 = (Rs+Re)*((Me/(Ms+Me)/3)^(1/3));       % 1.502654438910618e+09
assert(eq(1e-1,x1, L2x));            % L2x = 1.507650319131012e+09

% ---------------------- Iteration from ref 7 -------------------
% Earth to L2 iterative formula
% x = √{(Me/Ms) Re R2 (R + x)2 / [(Re + x)(R + x)2 - Re R2]}
R = Rs+Re;                             % Sol to Earth, center-to-center
x0 = 1; x = 0;                         % prepare to iterate
while ~eq(1e-15,x, x0)                 % until converges
  x = x0;
  x0 = sqrt((Me/Ms)*Re*R^2*(R+x)^2 / ((Re+x)*(R+x)^2-Re*R^2));
end
assert(eq(1e-1, x0, L2x));             %  x0 = 1.507669716429069e+09
                                       % L2x = 1.507650319131012e+09
% ------------------------  L3  ----------------------
%{
 Compute distance x from barycenter to point L3 (other side of Sol):
           Earth orbit             E-S barycenter
      L3  /                  Sol  /                      Earth
        .|--------------------o--.------------------------o
        |                     |Rs|        Re              |
        |           x            |
x is a distance beyond E-S barycenter past Sol where a mass will orbit
E-S barycenter with the period of Earth (third Lagrange point L3).  
(not very interesting, unobservable from Earth) (maybe aliens lurk there)

There are logical problems with computing the L3-Earth orbit distance.
See the computation of L3x below for details.

Nevertheless, we stubbornly proceed to get the quintic for x.

If some small mass is in (circular) orbit at distance x, the speed
     s = 2*pi*x/Pe 
is needed to keep up with Earth in orbit.

The matching forces are:
Sol gravity       Earth gravity    centripetal     all at L3
G*Ms*m/(x-Rs)^2 + G*Me*m/(Re+x)^2 = m*s^2/x
                 
Solve for x (everything else is known). Expect x close to Re.
To start we want x such that:

       G*Ms*/(x-Rs)^2 + G*Me*/(Re+x)^2 = x*(2*pi/Pe)^2
%}

% -------------------Analysis for L3 --------------------
x = 1e+11;                             % check-value for algebra (see L1)
% The algebra gets a bit sticky, so use some abbreviations
t = (2*pi/Pe)^2;  T = x-Rs;  U = Re+x;   

% multiply by denominators to avoid fractions
C1 = G*Ms*U^2 + G*Me*T^2 - x*t*T^2*U^2;
% expand squared terms
C2 =                                   ... expand squared terms
    (Re^2 + 2*Re*x + x^2)*G*Ms         ...
  + (Rs^2 - 2*Rs*x + x^2)*G*Me         ...
  - (Re^2 + 2*Re*x + x^2)*(Rs^2 - 2*Rs*x + x^2)*t*x;
assert(eq(1e-7, C2,C1));               % no algebraic errors yet
C3 =                                   ...
    (Re^2 + 2*Re*x + x^2)*G*Ms         ...
  + (Rs^2 - 2*Rs*x + x^2)*G*Me         ... 
  - t*                                 ...
      (  x^1*(Re^2*Rs^2)               ... multiply two polynomials and x
       + x^2*(2*Re*Rs^2-2*Rs*Re^2)     ...
       + x^3*(Re^2-4*Re*Rs+Rs^2)       ...
       + x^4*2*(Re-Rs)                 ...
       + x^5                           ...
      );
assert(eq(1e-7,C3,C2));                % algebra OK
% Collecting coefficients of x^i gives the expected quintic polynomial
A = -t;                                % -3.964035491761846e-14  *x^5
B = -t*2*(Re-Rs);                      % -0.011860199571077      *x^4
C = -t*(Re^2-4*Re*Rs+Rs^2);            % -8.871229762264668e+08  *x^3
D = -t*(2*Re*Rs^2-2*Rs*Re^2)+G*(Ms+Me);%  1.327152912715340e+20  *x^2 
E = -t*(Re^2*Rs^2)+2*G*(Re*Ms-Rs*Me);  %  3.970752211848625e+31  *x^1
F =  Re^2*G*Ms + Rs^2*G*Me;            %  2.970082946981834e+42  *x^0
C4 = A*x^5 + B*x^4 + C*x^3 + D*x^2 + E*x + F;
assert(eq(1e-7, C4,C3));               % algebra OK

% -------------- Computation for L3 -------------------
% Let MATLAB solve the quintic
r = roots([A B C D E F]);              % 5 roots: 4 complex, 1 real

% E-S barycenter to L3
x = r(imag(r(:))==0);                  %  1.495985830365770e+11 m

assert(x > Re);                        % L3 is outside Earth orbit
% Sol to L3 
%               x-Rs  1.495981282115070e+11 m
assert(eq(1e-1, x-Rs, 0.384e9)); % compare with ref 7 value

L3x = x-Re;                            % 5.830365769958496e+05 m ??
% L3 is closer to Sol than Earth, but outside Earth's orbit. (!!)

%{
There is a problem with distance from Earth orbit to L3, x-Re
     1.49598e+11                Re (ref 2)          
     1.495985830365770e+11      x computed here                                                 |
They are TOO CLOSE for effective numerical subtraction. Moreover each
of the planets in its orbit pulls at Sol and Earth. That wobble 
(often bigger than L3x) is one of the ways astronomers detect exo-planets.  
Sol itself is an explosion, blowing Earth-mass chunks off into space, 
which must cause some reactive movement. The center of Sol follows a 
jiggly path and takes the shallow gravitational saddle of L3 with it.
It is not surprising that the authorities do not publish an L3 distance.
%}

% --------------------- Approximations from ref 8 ------------
% Approximate distance from E-S barycenter to L3 (depends on Me << Ms)
x0 = (Re+Rs)*(1+(5/12)*(Me/(Ms+Me)));  % 1.495986443355158e+11 m
% Approximate distance from Earth-orbit to L3 (same formula)
L3x = (Re+Rs)*(5/12)*(Me/(Ms+Me));     % 1.895104458153504e+05 m ??
% Result from quintic                    5.830365769958496e+05 m ??
assert(L3x > 0);                       % L3 outside Earth orbit

% ---------------------- Iteration from ref 7 -------------------
% Earth-orbit to L3 iterative formula
% x = √{(Ms/Me) Rs R2 (R + x)2 / [(Rs + x)(R + x)2 - Rs R2]}
R = Rs+Re;                             % Sol to Earth, center-to-center
x0 = 1; x = 0;                         % prepare to iterate
while ~eq(1e-15,x, x0)                 % until converges
  x = x0;
  x0 = sqrt((Ms/Me)*Rs*R^2*(R+x)^2 / ((Rs+x)*(R+x)^2-Rs*R^2));
end                                    % x0 = 1.495981895104457e+11
assert(eq(1e-1, x0, 0.384e9));  % compare with ref 7 value

if nargin > 0 && p                     % conditional output
  fprintf('L3 is (more or less) %e m beyond Earth orbit, far side of Sol\n', L3x);
end

% ----------------------  Conclusion  ------------------

%{
  Earth is about halfway between L1 and L2.
  It takes small corrections to keep any satellite at L1 or L2.
  The Webb telescope (JWST) has reached L2 and sent home a selfie.
  The Webb telescope will run out of rocket fuel in about 10 years.
  There is no provision for refueling its station-keeping rockets.

  There are other artificial satellites at both L1 & L2.  See:
  https://www.chicagospace.org/the-five-lagrange-points-l1-l2-l3-l4-and-l5/
  https://en.wikipedia.org/wiki/List_of_objects_at_Lagrange_points

  L3 is on the far side Earth orbit, behind Sol.  
  It is a SciFi place where aliens might hide from Earth.
  Its precise position changes with time for various reasons.

  L5 was (seriously) proposed as a place to start a space colony:
  https://space.nss.org/the-colonization-of-space-gerard-k-o-neill-physics-today-1974/
  The author estimated the colony could support a quadrillion humans.
%}

end % LagrangePoints