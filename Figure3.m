%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Figure3.m
% 
%  Computes overapproximations of reachable set for tunnel diode oscillator
%  using algorithms described in "Reachability analysis of nonlinear
%  systems using matrix measures" by John Maidens and Murat Arcak
%
%  Used to produce parts (b) and (d) of Figure 3 from this paper
%
%  John Maidens
%  July 17, 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



flowType = 'linear';
accuracy = 'veryHigh';

[ data, g, data0 ] = reachTunnelDiode(flowType, accuracy);

