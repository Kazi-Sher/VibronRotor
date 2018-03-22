% Copyright (C) 2004-2011 David Legland <david.legland@grignon.inra.fr>
% Copyright (C) 2004-2011 INRA - CEPIA Nantes - MIAJ (Jouy-en-Josas)
% Copyright (C) 2016 Adapted to Octave by Juan Pablo Carbajal <ajuanpi+dev@gmail.com>
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%     1 Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%     2 Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS''
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%INTERSECTPOLYLINES Find the common points between 2 polylines
%
%   INTERS = intersectPolylines(POLY1, POLY2)
%   Returns the intersection points between two polylines. Each polyline is
%   defined by a N-by-2 array representing coordinates of its vertices: 
%   [X1 Y1 ; X2 Y2 ; ... ; XN YN]
%   INTERS is a NP-by-2 array containing coordinates of intersection
%   points.
%
%   INTERS = intersectPolylines(POLY1)
%   Compute self-intersections of the polyline.
%
%   Example
%   % Compute intersection points between 2 simple polylines
%     poly1 = [20 10 ; 20 50 ; 60 50 ; 60 10];
%     poly2 = [10 40 ; 30 40 ; 30 60 ; 50 60 ; 50 40 ; 70 40];
%     pts = intersectPolylines(poly1, poly2);
%     figure; hold on; 
%     drawPolyline(poly1, 'b');
%     drawPolyline(poly2, 'm');
%     drawPoint(pts);
%     axis([0 80 0 80]);
%
%   This function is largely based on the 'interX' function, found on the
%   FileExchange:
%   https://fr.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections
%   
%   See also
%   polygons2d, polylineSelfIntersections, intersectLinePolygon

function pts = intersectPolylines(poly1, varargin)

  % Check number of inputs
  narginchk(1, 2);

  % Specific init depending on number of inputs
  if nargin == 1
      % Compute self-intersections 
      % -> Avoid the inclusion of common points
      poly2 = poly1;
      hF = @lt;
  else
      % Compute intersections between distinct lines
      poly2 = varargin{1}; 
      hF = @le;
  end

  % Get coordinates of polyline vertices
  x1 = poly1(:,1);  
  x2 = poly2(:,1)';
  y1 = poly1(:,2);  
  y2 = poly2(:,2)';

  % differentiate coordinate arrays
  dx1 = diff(x1); dy1 = diff(y1);
  dx2 = diff(x2); dy2 = diff(y2);

  % Determine 'signed distances'
  S1 = dx1 .* y1(1:end-1) - dy1 .* x1(1:end-1);
  S2 = dx2 .* y2(1:end-1) - dy2 .* x2(1:end-1);

  C1 = feval(hF, D(bsxfun(@times,dx1,y2) - bsxfun(@times,dy1,x2), S1), 0);
  C2 = feval(hF, D((bsxfun(@times,y1,dx2) - bsxfun(@times,x1,dy2))', S2'), 0)';

  % Obtain the segments where an intersection is expected
  [i, j] = find(C1 & C2); 

  % Process case of no intersection
  if isempty(i)
      pts = zeros(0, 2);
      return;
  end

  % Transpose and prepare for output
  i=i'; dx2=dx2'; dy2=dy2'; S2 = S2';
  L = dy2(j).*dx1(i) - dy1(i).*dx2(j);

  % Avoid divisions by zero
  i = i(L~=0);
  j = j(L~=0);
  L = L(L~=0);

  % Solve system of eqs to get the common points
  res = [dx2(j).*S1(i) - dx1(i).*S2(j), dy2(j).*S1(i) - dy1(i).*S2(j)] ./ [L L];
  pts = unique(res, 'rows');

  % Innre function computing a kind of cross-product
  function u = D(x,y)
      u = bsxfun(@minus, x(:,1:end-1), y) .* bsxfun(@minus, x(:,2:end), y);
  end

end

%!demo
%! poly1 = [20 10 ; 20 50 ; 60 50 ; 60 10];
%! poly2 = [10 40 ; 30 40 ; 30 60 ; 50 60 ; 50 40 ; 70 40];
%! pts = intersectPolylines(poly1, poly2);
%! figure; hold on; 
%! drawPolyline(poly1, 'b');
%! drawPolyline(poly2, 'm');
%! drawPoint(pts);
%! axis([0 80 0 80]);
%! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%! % Compute intersection points between 2 simple polylines
