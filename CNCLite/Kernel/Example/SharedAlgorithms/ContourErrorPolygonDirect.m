function con = ContourErrorPolygonDirect( pts, ptsPolygon )
% Calculate the contour error from pts to a polygon.
% The polygon is defined by the ptsPolygon, and is not necessary convex.
% This function is different from the contourErrorPolygon() function in the
% sense that it calculates the difference from the point to polygon directly. 
% input:
%   pts, 3-D points, of which the contour errors need to be calculated.
%       each row is a point.
%   ptsPolygon, 3-D points defining the polygon.
%       each row is a vertex of the polygon
% output:
%   con, contour errors.
% HJ, 20161221

len1 = size(pts, 1); % number of points
len2 = size(ptsPolygon, 1); % number of polygon points
con = zeros(len1, 1);
flag = norm(ptsPolygon(1, :) - ptsPolygon(len2, :) ); % whether the curve is closed.
if 0 == flag
    numberSegment = len2 - 2; % if the polygon is closed, the last segment is null.
else
    numberSegment = len2 - 1;
end
for ii=1:len1
    pt = pts(ii, :);
    dis = norm(pt - ptsPolygon(1, :) ); % an initial distance.
    for jj = 1:numberSegment
        d = distancePoint2Line(pt, ptsPolygon(jj, :), ptsPolygon(jj+1, :) );
        if d < dis
            dis = d;
        end
    end
    con(ii) = dis;
end
end

function dis = distancePoint2Line(pt, pt0, pt1)
% calculate the distance from pt to a line start from pt0 to pt1.
% only applied to 3-D situation.
vec1 = pt - pt0;
vec2 = pt1 - pt0;
dis = norm(cross(vec1, vec2)) / norm(vec2);
end