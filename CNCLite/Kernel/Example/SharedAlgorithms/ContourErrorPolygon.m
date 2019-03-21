function con = ContourErrorPolygon( pts, ptsPolygon )
% calculate the contour error from pts to a polygon.
% the polygon is defined by the ptsPolygon.
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
flag = norm(ptsPolygon(1, :) - ptsPolygon(len2, :)); % whether the curve is closed.
for ii=1:len1
    pt = pts(ii, :);
    A = repmat(pt, [len2, 1]); % repeat pts(ii, :)
    disA = A - ptsPolygon; % disMatrix
    dis = sqrt(sum(disA.*disA, 2)); % minimum distance index
    [~, index] = min(dis);
    if 1 == index
        if(flag > 0)
            con(ii) = distancePoint2Line(pts(ii, :), ptsPolygon(index, :), ptsPolygon(index+1, :));
        else
            con(ii) = contourErrorByRegion(pts(ii, :), ptsPolygon(len2-1, :), ptsPolygon(index, :),...
                ptsPolygon(index+1, :)); % for closed curve, ptsPolygon(len2, :) = ptsPolygon(1, :)
        end
    elseif len2 == index
        if(flag > 0)
            con(ii) = distancePoint2Line(pts(ii, :), ptsPolygon(index-1, :), ptsPolygon(index, :));
        else
            con(ii) = distancePoint2Line(pts(ii, :), ptsPolygon(index-1, :), ptsPolygon(index, :),...
                ptsPolygon(2, :)); % for closed curve, ptsPolygon(len2, :) = ptsPolygon(1, :)
        end
    else
        con(ii) = contourErrorByRegion(pts(ii, :), ptsPolygon(index-1, :), ptsPolygon(index, :),...
            ptsPolygon(index+1, :));
    end
end
end


% contour error for the intermediate points. determine the region firstly. Refer to Erkorkmaz,
% 2006, IJMTM, virtual CNC system. part II
function con = contourErrorByRegion(pt, pt0, pt1, pt2)
% 3-D points
vec1 = pt1 - pt0;
vec2 = pt2 - pt1;
vec1 = vec1/norm(vec1);
vec2 = vec2/norm(vec2);
vec3 = vec1 + vec2;
vec = pt - pt1;
N_i = dot(vec1, vec);
N_i_1 = dot(vec2, vec);
B_i = dot(vec3, vec);
if norm(vec3) ~= 0
    if (B_i <0) && (N_i <= 0) flag = 1; end
    if(B_i >= 0) && (N_i_1 >= 0) flag = 2; end
    if (N_i >0) && (N_i_1 < 0) flag = 3; end
else
    if (N_i > 0) && (N_i_1 < 0)
        flag = 3;
    else
        flag = 2;
    end
end
switch flag
    case 1
        con = distancePoint2Line(pt, pt0, pt1);
    case 2
        con = distancePoint2Line(pt, pt1, pt2);
    case 3
        con = norm(vec); % distance between pt and pt1.
end
end


function dis = distancePoint2Line(pt, pt0, pt1)
% calculate the distance from pt to a line start from pt0 to pt1.
% only applied to 3-D situation.
vec1 = pt - pt0;
vec2 = pt1 - pt0;
dis = norm(cross(vec1, vec2)) / norm(vec2);
end

