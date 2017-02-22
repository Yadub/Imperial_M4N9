function [ texoutput ] = tabletop50( top50 )
%TABLETOP50 Produces a latex array out of the top50 array
%   Detailed explanation goes here

nums1 = 1:25;
nums2 = 26:50;
top25 = top50(1:25);
next25 = top50(26:50);

o1 = [nums1',top25];
o2 = [nums2',next25];
o = [o1,o2];

texoutput = latex(sym(o));

end

