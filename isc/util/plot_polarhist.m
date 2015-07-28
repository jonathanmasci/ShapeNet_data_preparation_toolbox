function plot_polarhist(HIST,rr,th,offset,c)

rr = [0 rr];
th = [0 th];
cc = HIST;

hold on
for r = 1:size(HIST,1)
    for t = 1:size(HIST,2)
        [X,Y] = drawpatch(th(t)+offset,th(t+1)+offset,rr(r),rr(r+1));
        patch(X,Y,cc(r,t),'LineWidth',3);
    end
end
if nargin==5
    colormap(c)
end

function [X,Y] = drawpatch(th1,th2,r1,r2)

X = [r1*cos(th1:(th2-th1)/10:th2), r2*cos(th2:-(th2-th1)/10:th1)];
Y = [r1*sin(th1:(th2-th1)/10:th2), r2*sin(th2:-(th2-th1)/10:th1)];
