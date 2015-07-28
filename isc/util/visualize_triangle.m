function ha= visualize_triangle(fig_id,points,u_start,shoot_dir,cl1,lw1,cl2,lw2,k);
%% recon is the shooting direction, as reconstructed within the
%% previous triangle.
%% We are now extending the path along this direction within the new triangle.

coordinates = points*u_start';
figure(fig_id);
if nargout==1,
    ha= axes;
end
not_in = find(~u_start);
in     = find(u_start);
if length(not_in)>length(in),
    [not_in,in] = swap(not_in,in);
end

hq = quiver3(coordinates(1),coordinates(2),coordinates(3),shoot_dir(1),shoot_dir(2),shoot_dir(3));
hold on,
if ~isempty(k)
text(coordinates(1)+ shoot_dir(1),coordinates(2)+ shoot_dir(2),coordinates(3)+ shoot_dir(3),num2str(k));
end
if ~exist('clr','var'),
    clr = 'k';
end

set(hq,'color',cl1,'linewidth',lw1);
hold on,
c_out = [points(1,not_in);points(2,not_in);points(3,not_in)];
c_in  = [points(1,in);points(2,in);points(3,in)];

pt_1 = c_in(:,1); pt_2 = c_in(:,2);

hold on
plot3([pt_1(1),c_out(1)],[pt_1(2),c_out(2)],[pt_1(3),c_out(3)],'color',cl2,'linewidth',lw2);
hold on
plot3([pt_2(1),c_out(1)],[pt_2(2),c_out(2)],[pt_2(3),c_out(3)],'color',cl2,'linewidth',lw2);
hold on
plot3([pt_2(1),pt_1(1)],[pt_2(2),pt_1(2)],[pt_2(3),pt_1(3)],'color',cl2,'linewidth',lw2);


function [a,b] = swap(a,b);
temp = a;
a= b;
b= temp;