
for i=1:N
    fitQ{i} = G.exp(mean_pos,time(i)*mean_vel);
    x_test(:,i) = fitQ{i};   
end    

figure('DefaultAxesFontSize',25,'DefaultAxesFontWeight','bold');
[l1,l2,l3]=sphere(20); 
surf(l1,l2,l3,'FaceColor','w','EdgeColor', [.7 .7 .7], 'FaceAlpha', 0);  % plot the sphere with green color
hold on;
plot3(x_noisy(1,:) , x_noisy(2,:), x_noisy(3,:), 'ks','Linewidth',2)
hold on;
plot3(x_test(1,:) , x_test(2,:), x_test(3,:), '-r','Linewidth',3)
axis equal
grid off
axis off
xlim([-1.2, 1.2])
ylim([-1.2, 1.2])
zlim([-1.2, 1.2])
%title('Fitting on Gr(3,1)')
view(-120,-50)
set(gcf,'color','w'); % Set background color to white
set(gca,'box','off'); % Turn off the axes box
%saveas(gcf,'fit_Gr31.png');