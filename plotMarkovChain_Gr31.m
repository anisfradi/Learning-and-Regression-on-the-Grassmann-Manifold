


figure('DefaultAxesFontSize',25,'DefaultAxesFontWeight','bold');
subplot(2,2,1)
plot(exp_pos,'b')
%title('Markov chain')
set(gcf,'color','w'); % Set background color to white
set(gca,'box','off'); % Turn off the axes box
%xlabel('iterations')
subplot(2,2,2)
plot(support_pos,pdf_pos,'b')
%title('Density estimation')
set(gcf,'color','w'); % Set background color to white
set(gca,'box','off'); % Turn off the axes box
subplot(2,2,3)
plot(exp_vel,'b')
%title('Markov chain')
set(gcf,'color','w'); % Set background color to white
set(gca,'box','off'); % Turn off the axes box
%title('Density estimation')
%xlabel('iterations')
subplot(2,2,4)
plot(support_vel,pdf_vel,'b')
%title('Density estimation')
set(gcf,'color','w'); % Set background color to white
set(gca,'box','off'); % Turn off the axes box
%saveas(gcf,'fit_MCMC_Gr31.png');


