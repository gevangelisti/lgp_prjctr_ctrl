cosc_mat = zeros(size(x(:,1:Ncc)));
ind_neq_0 = (x(:,1:Ncc)~=0);
q_mat = x(:,1:Ncc);
cosc_mat(ind_neq_0) = (1 - cos(q_mat(ind_neq_0))) ./ q_mat(ind_neq_0);
xy_pcc = zeros(2,Ncc,length(t));
for i = 1:length(t)
    xy_pcc(:,:,i) = xy_func(cosc_mat(i,:)',q_mat(i,:)',sinc(q_mat(i,:)'./pi));
end

%%
figure
hold on
xy_pcc_tp = zeros(2,Ncc,length(t_p_vec));
q_mat_tp = zeros(length(t_p_vec),Ncc);
for i = 1:Ncc
    xy_pcc_tp(1,i,:) = interp1(t,squeeze(xy_pcc(1,i,:)),t_p_vec);
    xy_pcc_tp(2,i,:) = interp1(t,squeeze(xy_pcc(2,i,:)),t_p_vec);
    q_mat_tp(:,i) = interp1(t,q_mat(:,i),t_p_vec);
end

handles = [];
for j = 1:length(t_p_vec)
    x_vec = [0,xy_pcc_tp(1,:,j)];
    y_vec = [0,xy_pcc_tp(2,:,j)];
    plot(lxy_fem(1,:,j),lxy_fem(2,:,j),':','LineWidth',mlw,'Color',c(j,:))
    plot(x_vec,y_vec,'--','Color',c(j,:),'LineWidth',mlw)
    q_vec = [0,q_mat_tp(j,:)];
    if j > 1
        for i = 1:Ncc
            r = L(i)/q_mat_tp(j,i);
            phi = linspace(0,q_vec(i+1),1e2);
            x_loc = r*sin(phi);
            y_loc = r*(1-cos(phi));
            han = plot(x_vec(i)+cos(sum(q_vec(1:i)))*x_loc-sin(sum(q_vec(1:i)))*y_loc,y_vec(i)+sin(sum(q_vec(1:i)))*x_loc+cos(sum(q_vec(1:i)))*y_loc,'Color',c(j,:),'LineWidth',mlw);
        end
    else
        han = plot(linspace(0,1,1e2),zeros(1,1e2),'Color',c(j,:),'LineWidth',mlw);
    end
    handles = [handles,han];
end
xlim([0 0.71])
ylim([0 0.71])
legend(handles,['$t = ',num2str(t_p_vec(1)),'$ s'],['$t = ',num2str(t_p_vec(2)),'$ s'],['$t = ',num2str(t_p_vec(3)),'$ s'],['$t = ',num2str(t_p_vec(4)),'$ s'],['$t = ',num2str(t_p_vec(5)),'$ s'], 'Interpreter', 'Latex','FontSize',ax_fs,'Location','Best')
xlabel('x (m)', 'Interpreter', 'Latex','FontSize',ax_fs)
ylabel('y (m)', 'Interpreter', 'Latex','FontSize',ax_fs)
grid
%title(['FEM: step response for a = ',num2str(a),' Nm'])
set(gca,'FontSize',ax_fs)

axes('position',[0.4 .71 .5 .2])
box on % put box around new pair of axes
hold on
for j = 1:length(t_p_vec)
    x_vec = [0,xy_pcc_tp(1,:,j)];
    y_vec = [0,xy_pcc_tp(2,:,j)];
    plot(lxy_fem(1,:,j),lxy_fem(2,:,j),':','LineWidth',mlw,'Color',c(j,:))
    plot(x_vec,y_vec,'--','Color',c(j,:),'LineWidth',mlw)
    q_vec = [0,q_mat_tp(j,:)];
    if j > 1
        for i = 1:Ncc
            r = L(i)/q_mat_tp(j,i);
            phi = linspace(0,q_vec(i+1),1e2);
            x_loc = r*sin(phi);
            y_loc = r*(1-cos(phi));
            han = plot(x_vec(i)+cos(sum(q_vec(1:i)))*x_loc-sin(sum(q_vec(1:i)))*y_loc,y_vec(i)+sin(sum(q_vec(1:i)))*x_loc+cos(sum(q_vec(1:i)))*y_loc,'Color',c(j,:),'LineWidth',mlw);
        end
    else
        han = plot(linspace(0,1,1e2),zeros(1,1e2),'Color',c(j,:),'LineWidth',mlw);
    end
    handles = [handles,han];
end
axis tight
%title(['FEM: step response for a = ',num2str(a),' Nm'])
grid
set(gca,'FontSize',ax_fs)