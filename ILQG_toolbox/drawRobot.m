function drawRobot(pos,offset_th,L)

Pf= [pos(2) pos(3)]';
thetaf	= pos(1);
XYf	= [ 0 L nan 0  0 nan   -0.5*L  0.5*L nan -0.5*L -0.5*L nan -0.5*L  0.5*L nan 0.5*L L nan L 0.5*L nan -0.3*L 0.3*L nan -0.3*L 0.3*L; 
        0 0 nan -0.7*L 0.7*L nan 0.5*L 0.5*L nan -0.5*L 0.5*L nan -0.5*L -0.5*L nan -0.5*L 0 nan 0 0.5*L nan -0.7*L -0.7*L nan 0.7*L 0.7*L];
RMatf	= [cos(thetaf) -sin(thetaf); sin(thetaf) cos(thetaf)];
lidar = [cosd(offset_th) -sind(offset_th); sind(offset_th) cosd(offset_th)]*[0 2*L;...
                                                                             0 0];

for i=1:26
  XYf(:,i) = RMatf * XYf(:,i) + Pf;
end
for i = 1:2
    lidar(:,i) = RMatf * lidar(:,i) + Pf;
end

plot(XYf(1,:),XYf(2,:),'b-','Linewidth',2); hold on;
plot(lidar(1,:),lidar(2,:),'r-','Linewidth',2.3);

end