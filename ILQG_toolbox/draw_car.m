function draw_car(pos)

a = 0.1;
e = 0.02;

x = pos(1);
y = pos(2);
th = pos(3);

P0loc = [0; 0; 1];
P1loc = [a; 0; 1];
P2loc = [0; e/2; 1];
P3loc = [0; -e/2; 1];

P0 = SE2Mat(th,[x;y])*P0loc;
P1 = SE2Mat(th,[x;y])*P1loc;
P2 = SE2Mat(th,[x;y])*P2loc;
P3 = SE2Mat(th,[x;y])*P3loc;

plot([P0(1) P1(1)],[P0(2) P1(2)],'r-o','Linewidth',1.5);
plot([P2(1) P3(1)],[P2(2) P3(2)],'r-o','Linewidth',1.5);

end