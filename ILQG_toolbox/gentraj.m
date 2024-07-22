function X = gentraj(init_state,U,W,dt)
% starting from an initial state "init_state" and given a set of
% command inputs "U" and noises "W" -> generate a trajectory with the timestep
% "dt"

% inputs____________________________
% init_state = [theta;x;y] (rad,m,m)
% U = set of command inputs, array 3xN where N is the number of timestep (rad/s;m/s;m/s)
% W = set of noises, array 3xN
% dt = timestep (sec)

% outputs__________________________
% X = set of poses (i.e trajectory) array 3xN (rad;m;m)

X = zeros(3,size(U,2));
X(:,1) = init_state;
for t = 2:size(U,2)
    X(:,t) = f(X(:,t-1),U(:,t-1),W(:,t-1),dt);
end