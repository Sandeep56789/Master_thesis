function position = eval_receiver_position(t,posini,velocity)

%
% calculate the position of the receiver along the time t according to the given
% initial position and velocity.
% the mouvement is assumed rectiligne and uniform
%
% input :   - t = vector of time at which the positions are evaluated
%           - posini = vector of initial position in (x,y,z)
%           - velocity = vector of velocity in (x,y,z)
%
%
% output :  - position = matrix 3*Nt of position in (x,y,z)
%

% uniform and rectilinear motion %
%===============================%
position(1,:) = posini(1) + velocity(1)*t;
position(2,:) = posini(2) + velocity(2)*t;
position(3,:) = posini(3) + velocity(3)*t;
