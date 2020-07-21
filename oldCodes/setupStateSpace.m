% setupStateSpace
myNetwork=0; % assign dummy value

% A, B, and F matrices are (2n)x(2n)
n=8; % n=length(feeder.edges) for example, set as the number of nodes
A=eye(2*n) %(2x)x(2n) dims
IndicMat=zeros(2*n,2*n) %(2x)x(2n) dims
[R,X]=createRXmatrices(myNetwork);
B=[X R; (-1/2)*R (1/2)*X]; % (2x)x(2n) dims