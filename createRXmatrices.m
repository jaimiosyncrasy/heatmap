function [R,X] = createRX(myNetwork)
    % paper ref: https://scholarspace.manoa.hawaii.edu/bitstream/10125/64132/0316.pdf

    % Input: feeder object
    % output: R matrix, X matrix
    % dims: R and X are each nxn where n=number of nodes in feeder including slack bus

    % Python code seems to access impedances with:
    %  impedance = feeder.network.get_edge_data(pred_list_max[0], max_depth_bus, default=None)['connector']
    % impedance.Z

    n=8; % for our example
    clear j; % reserve for imaginary symbol
    % MATLAB indices start at 1 so node 0 to 7 is represented with matrix
    % indices 1 to 8
    z(1,2)=1+j; z(2,1)=z(1,2);  % line 01
    z(2,3)=1+j; z(3,2)=z(2,3); % line 12
    z(3,4)=1+j; z(4,3)=z(3,4); % line 23
    z(4,5)=1+j; z(5,4)=z(4,5); % line 34
    z(3,6)=1+j; z(6,3)=z(3,6); % line 34
    z(6,7)=1+j; z(7,6)=z(6,7); % line 56
    z(6,8)=1+j; z(8,6)=z(6,8); % line 57

    %% Compute P{i}
    % Pi is the line path from node i to substation node,node0
    % each row is a tuple, representing an edge

    P=cell(n,1); % each cell holds matrix representing a set of lines

    %Modify func get_total_impedance_between_two_buses(feeder, node_name_1, node_name_2)
    % to return not only the impedance but the set of edges on the path between
    % node1 and node2
    % Then can simply setup P{i} in the following way:
    % for i = 1:n
    %     P{i}=get_total_impedance_between_two_buses(feeder,0,i)
    % end

    %% Example of R23 and X23, elements of R and X matrices
    P{3}=[4 3; 3 2; 2 1]; % P3 is set of lines 32,21,10
    P{2}=[3 2; 2 1]; % P2 is set of lines 21,10
    intersection = intersect(P{2},P{3}, 'rows') % returns the rows common to P2 and P3
    Rsum=0; Xsum=0;
    for i = 1:size(intersection,1) % iterate down the rows
        Rsum=Rsum + real(z(intersection(i,1), intersection(i,2)));
        Xsum=Xsum + imag(z(intersection(i,1), intersection(i,2)));
    end
    R23=2*Rsum
    X23=2*Xsum


    %% Generalizing above example to populate all of R and X matrices

    P{2}=[2 1]; % P1 is set of lines 10
    P{3}=[3 2; 2 1]; % P2 is set of lines 21,10
    P{4}=[4 3; 3 2; 2 1]; % P3 is set of lines 32,21,10
    P{5}=[5 4; 4 3; 3 2; 2 1]; % P4 is set of lines 43,32,21,10
    P{6}=[6 3; 3 2; 2 1]; % P5 is set of lines 52,21,10
    P{7}=[7 6; 6 3; 3 2; 2 1]; % P6 is set of lines 65,52,21,10
    P{8}=[8 6; 6 3; 3 2; 2 1]; % P7 is set of lines 75,52,21,10

    for r=2:n
        for s=2:n
            intersection = intersect(P{r},P{s}, 'rows'); % returns the rows common to P2 and P3
            Rsum=0; Xsum=0;
            for k = 1:size(intersection,1) % iterate down the rows
                Rsum=Rsum + real(z(intersection(k,1), intersection(k,2)));
                Xsum=Xsum + imag(z(intersection(k,1), intersection(k,2)));
            end
            R(r,s)=2*Rsum; R(s,r)=R(r,s);
            X(r,s)=2*Xsum; X(s,r)=X(r,s);
        end
    end

    % Print results
    fprintf('\n')
    disp('R='); disp(R);
    disp('X='); disp(X);
end

