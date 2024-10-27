% Removing rows and columns of matrix A
% To remove corresponding rows and columns in the considered matrix, \texttt{Matlab} function

function R=removing(A, NX, NY, RorC)
    U= 1:NX; % position of upper boundary conditions
    R= NX:NX:NX*NY; %position of right boundary conditions
    L= 1:NX:NX*NY; % position of left boundary conditions
    D= NX*(NY-1)+1:NX*NY; % position of down boundary conditions
    k = [L R U D];
    if RorC==1 % Removing rows
        A(k,:)=[];
    elseif RorC==0 % Removing columns
        A(:,k)=[];
    end
    R=A;
end
