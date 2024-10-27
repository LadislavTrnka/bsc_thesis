% Removing rows and columns of matrix A
% To remove corresponding rows and columns in the considered matrix in stress problems

function R=removing_stress(A, NX, NY, RorC, type)
    U= 1:NX; % position of upper boundary conditions
    R= NX:NX:NX*NY; % position of right boundary conditions
    L= 1:NX:NX*NY; % position of left boundary conditions
    D= NX*(NY-1)+1:NX*NY; % position of down boundary conditions
    if type==1
        k = [L R]; % left and right boundary conditions
    elseif type==2
        k = [U D]; % up and down boundary conditions
    elseif type==3
        k = [L R U D]; % all boundary conditions
    end
    if RorC==1
        A(k,:)=[];
    elseif RorC==0
        A(:,k)=[];
    end
    R=A;
end
