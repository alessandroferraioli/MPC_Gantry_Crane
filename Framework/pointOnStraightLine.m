function result = pointOnStraightLine(x1,y1,x2,y2,xp,yp) 
        line = yp-y1*(y2-y1)*(xp-x1)/(x2-x1);
        if(line<10^-7)
            result = 1;
        else
            result = 0;
        end
        
end

%CHEK WHERE IS THE TARGET RESPECT TO THE RECTANGLE 
if(pointOnStraightLine(x1,y1,x2,y


if(caseConst == 1)

    respect_rectangles = (xTar - x5)*(y5 - y6) - (yTar -y5)*(x5-x6);
    
    if(respect_rectangles >0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    
end


if(caseConst == 2)
respect_rectangles = (xTar - x5)*(y5 - y6) - (yTar -y6)*(x5-x6);
    
    if(respect_rectangles > 0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    
end


if(caseConst == 3)
    %NB we have - cause this shape is obtained by a rotation of pi/2
    respect_rectangles = -((xTar - x5)*(y5 - y6) - (yTar -y5)*(x5-x6));
    disp(respect_rectangles);
    if(respect_rectangles < 0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
  
end

if(caseConst == 4)
    respect_rectangles = -((xTar - x5)*(y5 - y6) - (yTar -y5)*(x5-x6))
    
    if(respect_rectangles <0)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    

    
end

if(caseConst == 5)
    
    if(yTar < y5 && xTar>x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
    
    
end


if(caseConst == 6)
    if(yTar > y5 && xTar<x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end

  
end


if(caseConst == 7)
    
    if(yTar < y5 && xTar<x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end
 
end


if(caseConst == 8)
    
    if(yTar > y5 && xTar>x5)
        %second rect
        
        param.xTarSigned = (x5 + x2)*0.5;
        param.yTarSigned = (y5 + y2)*0.5;
        param.whichRectTarget = 2;
        disp('TARGET IS IN THE SECOND RECT');
    else
        %first rect
        disp('TARGET IS IN THE FIRST RECT');
        param.whichRectTarget = 1;
        param.xTarSigned = xTar;
        param.yTarSigned = yTar;
    end

end





% %-----------------using middle point
% 
% param.xTarSigned = (x5 + x2)*0.5;
% param.yTarSigned = (y5 + y2)*0.5;

%-----------------using intersection 

%m23 = (y2 - y3)/(x2 - x3);

% if(x2 == x5)
%     pointTarSigned = [x2 ; x2*m23-m23*xTar+yTar];
%     disp('vertical');
%     
% else
%     m25 = (y2-y5)/(x2-x5);
%     ATarSigned = [[-m25 1];[-m23 1]];
%     bTarSigned = [-m25*x5 + y5 ; -m23*xTar + yTar];
% 
%     pointTarSigned = ATarSigned\bTarSigned;
% 
% end
% param.xTarSigned = pointTarSigned(1);
% param.yTarSigned = pointTarSigned(2);


if(caseConst == 1 || caseConst  == 2|| caseConst  == 3|| caseConst  == 4)

    D = zeros(6,8);
    
    D(1,1) = -m2;
    D(1,3) = 1;
    
    D(2,1) = -m1;
    D(2,3) = 1;
    
    D(3,1) = -m2;
    D(3,3) = 1;
    D(3,5) = -len*m2;
    D(3,7) = len;
    
    D(4,1) = -m1;
    D(4,3) = 1;
    D(4,5) = -len*m1;
    D(4,7) = len;


end

    
    
D(5,5) = 1;
D(6,7) = 1;
D2 = D;


D_check = zeros(4,2);
D_check(1,1) = -m2;
D_check(1,3) = 1;

D_check(2,1) = -m1;
D_check(2,3) = 1;

D_check(3,1) = -m2;
D_check(3,3) = 1;

D_check(4,1) = -m1;
D_check(4,3) = 1;

check_const1  = [c1Lower;c2Lower;c1Upper;c2Upper]';

if(D_check * [xTar;yTar] < check_const1)
            param.whichRectTarget = 1;
            disp('TARGET IN FIRST RECT');
else
    
            param.whichRectTarget = 2;
            disp('TARGET IN SECOND RECT');
end



angleDegreeConst = deg2rad(1);
angleConstraint = [-angleDegreeConst , angleDegreeConst];


cl1 = [c1Lower;c2Lower;c1Lower;c2Lower;angleConstraint(1);angleConstraint(1)];
ch1 = [c1Upper;c2Upper;c1Upper;c2Upper;angleConstraint(2);angleConstraint(2)];
param.cl1 = cl1;
param.ch1 = ch1;
param.D = D;

%1-2 chart
%3-4 mass
%5-6 angle
cl2 = [c1Lower2;c2Lower2;c1Lower2;c2Lower2;angleConstraint(1);angleConstraint(1)];
ch2 = [c1Upper2;c2Upper2;c1Upper2;c2Upper2;angleConstraint(2);angleConstraint(2)];
param.cl2 = cl2;
param.ch2 = ch2;
param.D2 = D2;




%+++++++++++++++++++++++++++++++++++++++++++matrix 1
% Declare penalty matrices and tune them here:
Q=eye(8);
% Q(1,1) = 3;
% Q(2,2) = 1;
% Q(3,3) = 3;
% Q(2,2) = 1;

Q(1,1) = 10;
Q(3,3) = 10;
Q(5,5) = 5;
Q(7,7) = 5;

weightInput = 0.5;
R=eye(2)*weightInput; 
P=Q; % terminal weight


% First rect
[Dt,Et,bt]=genStageConstraints(A,B,D,cl1,ch1,ul,uh);
[DD,EE,bb]=genTrajectoryConstraints(Dt,Et,bt,N);
[Gamma,Phi] = genPrediction(A,B,N);
[F,J,L]=genConstraintMatrices(DD,EE,Gamma,Phi,N);
[H,G] = genCostMatrices(Gamma,Phi,Q,R,P,N);       
H = chol(H,'lower');
H=(H'\eye(size(H)))';

param.Gamma = Gamma;
param.Phi = Phi;


% second rect
[Dt2,Et2,bt2]=genStageConstraints(A,B,D2,cl2,ch2,ul,uh);
[DD2,EE2,bb2]=genTrajectoryConstraints(Dt2,Et2,bt2,N);
[F2,J2,L2]=genConstraintMatrices(DD2,EE2,Gamma,Phi,N);



param.H = H;
param.G = G;
param.F = F;
param.J = J;
param.L = L;
param.Gamma = Gamma;
param.Phi = Phi;
param.EE = EE;
param.bb = bb;



param.F2 = F2;
param.J2 = J2;
param.L2 = L2;
param.bb2 = bb2;
param.EE2 = EE2;




%BACKUP : Controller LQR close to the solution
R_backup = 0.02*eye(2);
Q_backup = C' * C;

param.K_LQR = dlqr(A,B,Q_backup,R_backup);




end % End of mySetup

%% ##################################################################################################




