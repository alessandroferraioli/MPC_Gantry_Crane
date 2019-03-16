function [ shape ] = generateShape( )
%%GENERATESHAPE Create the shape to test with

            
shape.eps_r  = 0.005;
shape.eps_t  = 0.005;

            
% shape.eps_r  = 0.02;
% shape.eps_t  = 0.02;



%default  1
shape.c      = [0.00, 0.05;
                0.25, 0.30;
                0.50, 0.05;
                0.45, 0.00;
                0.25, 0.20;
                0.05, 0.00];
         

shape.start  = [0.05, 0.05];
shape.target = [0.45, 0.05];

%first_rect , before middle
%shape.target= [0.2 , 0.2];
%first_rect , after middle
%shape.target= [0.28 , 0.26];


%Calculate the other constraints with rotational matrix

%tested angle  0 , +-pi/4 , +-pi/2 , +- pi/6 , +-pi , +- (pi/2+pi/4) , +-
%(pi/2 + pi/6)

theta = 0;
theta = -theta;%rotation from default one
R = [[cos(theta) , -sin(theta)];[sin(theta),cos(theta)]]';
for i=1:1:6
    
   oldPoint = [shape.c(i,1),shape.c(i,2)]';
   newPoint = R*oldPoint;
   shape.c(i,1) = newPoint(1);
   shape.c(i,2) = newPoint(2);%shift in order to don't have negative 
   
end

shiftY = 0;
for i=1:1:6
   if(abs(shape.c(i,2)) > shiftY && shape.c(i,2)<0)
    shiftY = abs(shape.c(i,2));
   end
end
  
shiftX = 0;
for i=1:1:6
   if(abs(shape.c(i,1)) > shiftX && shape.c(i,1)<0)
    shiftX = abs(shape.c(i,1));
   end
end  


shape.c(:,2) = shape.c(:,2)+ones(6,1)*shiftY;
shape.c(:,1) = shape.c(:,1)+ones(6,1)*shiftX;

shape.start  = R*shape.start';
shape.start = shape.start' + [shiftX ,shiftY];

shape.target = R*shape.target';
shape.target    = shape.target'+ [shiftX ,shiftY];    



%% UP RIGHT
% shape.c      = [0.05, 0.4;
%                 0.40, 0.40;
%                 0.4, 0.05;
%                 0.33, 0.05;
%                 0.33, 0.32;
%                 0.05, 0.32];
%          
% 
% shape.start  = [0.07, 0.38];
% shape.target = [0.37, 0.1];
% % %first rect before middle
% %shape.target = [0.3, 0.35];
% % %first rect after middle
% %shape.target = [0.38, 0.35];

 
            

% 
%% BOTTOM RIGHT
%  shape.c      = [0.40, 0.40;
%                 0.40, 0.05;
%                 0.05, 0.05;
%                 0.05, 0.12;
%                 0.33, 0.12;
%                 0.33, 0.40];
%          
% 
% shape.start  = [0.35, 0.38];
% shape.target = [0.07, 0.1];
% 
% % %first rect before middle
% %shape.target = [0.35, 0.18];
% % %first rect after middle
% shape.target = [0.38, 0.08];


%             
% 
%% % BOTTOM LEFT
% shape.c      = [0.40, 0.05;
%                 0.05, 0.05;
%                 0.05, 0.4;
%                 0.12, 0.4;
%                 0.12, 0.12;
%                 0.4, 0.12];
%            %default
% shape.start  = [0.38, 0.07];
% shape.target = [0.07, 0.38];
% % %first rect before middle
%  shape.target = [0.2, 0.11];
% % %first rect after middle
% shape.target = [0.10, 0.11];


% 
% 
%% UP LEFT
% shape.c      = [ 0.05, 0.05;
%                 0.05, 0.4;
%                 0.4, 0.4;
%                 0.4, 0.33;
%                 0.12, 0.33;
%                 0.12, 0.05];
%          
% 
% shape.start  = [0.07, 0.07];
% %second rect
% shape.target = [0.37, 0.37];
% % %first rect before middle
% shape.target = [0.09,0.30]
% % %first rect after middle
% %shape.target = [0.09,0.37]
% 
% 
% 
% shape.c      = [0.10, 0.00;
%                 0.0, 0.2;
%                 0.3, 0.35;
%                 0.35, 0.25;
%                 0.15, 0.15;
%                 0.2, 0.05];
% shape.eps_r  = 0.1*0.02;
% shape.eps_t  = 0.1*0.02;
% shape.start  = [0.2, 0.2];
% shape.target = [0.15, 0.05];
%             
%             
% 
% 
% 
% 
% end

% 
%% different start target pos
% shape.c      = [ 0.03, 0.00;
%                 0.00, 0.06;
%                 0.06, 0.09;
%                 0.07, 0.07;
%                 0.03, 0.05;
%                 0.05, 0.01];
%        %  1-2,1-4,1-6,3-2,3-4,3-6,5-2,5-4,5-6
% points = [ 0.03, 0.01;
%             0.04, 0.02;
%             0.032, 0.052;
%             0.06, 0.08;
%             0.01, 0.06;
%             0.02, 0.05];
% shape.start  = points(1,:);
% shape.target = points(3,:);
