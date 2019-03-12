function [ shape ] = generateShape( )
%GENERATESHAPE Create the shape to test with

%default  1
shape.c      = [0.00, 0.05;
                0.25, 0.30;
                0.50, 0.05;
                0.45, 0.00;
                0.25, 0.20;
                0.05, 0.00];
         

shape.start  = [0.05, 0.05];
shape.target = [0.45, 0.05];




theta =0 ; %rotation from default one
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


% 
%  % UP RIGHT
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
% 
%       
            

% 
% % BOTTOM RIGHT
% shape.c      = [0.40, 0.40;
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
% 
% % BOTTOM LEFT
% shape.c      = [0.40, 0.05;
%                 0.05, 0.05;
%                 0.05, 0.4;
%                 0.12, 0.4;
%                 0.12, 0.12;
%                 0.4, 0.12];
%          
% 
% shape.start  = [0.38, 0.07];
% shape.target = [0.07, 0.38];
% 
% 
% % UP LEFT
% shape.c      = [ 0.05, 0.05;
%                 0.05, 0.4;
%                 0.4, 0.4;
%                 0.4, 0.33;
%                 0.12, 0.33;
%                 0.12, 0.05];
%          
% 
% shape.start  = [0.07, 0.07];
% shape.target = [0.37, 0.37];

            
shape.eps_r  = 0.02;
shape.eps_t  = 0.02;



            
            




end

