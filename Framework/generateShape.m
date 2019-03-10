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
            
            
            
%upside down    2                
% shape.c      = [0.6, 0.25;
%                 0.35, 0.00;
%                 0.10, 0.25;
%                 0.15, 0.30;
%                 0.35, 0.10;
%                 0.55, 0.30];
%             
% shape.start  = [0.15, 0.25];
% shape.target = [0.55, 0.25];



theta =-pi/4; %rotation from default one
R = [[cos(theta) , -sin(theta)];[sin(theta),cos(theta)]]';
for i=1:1:6
    
   oldPoint = [shape.c(i,1),shape.c(i,2)]';
   newPoint = R*oldPoint;
   shape.c(i,1) = newPoint(1);
   shape.c(i,2) = newPoint(2);
   
end
shape.start  = R*shape.start';
shape.start = shape.start';
shape.target = R*shape.target';
shape.target    = shape.target';     



            
shape.eps_r  = 0.02;
shape.eps_t  = 0.02;





end

