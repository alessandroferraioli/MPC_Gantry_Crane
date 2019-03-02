function [ shape ] = generateShape( )
%GENERATESHAPE Create the shape to test with

shape.c      = [0.00, 0.05;
                0.25, 0.30;
                0.50, 0.05;
                0.45, 0.00;
                0.25, 0.20;
                0.05, 0.00];
            
            
shape.eps_r  = 0.02;
shape.eps_t  = 0.02;
shape.start  = [0.1, 0.1];
shape.start  = [0.05, 0.05];

shape.target = [0.45, 0.05];

end

