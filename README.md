# izp-second-project

The program goal is to find object groups with the nearest neighbor method based on their x and y coordinates. 

The input file must have the format of objekty file:
- count=N in the beginning
- N is the number of objects in the input file
- following objects has the following format ```ID X Y```; 
- X and Y are values in between 0 and 1000 including that represents coordinates of an object;
- ID is positive number

The program is launched in the following form ```./cluster FILE [N]```, where N is the final group count. If N is not specified then the final cluster count is 1.
