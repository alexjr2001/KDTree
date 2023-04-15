# KDTree
Our C++ implementation of the KD Tree has a template for managing the data type and the number of dimensions used to represent our points. This KD Tree follows what is described in the following image:

<p align = "center"><img width="326" alt="image" src="https://user-images.githubusercontent.com/63054183/232251319-8af75596-0e74-440e-89ba-8764096c29ff.png"></p

Our KD Tree has a lot of added functions such as rangeQueries, it can look for the closest neighbour or K-closest neighbours. We can edit our initial tree with insertions, deletions and updates; all of them supported by a search function. 

### Libraries used  
```
#include<iostream>
#include<algorithm> 
#include<vector>
#include<cmath>
#include<math.h>
```

### Running the project
Firstly, we compile it in the path where is stored [myKD.cpp]():
```
$ g++ myKD.cpp -o kd.exe 
```
Secondly, we run it:
```
$ ./kd
```


Let's see an example:
1. We enter the number of neighbours to find (4)
2. We write the coordinates of the point for looking for its KNN (2 10)
3. It prints the result
4. Insert the coordinates of the node to delete (1 10)
5. Insert the new coordinates (20 2) and then it'll print all of them

<p align = "center"><img width="326" alt="image" src="https://user-images.githubusercontent.com/63054183/232250898-497a6892-88b3-4eaa-9d15-fb4b9a479f8f.png"></p>
