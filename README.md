# NodeEmbedding

A graph representation learning algorithm extracts the node embeddings of a complex network from the information cascades. Developed from scratch by C++ and OpenMPI. (The MPI version will be released soon)

The input file starts with `NodeID,NodeName` per line, followed by an empty line and then each cascade per line
```
CascadeID;NodeID1,Timestamp1,NodeID2,Timestamp2
```
where `NodeID` and `CascadeID` need to be integers. For example, a valid input file looks like
```
1,1
2,2
3,3

1;1,0.1,2,0.2,3,0.3
```

Compile the code with OpenMPI support and execute:
```
mpirun -n=2 ./sample input=<path to input file> d=10 max_iterations=20
```

Our paper ["Scalable Prediction of Global Online Media News Virality"](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8438324) is published at IEEE Transactions on Computational Social Systems.
```
@article{lu2018scalable,
  title={Scalable Prediction of Global Online Media News Virality},
  author={Lu, Xiaoyan and Szymanski, Boleslaw K},
  journal={IEEE Transactions on Computational Social Systems},
  number={99},
  pages={1--13},
  year={2018},
  publisher={IEEE}
}
```
