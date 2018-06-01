# homomorphic_LBP
Computation of Local Binary Pattern (LBP) codes of encrypted pixels, using HElib.

This project tries to implement a face recognition method for encrypted images. The face recognition method is based on Local 
Binary Patterns (LBP), more about LBP can be found here http://www.scholarpedia.org/article/Local_Binary_Patterns/. The project
uses HElib (homomorphic encryption library) to deal with processing of encrypted data, more about HElib can be found here 
https://github.com/shaih/HElib.

So far, the project runs very slowly because the techniques used were the raw adaptation of the LBP methods to encrypted data,
without deep modifications of the underlying HE scheme. 

Happy coding!
