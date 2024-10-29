# zkFormer: Verifiable and Private Transformer Inference

This is the code for zkFormer: Verifiable and Private Transformer Inference. This code is built upon [libsnark](https://github.com/scipr-lab/libsnark): a C++ library for zkSNARK proofs.

## Overview and installation

For better understanding of the code structure and functions, the idea of zkSNARKs is briefly introduced. A prover who knows the witness for the NP statement (i.e., a satisfying input/assignment) can produce a short proof attesting to the truth of the NP statement. This proof can be verified by anyone, and offers the following properties.

-   __Zero knowledge:__
    the verifier learns nothing from the proof beside the truth of the statement (i.e., the value _qux_, in the above examples, remains secret).
-   __Succinctness:__
    the proof is short and easy to verify.
-   __Non-interactivity:__
    the proof is a string (i.e. it does not require back-and-forth interaction between the prover and the verifier).
-   __Soundness:__
    the proof is computationally sound (i.e., it is infeasible to fake a proof of a false NP statement). Such a proof system is also called an _argument_.
-   __Proof of knowledge:__
    the proof attests not just that the NP statement is true, but also that the
    prover knows why (e.g., knows a valid _qux_).

We refer to [libsnark](https://github.com/scipr-lab/libsnark) for more details. Please install libsnark properly.

## Testing

To test the CRPC and PSQ in zkFormer, please run tests for ```libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/tests/test_r1cs_gg_ppzksnark_conv.cpp```. To test a matrix multiplication (Figure 7 as well as Table 3), you can specify the shape of the matrix multiplication. Then,under the ```./build``` directory, you can build and run the test by:

    make && make check

To check the runtime, memory usage and breakdown, check the log files in ```build/Testing/Temporary/LastTest.log```.