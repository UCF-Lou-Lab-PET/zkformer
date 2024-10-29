/** @file
 *****************************************************************************
 * @author     This file is part of libsnark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SIMPLE_EXAMPLE_TCC_
#define SIMPLE_EXAMPLE_TCC_

#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libsnark/gadgetlib1/gadgets/net_layers/relu_gadget.hpp>

namespace libsnark {

/* NOTE: all examples here actually generate one constraint less to account for soundness constraint in QAP */

template<typename FieldT>
r1cs_example<FieldT> gen_r1cs_example_from_protoboard(const size_t num_constraints)
{
    const size_t new_num_constraints = num_constraints - 1;

    /* construct dummy example: inner products of two vectors */
    protoboard<FieldT> pb;
    pb_variable_array<FieldT> A;
    pb_variable_array<FieldT> B;
    pb_variable<FieldT> res;

    // the variables on the protoboard are (ONE (constant 1 term), res, A[0], ..., A[num_constraints-1], B[0], ..., B[num_constraints-1])
    // followed by S[0] S[1] ... S[new_num_cs - 1] (which is declared by gen_r1cs_constraints() and assigned by gen_r1cs_witness() )
    res.allocate(pb, "res");
    A.allocate(pb, new_num_constraints, "A");
    B.allocate(pb, new_num_constraints, "B");

    inner_product_gadget<FieldT> compute_inner_product(pb, A, B, res, "compute_inner_product");
    compute_inner_product.generate_r1cs_constraints();

    /* fill in random example */
    for (size_t i = 0; i < new_num_constraints; ++i)
    {
        // pb.val(A[i]) = FieldT::random_element();
        // pb.val(B[i]) = FieldT::random_element();
        pb.val(A[i]) = FieldT(i+1);
        pb.val(B[i]) = FieldT(i+1);
    }

    size_t primary_input_size =  2 * new_num_constraints;
    pb.set_input_sizes(0);

    compute_inner_product.generate_r1cs_witness();
    // variable number is n; but there is a one in index0, so n+1 terms
    // check the full_assigments
    // i.e. 30= (1,2,3,4) * (1,2,3,4)T [one, out, a0~a3, b0~b3, product0~2]
    for(size_t i=0; i < pb.num_variables() + 1; i++){
        std::cout<<pb.val(i).as_ulong() <<" ";
    }
    std::cout<<std::endl;
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    // assert(cs.is_satisfied(pb.primary_input(),pb.auxiliary_input()));
    return r1cs_example<FieldT>(cs, pb.primary_input(), pb.auxiliary_input());
}

template<typename FieldT>
r1cs_example<FieldT> gen_r1cs_example_from_relu(const size_t num_constraints)
{
    protoboard<FieldT> pb;
    pb_variable<FieldT> out;
    pb_variable<FieldT> in;
    pb_variable<FieldT> threshold;

    out.allocate(pb, "out");
    in.allocate(pb, "in");
    threshold.allocate(pb, "threshold");
    

    // trying to prove out = max{in, 0}, where the threshold is 0
    ReLU_gadget<FieldT> relu(pb, out, in, threshold, "out ? max{in, 0}");
    relu.generate_r1cs_constraints();

    /* fill in random example */
    // assign wavlues for variables
    pb.val(out) = FieldT(3);
    pb.val(in) = FieldT(3);
    // pb.val(threshold) = FieldT::zero();
    pb.val(threshold) = FieldT(0);


    // size_t primary_input_size =  2 * new_num_constraints;
    // pb.set_input_sizes(0);

    relu.generate_r1cs_witness();
    // variable number is n; but there is a one in index0, so n+1 terms
    // check the full_assigments
    // i.e. 30= (1,2,3,4) * (1,2,3,4)T [one, out, a0~a3, b0~b3, product0~2]
    // for(size_t i=0; i < pb.num_variables() + 1; i++){
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    for(size_t i=0; i < pb.num_variables() + 1; i++){
        if(i == pb.num_variables()){
            std::cout<< pb.val(i).inverse().as_ulong() <<" ";
            break;
        }
        std::cout<<pb.val(i).as_ulong() <<" ";
    }
    std::cout<<std::endl;
    // pb.dump_variables();
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    assert(cs.is_satisfied(pb.primary_input(),pb.auxiliary_input()));
    return r1cs_example<FieldT>(cs, pb.primary_input(), pb.auxiliary_input());
}

template<typename FieldT>
r1cs_example<FieldT> gen_r1cs_example_from_relu_vector(const size_t vector_size)
{
    libff::enter_block("Call to generate_r1cs_example_from_relu_vector");
    protoboard<FieldT> pb;
    pb_variable_array<FieldT> out;
    pb_variable_array<FieldT> in;
    pb_variable<FieldT> threshold;

    out.allocate(pb, vector_size, "out");
    in.allocate(pb, vector_size, "in");
    threshold.allocate(pb, "threshold");
    
    libff::enter_block("set constraints");
    // trying to prove out = max{in, 0}, where the threshold is 0
    ReLU_vector_gadget<FieldT> relu(pb, out, in, threshold, vector_size, "out ? max{in, 0}");
    // try with pvcnn's code
    // ReLU_vector_gadget_wmc<FieldT> relu(pb, out, in, threshold, vector_size, "out ? max{in, 0}");
    relu.generate_r1cs_constraints();
    libff::leave_block("set constraints");

    libff::enter_block("set variables");
    // assign wavlues for variables

    // value can't be too large, ~8 bit
    std::vector<int> temp_vector;
    for(size_t i=0; i<vector_size; i++){
        temp_vector.push_back((std::rand() % 1000) - 300);
        // temp_vector.push_back(i+1);
        // std::cout<< temp_vector[i] <<" ";
    }



    for(size_t i=0; i<vector_size; i++){
        pb.val(out[i]) =  FieldT (temp_vector[i] > 0 ? temp_vector[i] : 0);
    }
    
    // set a wrong assignment 
    // pb.val(out[2]) = FieldT(3);

    for(size_t i=0; i<vector_size; i++){
        pb.val(in[i]) = FieldT(temp_vector[i]);
    }
    // pb.val(threshold) = FieldT::zero();
    pb.val(threshold) = FieldT(0);


    // size_t primary_input_size =  2 * new_num_constraints;
    // pb.set_input_sizes(0);

    relu.generate_r1cs_witness();
    libff::leave_block("set variables");
    // // variable number is n; but there is a one in index0, so n+1 terms
    // // check the full_assigments

    // size_t i=0;

    // std::cout<<"one: "<<pb.val(i++).as_ulong()<<std::endl;

    // std::cout<<"output: ";
    // for(; i < vector_size + 1; i++){
    //     std::cout<<"x_" << i <<": ";
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"ipput: ";
    // for(; i < 2 * vector_size + 1; i++){
    //     std::cout<<"x_" << i <<": ";
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    // std::cout<<std::endl;
    
    // std::cout<<"theshold: "<<pb.val(i++).as_ulong()<<std::endl;
    
    // std::cout<<"sym1: ";
    // for(; i < 3 * vector_size + 2; i++){
    //     std::cout<<"x_" << i <<": ";
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"sym2: ";
    // for(; i < 4 * vector_size + 2; i++){
    //     std::cout<<"x_" << i <<": ";
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"less: ";
    // for(; i < 5 * vector_size + 2; i++){
    //     std::cout<<"x_" << i <<": ";
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"less_eq: ";
    // for(; i < 6 * vector_size + 2; i++){
    //     std::cout<<"x_" << i <<": ";
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    // std::cout<<std::endl;

    // std::cout<<"in the compare gadget: "<<std::endl;
    // size_t number_of_out = i;
    // size_t variable_size_in_bit = 10;
    // for(; i < pb.num_variables() + 1; i++){
    //     std::cout<<"x_" << i <<": ";
    //     if( (i + 1- number_of_out) % (variable_size_in_bit+3) ==0 ){
    //         if(pb.val(i).as_ulong() == 0){
    //             std::cout<<pb.val(i).as_ulong();
    //         }else{
    //             std::cout<<pb.val(i).inverse().as_ulong();
    //         }
    //         std::cout<<std::endl;
    //         // std::cout<<std::endl;
    //         continue;
    //     }
    //     std::cout<<pb.val(i).as_ulong() <<" ";
    // }
    // std::cout<<std::endl;
    // // end of printing assignments

    // pb.dump_variables();
    r1cs_constraint_system<FieldT> cs(pb.get_constraint_system());
    // this is sanity check for the generator
    assert(cs.is_satisfied(pb.primary_input(),pb.auxiliary_input()));

    libff::leave_block("Call to generate_r1cs_example_from_relu_vector");
    return r1cs_example<FieldT>(cs, pb.primary_input(), pb.auxiliary_input());
}




} // libsnark
#endif // R1CS_EXAMPLES_TCC_
