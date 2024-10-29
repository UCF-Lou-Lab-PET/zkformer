#include <iostream>

#include <libsnark/gadgetlib1/gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/basic_gadgets.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include "libff/algebra/fields/field_utils.hpp"


/**
 * y = max{x, 0}
 * equals to
 * y = x*(!(x<0))
 * which is
 *      if x>=0: y=x*1=x
 *      if x<0: y=x*0=0
 */
// ReLU max(0,x) comparison_gadget

namespace libsnark {

template <typename FieldT>
class ReLU_gadget : public gadget<FieldT>
{
private:
    // comparison_gadget<FieldT> comparison;
    // pb_variable<FieldT> x;
    // pb_variable<FieldT> threshold;
    pb_variable<FieldT> less_or_eq;
    std::shared_ptr<comparison_gadget<FieldT>> comparison;

public:
    pb_variable<FieldT> out;
    pb_variable<FieldT> in;
    pb_variable<FieldT> threshold;
    pb_variable<FieldT> less;
    ReLU_gadget(){};
    ReLU_gadget(protoboard<FieldT> &pb,
                const pb_variable<FieldT> &out,
                const pb_variable<FieldT> &in,
                const pb_variable<FieldT> &threshold,
                const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold)
    {
        less.allocate(pb, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, FMT(this->annotation_prefix, "less_or_eq"));

    // what is 10? why set n to 10???
    // quantize bit  == 10, the bit size of AB
    // test (in < threshold // in <= threshold), result in less / less_or_eq
    // n ==10 ,n==32
        comparison.reset(new comparison_gadget<FieldT>(pb, 10, in, threshold, less, less_or_eq,
                                                       FMT(this->annotation_prefix, "comparison")));
    };

    // varaibles in the pb: [one, out, in, threshold, sym_1...] 
    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void ReLU_gadget<FieldT>::generate_r1cs_constraints()
{
    comparison->generate_r1cs_constraints();

    this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in, FieldT::one() - less, out),
                                 FMT(this->annotation_prefix, "y = x * (! (x<0) )"));
}
/*
compute the value of intermediate variables based on given in/out
*/
template <typename FieldT>
void ReLU_gadget<FieldT>::generate_r1cs_witness()
{
    comparison->generate_r1cs_witness();
}


/**
 * relu gadget for arrays/vectors
 * relu(pb, OUT, IN, threshold, "out ? max{in, 0}");
 * for the relu computation, we only care about whether x<0 or not, less_or_eq is not needed
*/
template <typename FieldT>
class ReLU_vector_gadget : public gadget<FieldT>
{
private:
    pb_variable_array<FieldT> less_or_eq;
    std::vector<comparison_gadget<FieldT>> comparison;

public:
    size_t n;
    pb_variable_array<FieldT> out;
    pb_variable_array<FieldT> in;
    pb_variable<FieldT> threshold;
    pb_variable_array<FieldT> less;

    ReLU_vector_gadget(protoboard<FieldT> &pb,
                        pb_variable_array<FieldT> &out,
                        pb_variable_array<FieldT> &in,
                        const pb_variable<FieldT> &threshold,
                        size_t n,
                        const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold), n(n)
    {
        less.allocate(pb, n, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, n, FMT(this->annotation_prefix, "less_or_eq"));
        for (size_t i=0; i < n; i++)
        {
            comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], threshold, less[i], less_or_eq[i], "comparison"));
        }
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void ReLU_vector_gadget<FieldT>::generate_r1cs_constraints()
{

    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_constraints();
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], out[i]),
                                    FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) )"));
    }
}

template <typename FieldT>
void ReLU_vector_gadget<FieldT>::generate_r1cs_witness()
{
    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_witness();
    }
}




/**
 * create a look up table in relu gadget
 * st. all comparison gadget share the same space for alpha
*/
// template <typename FieldT>
// class ReLU_vector_gadget_opt : public gadget<FieldT>
// {
// private:
//     pb_variable_array<FieldT> less_or_eq;
//     pb_variable_array<FieldT> alpha_table;
//     // std::vector<comparison_gadget<FieldT>> comparison;
//     std::vector<comparison_gadget_opt<FieldT>> comparison;

// public:
//     size_t n;
//     pb_variable_array<FieldT> out;
//     pb_variable_array<FieldT> in;
//     pb_variable<FieldT> threshold;
//     pb_variable_array<FieldT> less;

//     ReLU_vector_gadget_opt(protoboard<FieldT> &pb,
//                         pb_variable_array<FieldT> &out,
//                         pb_variable_array<FieldT> &in,
//                         const pb_variable<FieldT> &threshold,
//                         size_t n,
//                         const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold), n(n)
//     {
//         less.allocate(pb, n, FMT(this->annotation_prefix, "less"));
//         less_or_eq.allocate(pb, n, FMT(this->annotation_prefix, "less_or_eq"));
//         alpha_table.allocate(pb, 2*n, FMT(this->annotation_prefix, "alpha look up table"));
//         for (size_t i=0; i < n; i++)
//         {
//             // method1: let comparison gadget creat alpha
//             // comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], threshold, less[i], less_or_eq[i], "comparison"));
//             // method2: use pre-defined alpha
//             comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], threshold, less[i], less_or_eq[i], alpha_table, "comparison"));
//         }
//     };

//     void generate_r1cs_constraints();
//     void generate_r1cs_witness();
// };

// template <typename FieldT>
// void ReLU_vector_gadget_opt<FieldT>::generate_r1cs_constraints()
// {
//     for(size_t i = 0; i < n; i++){
//         comparison[i].generate_r1cs_constraints();
//         this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], out[i]),
//                                     FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) )"));
//     }
// }

// template <typename FieldT>
// void ReLU_vector_gadget_opt<FieldT>::generate_r1cs_witness()
// {
//     for(size_t i = 0; i < 2*n; i+=2){
//         pb.val(alpha_table[i]) = FieldT::zero();
//         pb.val(alpha_table[i+1]) = FieldT::one();
//     }
//     for(size_t i = 0; i < n; i++){
//         comparison[i].generate_r1cs_witness();
//     }
// }






/**
 * ReLU_vector_gadget with more constraints
 * namely, y = x * (1-(x<a)) + a * (x<a)
 *  if x<a, y == a
 *  else y==x
 * where in relu, a(threshold) is set to 0 
*/
template <typename FieldT>
class ReLU_vector_gadget_wmc : public gadget<FieldT>
{
private:
    pb_variable_array<FieldT> less_or_eq;
    std::vector<comparison_gadget<FieldT>> comparison;

public:
    size_t n;
    pb_variable_array<FieldT> out;
    pb_variable_array<FieldT> in;
    pb_variable<FieldT> threshold;
    pb_variable_array<FieldT> sym1;
    pb_variable_array<FieldT> sym2;
    pb_variable_array<FieldT> less;

    ReLU_vector_gadget_wmc(protoboard<FieldT> &pb,
                        pb_variable_array<FieldT> &out,
                        pb_variable_array<FieldT> &in,
                        const pb_variable<FieldT> &threshold,
                        size_t n,
                        const std::string &annotation_prefix = "") : gadget<FieldT>(pb, annotation_prefix), out(out), in(in), threshold(threshold), n(n)
    {
        // intermediate syms
        sym1.allocate(pb, n, FMT(this->annotation_prefix, "x*(!x<a)"));
        sym2.allocate(pb, n, FMT(this->annotation_prefix, "a*(x<a)"));
        less.allocate(pb, n, FMT(this->annotation_prefix, "less"));
        less_or_eq.allocate(pb, n, FMT(this->annotation_prefix, "less_or_eq"));
        for (size_t i=0; i < n; i++)
        {
            comparison.push_back(comparison_gadget<FieldT>(pb, 10, in[i], threshold, less[i], less_or_eq[i], "comparison"));
        }
    };

    void generate_r1cs_constraints();
    void generate_r1cs_witness();
};

template <typename FieldT>
void ReLU_vector_gadget_wmc<FieldT>::generate_r1cs_constraints()
{

    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_constraints();
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], sym1[i]),
                                    FMT(this->annotation_prefix, "x[i] * (! (x[i]<a) )"));
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(threshold, less[i], sym2[i]),
                                    FMT(this->annotation_prefix, "a * (x[i]<0) "));
        this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(sym1[i] + sym2[i], 1, out[i]),
                                    FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) ) + a * (x[i]<0)"));
        // this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(in[i], FieldT::one() - less[i], out[i]),
        //                             FMT(this->annotation_prefix, "y[i] = x[i] * (! (x[i]<0) )"));
    }
}

template <typename FieldT>
void ReLU_vector_gadget_wmc<FieldT>::generate_r1cs_witness()
{
    for(size_t i = 0; i < n; i++){
        comparison[i].generate_r1cs_witness();
        this->pb.val(sym1[i]) = this->pb.val(in[i]) * (FieldT::one() - this->pb.val(less[i]));
        this->pb.val(sym2[i]) = this->pb.val(threshold) * this->pb.val(less[i]);
    }
}



}// libsnark