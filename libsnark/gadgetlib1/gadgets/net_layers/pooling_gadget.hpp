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

}// libsnark