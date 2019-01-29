#define BOOST_TEST_MODULE dcgp_mutation_perf
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>

#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>

void perform_mutations_at_once(unsigned in, unsigned out, unsigned rows, unsigned columns, unsigned levels_back,
                               unsigned arity, unsigned N, std::vector<dcgp::kernel<double>> kernel_set,
                               std::string kind = "active")
{
    // Instatiate the expression
    dcgp::expression<double> ex(in, out, rows, columns, levels_back, arity, kernel_set, 123);
    std::cout << "Performing " << N << " " << kind << " mutations, in:" << ex.get_n() << " out:" << ex.get_m()
              << " rows:" << ex.get_r() << " columns:" << ex.get_c() << " levels-back:" << ex.get_l()
              << " active genes: " << ex.get_active_genes().size() << std::endl;
    {
        if (kind == "active") {
            boost::timer::auto_cpu_timer t;
            ex.mutate_active(N);
        } else if (kind == "connection") {
            boost::timer::auto_cpu_timer t;
            ex.mutate_active_cgene(N);
        } else if (kind == "function") {
            boost::timer::auto_cpu_timer t;
            ex.mutate_active_fgene(N);
        } else if (kind == "random") {
            boost::timer::auto_cpu_timer t;
            ex.mutate_random(N);
        }
    }
}

void perform_mutations_one_by_one(unsigned in, unsigned out, unsigned rows, unsigned columns, unsigned levels_back,
                                  unsigned arity, unsigned N, std::vector<dcgp::kernel<double>> kernel_set,
                                  std::string kind = "active")
{
    // Instatiate the expression
    dcgp::expression<double> ex(in, out, rows, columns, levels_back, arity, kernel_set, 123);
    std::cout << "Performing " << N << " " << kind << " mutations, in:" << ex.get_n() << " out:" << ex.get_m()
              << " rows:" << ex.get_r() << " columns:" << ex.get_c() << " levels-back:" << ex.get_l()
              << " active genes: " << ex.get_active_genes().size() << std::endl;
    {
        if (kind == "active") {
            boost::timer::auto_cpu_timer t;
            for (auto i = 0u; i < N; ++i) {
                ex.mutate_active();
            }
        } else if (kind == "connection") {
            boost::timer::auto_cpu_timer t;
            for (auto i = 0u; i < N; ++i) {
                ex.mutate_active_cgene();
            }
        } else if (kind == "function") {
            boost::timer::auto_cpu_timer t;
            for (auto i = 0u; i < N; ++i) {
                ex.mutate_active_fgene();
            }
        } else if (kind == "random") {
            boost::timer::auto_cpu_timer t;
            for (auto i = 0u; i < N; ++i) {
                ex.mutate_random();
            }
        } else {
            std::cout << "Should not end up here, EVER" << std::endl;
        }
    }
}

BOOST_AUTO_TEST_CASE(mutate_active_speed)
{
    dcgp::kernel_set<double> basic_set({"sum", "diff", "mul", "div"});

    std::cout << "Active Mutations: ONE BY ONE" << std::endl;
    perform_mutations_one_by_one(2, 4, 2, 3, 4, 2, 1000, basic_set());
    perform_mutations_one_by_one(2, 4, 10, 10, 1, 2, 1000, basic_set());
    perform_mutations_one_by_one(2, 4, 20, 20, 1, 2, 1000, basic_set());
    perform_mutations_one_by_one(1, 1, 100, 100, 1, 2, 100, basic_set());
    perform_mutations_one_by_one(5000, 1, 100, 100, 1, 2, 100, basic_set());
    std::cout << "\nActive Mutations: AT ONCE" << std::endl;
    perform_mutations_at_once(2, 4, 2, 3, 4, 2, 1000, basic_set());
    perform_mutations_at_once(2, 4, 10, 10, 1, 2, 1000, basic_set());
    perform_mutations_at_once(2, 4, 20, 20, 1, 2, 1000, basic_set());
    perform_mutations_at_once(1, 1, 100, 100, 1, 2, 100, basic_set());
    perform_mutations_at_once(5000, 1, 100, 100, 1, 2, 100, basic_set());
    std::cout << "------------------" << std::endl;
}

BOOST_AUTO_TEST_CASE(mutate_connections_speed)
{
    dcgp::kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::string kind("connection");

    std::cout << "Mutating active connections: ONE BY ONE" << std::endl;
    perform_mutations_one_by_one(2, 4, 2, 3, 4, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(2, 4, 10, 10, 1, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(2, 4, 20, 20, 1, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(1, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    perform_mutations_one_by_one(5000, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    std::cout << "\nMutating active connections: AT ONCE" << std::endl;
    perform_mutations_at_once(2, 4, 2, 3, 4, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(2, 4, 10, 10, 1, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(2, 4, 20, 20, 1, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(1, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    perform_mutations_at_once(5000, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    std::cout << "------------------" << std::endl;
}

BOOST_AUTO_TEST_CASE(mutate_function_speed)
{
    dcgp::kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::string kind("function");

    std::cout << "Mutating active function genes: ONE BY ONE" << std::endl;
    perform_mutations_one_by_one(2, 4, 2, 3, 4, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(2, 4, 10, 10, 1, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(2, 4, 20, 20, 1, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(1, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    perform_mutations_one_by_one(5000, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    std::cout << "\nMutating active function genes: AT ONCE" << std::endl;
    perform_mutations_at_once(2, 4, 2, 3, 4, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(2, 4, 10, 10, 1, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(2, 4, 20, 20, 1, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(1, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    perform_mutations_at_once(5000, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    std::cout << "------------------" << std::endl;
}


BOOST_AUTO_TEST_CASE(mutate_random_speed)
{
    dcgp::kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::string kind("random");

    std::cout << "Mutating random genes: ONE BY ONE" << std::endl;
    perform_mutations_one_by_one(2, 4, 2, 3, 4, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(2, 4, 10, 10, 1, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(2, 4, 20, 20, 1, 2, 1000, basic_set(), kind);
    perform_mutations_one_by_one(1, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    perform_mutations_one_by_one(5000, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    std::cout << "\nMutating random genes: AT ONCE" << std::endl;
    perform_mutations_at_once(2, 4, 2, 3, 4, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(2, 4, 10, 10, 1, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(2, 4, 20, 20, 1, 2, 1000, basic_set(), kind);
    perform_mutations_at_once(1, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    perform_mutations_at_once(5000, 1, 100, 100, 1, 2, 100, basic_set(), kind);
    std::cout << "------------------" << std::endl;
}