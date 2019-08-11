#ifndef GP_ARCHIVE_HPP
#define GP_ARCHIVE_HPP

#include <limbo/stat/stat_base.hpp>

namespace limbo {
    namespace stat {
        /// @ingroup stat
        /// filename: `gp_archive_<iteration>.dat`
        template <typename Params>
        struct GPArchive : public StatBase<Params> {
            template <typename BO, typename AggregatorFunction>
            void operator()(const BO& bo, const AggregatorFunction& afun)
            {
                if (!bo.stats_enabled())
                    return;
                std::string fname = bo.res_dir() + "/" + "gp_archive_" + std::to_string(bo.total_iterations()) + ".dat";
                std::ofstream ofs(fname.c_str());
                int gp_in = bo.model().dim_in();
                int gp_out = bo.model().dim_out();
                ofs << "#Point[" << gp_in << "d] mu[" << gp_out << "d] sigma[1d] acquisition[1d]" << std::endl;
                // for each point from iterator, sample query bo
                typedef typename Params::archiveparams::archive_t::const_iterator archive_it_t;

                for (archive_it_t it = Params::archiveparams::archive.begin(); it != Params::archiveparams::archive.end(); ++it) {
                    Eigen::VectorXd point(it->first.size());
                    for (size_t i = 0; i < it->first.size(); i++)
                        point[i] = it->first[i];

                    auto q = bo.model().query(point);

                    double acqui = opt::fun(typename BO::acquisition_function_t(bo.model(), bo.current_iteration())(point, afun, false));

                    ofs << point.transpose() << " "
                        << std::get<0>(q).transpose() << " "
                        << std::get<1>(q) << " "
                        << acqui << std::endl;
                }
            }
        };
    }
}

#endif // GP_ARCHIVE_HPP
