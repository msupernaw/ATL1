#ifndef ATL_FUNCTION_MINIMIZER
#define ATL_FUNCTION_MINIMIZER

#include "../AutoDiff/Variable.hpp"
#include <valarray>
#include <vector>
#include <cmath>
namespace atl {

    template<typename REAL_T>
    class FunctionMinimizer {
    public:

        enum Routine {
            LBFGS = 0,
            BFGS,
            NEWTON,
            NEWTON_CG
        };
    protected:
        std::vector<atl::Variable<REAL_T>* > parameters;
        std::vector<atl::Variable<REAL_T>* > active_parameters;
    private:
        std::valarray<REAL_T> x;
        std::valarray<REAL_T> best;
        std::valarray<REAL_T> gradient;
        std::valarray<std::valarray<REAL_T> > hessian;
        std::valarray<std::valarray<REAL_T> > identity;
        uint32_t iteration;
        uint32_t max_iterations;
        uint32_t max_line_searches;


        std::vector<unsigned int> phases;
        uint32_t max_phase;
        uint32_t current_phase;
        size_t sum_time_in_user_function;
        size_t average_time_in_user_function;
        size_t sum_time_in_grad_calc;
        size_t average_time_in_grad_calc;
        size_t function_calls;
        size_t gradient_calls;

        bool verbose;
        bool min_verbose;
        bool print_internal;
        uint32_t iprint;
        REAL_T maxgc;

        size_t max_history;
        uint32_t unrecorded_calls;


        REAL_T tolerance;
        std::vector<REAL_T> phase_tolerances;
        std::vector<uint32_t> phase_max_iterations;
        std::vector<uint32_t> phase_max_line_searches;
        std::vector<Routine> phase_routines;

        REAL_T function_value;
        REAL_T step;
        Routine routine;

    public:

        FunctionMinimizer() :
        routine(LBFGS),
        step(1.0),
        max_iterations(1000),
        max_line_searches(5000),
        max_phase(1),
        verbose(true),
        min_verbose(false),
        iprint(10),
        print_internal(false),
        maxgc(std::numeric_limits<REAL_T>::min()),
        max_history(1000),
        unrecorded_calls(0)/*,
        gradient_structure(&atl::Variable<REAL_T>::gradient_structure_g)*/ {

        }

        virtual void Initialize() {

        }

        virtual void ObjectiveFunction(atl::Variable<REAL_T>& f) {

        }

        virtual void Finalize() {

        }

        virtual void TransitionPhase() {

        }

        void SetRoutine(Routine r, uint32_t phase) {
            if (phase <= this->phase_routines.size) {
                this->phase_routines[(phase - 1)] = r;
            } else {
                this->phase_routines.resize(phase, LBFGS);
                this->phase_routines[(phase - 1)] = r;
            }
        }

        void Register(atl::Variable<REAL_T>& var, uint32_t phase = 1, std::string name = "") {
            this->parameters.push_back(&var);
            this->phases.push_back(phase);
            if (max_phase < phase) {
                max_phase = phase;
            }

            if (name != std::string("")) {
                var.SetName(name);
            }
        }

        void Register(atl::Matrix<atl::Variable<REAL_T> >&m, uint32_t phase = 1, std::string name = "") {
            if (max_phase < phase) {
                max_phase = phase;
            }
            if (name == std::string("")) {
                for (int i = 0; i < m.Size(0); i++) {
                    for (int j = 0; j < m.Size(1); j++) {
                        this->parameters.push_back(&m(i, j));
                        this->phases.push_back(phase);
                    }
                }
            } else {
                std::stringstream ss;
                for (int i = 0; i < m.Size(0); i++) {
                    for (int j = 0; j < m.Size(1); j++) {
                        ss << name << "(" << i << "," << j << ")";
                        m(i, j).SetName(ss.str());
                        this->parameters.push_back(&m(i, j));
                        this->phases.push_back(phase);
                        ss.str("");
                    }
                }
            }
        }

        void Register(atl::Vector<atl::Variable<REAL_T> >&m, uint32_t phase = 1, std::string name = "") {
            if (max_phase < phase) {
                max_phase = phase;
            }
            if (name == std::string("")) {
                for (int i = 0; i < m.Size(0); i++) {

                    this->parameters.push_back(&m(i));
                    this->phases.push_back(phase);

                }
            } else {
                std::stringstream ss;
                for (int i = 0; i < m.Size(0); i++) {
                    ss << name << "(" << i << ")";
                    m(i).SetName(ss.str());
                    this->parameters.push_back(&m(i));
                    this->phases.push_back(phase);
                    ss.str("");

                }
            }
        }

        void Register(atl::Array<atl::Variable<REAL_T> >&m, uint32_t phase = 1, std::string name = "") {
            if (max_phase < phase) {
                max_phase = phase;
            }
            size_t dims = m.Dimensions();

            if (name == std::string("")) {

                switch (dims) {

                    case 1:

                        break;

                    case 2:

                        break;

                    case 3:

                        break;

                    case 4:

                        break;

                    case 5:

                        break;

                    case 6:

                        break;

                    case 7:

                        break;

                }


            } else {
                std::stringstream ss;

            }
        }

        void Run() {
            for (int p = 1; p <= this->max_phase; p++) {

                Routine r = routine;
                set_defaults();
                this->active_parameters.resize(0);
                for (int i = 0; i < this->parameters.size(); i++) {
                    this->current_phase = p;
                    if (phases[i] <= p) {
                        this->active_parameters.push_back(parameters[i]);
                    }
                }

                if (this->phase_max_iterations.size() > 0 && this->phase_max_iterations.size() >= p) {
                    this->max_iterations = this->phase_max_iterations[(p - 1)];

                }

                if (this->phase_max_line_searches.size() > 0 && this->phase_max_line_searches.size() >= p) {
                    this->max_line_searches = this->phase_max_line_searches[(p - 1)];
                }

                if (this->phase_tolerances.size() > 0 && this->phase_tolerances.size() >= p) {
                    this->tolerance = this->phase_tolerances[(p - 1)];
                }


                if (this->phase_routines.size() > 0 && this->phase_routines.size() >= p) {
                    routine = this->phase_routines[(p - 1)];
                } else {
                    routine = LBFGS;
                }

                switch (routine) {
                    case LBFGS:
                        lbfgs(this->tolerance);
                        break;
                    case NEWTON:

                        newton(this->tolerance);
                        break;
                    case NEWTON_CG:
                        newton_cg(this->tolerance);
                        break;
                    default:
                        lbfgs(this->tolerance);
                        break;

                }
                this->TransitionPhase();
                routine = r;
            }
        }

        std::vector<uint32_t>& GetPhaseMaxIterations() {
            return phase_max_iterations;
        }

        void SetPhaseMaxIterations(std::vector<uint32_t> phase_max_iterations) {
            this->phase_max_iterations = phase_max_iterations;
        }

        std::vector<uint32_t>& GetPhaseMaxLineSearches() {
            return phase_max_line_searches;
        }

        void SetPhaseMaxLineSearches(std::vector<uint32_t> phase_max_line_searches) {
            this->phase_max_line_searches = phase_max_line_searches;
        }

        std::vector<REAL_T>& GetPhaseTolerances() {
            return phase_tolerances;
        }

        void SetPhaseTolerances(std::vector<REAL_T> phase_tolerances) {
            this->phase_tolerances = phase_tolerances;
        }

        std::vector<Routine>& GetPhaseRoutines() {
            return phase_routines;
        }

        void SetMaxPhase(uint32_t max) {
            this->max_phase = max;
        }

        void SetPhaseRoutines(std::vector<Routine> phase_routines) {
            this->phase_routines = phase_routines;
        }
    private:

        void set_defaults() {
            max_iterations = (1000);
            max_line_searches = (1000);
            this->tolerance = 1e-4;


        }

        bool lbfgs(REAL_T tolerance = 1e-4) {
            atl::Variable<REAL_T>::SetRecording(true);
            atl::Variable<REAL_T>::gradient_structure_g.Reset();
            int nops = this->active_parameters.size();
            this->x.resize(nops);
            this->best.resize(nops);
            this->gradient.resize(nops);
            this->hessian.resize(nops, std::valarray<REAL_T>(nops));

            for (int i = 0; i < nops; i++) {
                if (this->active_parameters[i]->IsBounded()) {
                    this->x[i] = this->active_parameters[i]->GetInternalValue();
                } else {
                    this->x[i] = this->active_parameters[i]->GetValue();
                }
                this->gradient[i] = 0;
                for (int j = 0; j < nops; j++) {
                    this->hessian[i][j] = 0;
                }
            }
            //
            //
            std::valarray<REAL_T> wg(nops);
            std::valarray<REAL_T> nwg(nops);
            std::valarray<REAL_T> ng(nops);

            //initial evaluation
            atl::Variable<REAL_T> fx(0.0);
            this->call_objective_function(fx);
            this->function_value = fx.GetValue();
            //
            //Historical evaluations
            std::valarray<REAL_T> px(nops);
            std::valarray<REAL_T> pg(nops);
            std::valarray<std::valarray<REAL_T> > dxs(std::valarray<REAL_T > (max_history), nops);
            std::valarray<std::valarray<REAL_T> > dgs(std::valarray<REAL_T > (max_history), nops);
            //search direction
            std::valarray<REAL_T> z(nops);

            this->get_gradient();

            step = 1.0; //0.0001;
            REAL_T relative_tolerance;
            REAL_T norm_g;
            std::valarray<REAL_T> p(max_history);
            std::valarray<REAL_T>a(max_history);

            int i;
            for (this->iteration = 0; this->iteration < this->max_iterations; this->iteration++) {
                i = this->iteration;
                for (int j = 0; j < nops; j++) {
                    wg[j] = active_parameters[j]->GetScaledGradient(active_parameters[j]->GetValue()) * gradient[j];
                }

                //
                if (this->verbose) {
                    if ((i % this->iprint) == 0) {
                        this->Print2(fx, this->gradient, this->active_parameters);
                    }
                }

                if (this->maxgc < tolerance) {
                    //                    if (this->verbose) {
                    atl::Variable<REAL_T> fx(0.0);
                    this->call_objective_function(fx);
                    this->function_value = fx.GetValue();
                    this->get_gradient();
                    this->Print2(fx, this->gradient, this->active_parameters);
                    //                    }
                    return true;
                }

                z = wg;

                if (i > 0 && max_history > 0) {

                    size_t h = std::min<size_t > (i, max_history);
                    size_t end = (i - 1) % h;

                    //update histories
                    for (size_t r = 0; r < nops; r++) {
                        // std::cout<<r<<" "<<g[r] <<" - "<< pg[r]<<"\n";
                        dxs[r][end] = active_parameters[r]->GetInternalValue() - px[r];
                        dgs[r][end] = wg[r] - pg[r];
                    }



                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end - j + h) % h;
                        p[k] = 1.0 / Dot(Column(dxs, k), Column(dgs, k));

                        a[k] = p[k] * Dot(Column(dxs, k), z);
                        z -= a[k] * this->Column(dgs, k);
                    }
                    // Scaling of initial Hessian (identity matrix)
                    z *= Dot(this->Column(dxs, end), Column(dgs, end)) / Dot(Column(dgs, end), Column(dgs, end));

                    for (size_t j = 0; j < h; ++j) {
                        const size_t k = (end + j + 1) % h;
                        const REAL_T b = p[k] * Dot(Column(dgs, k), z);
                        z += this->Column(dxs, k) * (a[k] - b);
                    }

                }//end if(i>0)

                for (size_t j = 0; j < nops; j++) {

                    px[j] = active_parameters[j]->GetInternalValue();
                    x[j] = px[j];
                    pg[j] = wg[j];

                }//end for

                if (!this->line_search(z, wg, i)) {
                    std::cout << "Max line serches.";
                    return false;
                }

            }
            return true;
        }

        bool bfgs(REAL_T tolerance = 1e-4) {

            atl::Variable<REAL_T>::SetRecording(true);
            atl::Variable<REAL_T>::gradient_structure_g.Reset();
            atl::Variable<REAL_T>::gradient_structure_g.derivative_trace_level = atl::GRADIENT_AND_HESSIAN;


            int nops = this->parameters.size();

            std::valarray<std::valarray<REAL_T> > system(nops);
            //            std::valarray<REAL_T> update;
            std::valarray<REAL_T> wg(nops);
            std::valarray<REAL_T> neg_g(nops);
            REAL_T old_fval;
            REAL_T old_old_fval;
            REAL_T eta;
            std::valarray<std::valarray<REAL_T> > B;
            this->x.resize(nops);
            this->best.resize(nops);
            this->gradient.resize(nops);
            this->hessian.resize(nops, std::valarray<REAL_T>(nops));
            B.resize(nops, std::valarray<REAL_T>(nops));


            for (int i = 0; i < nops; i++) {
                system[i].resize(nops);
                if (this->active_parameters[i]->IsBounded()) {
                    this->x[i] = this->active_parameters[i]->GetInternalValue();
                } else {
                    this->x[i] = this->active_parameters[i]->GetValue();
                }
                this->gradient[i] = 0;
                for (int j = 0; j < nops; j++) {
                    this->hessian[i][j] = 0;
                }
            }
            //initial evaluation
            atl::Variable<REAL_T> fx(0.0);
            this->call_objective_function(fx);
            this->function_value = fx.GetValue();
            this->get_gradient_and_hessian();

            for (int i = 0; i < nops; i++) {
                for (int j = 0; j < nops; j++) {
                    B[i][j] = this->hessian[i][j];
                }
            }

            int i;
            for (this->iteration = 0; this->iteration < this->max_iterations; this->iteration++) {
                i = this->iteration;


                //
                if (this->verbose) {
                    if ((i % this->iprint) == 0) {
                        this->Print2(fx, this->gradient, this->active_parameters);
                    }
                }


                for (int j = 0; j < nops; j++) {
                    wg[j] = active_parameters[j]->GetScaledGradient(active_parameters[j]->GetValue()) * gradient[j];
                    neg_g[j] = -wg[j];
                }



            }



        }

        bool newton(REAL_T tolerance = 1e-4) {
            atl::Variable<REAL_T>::SetRecording(true);
            atl::Variable<REAL_T>::gradient_structure_g.Reset();
            atl::Variable<REAL_T>::gradient_structure_g.derivative_trace_level = atl::GRADIENT_AND_HESSIAN;


            int nops = this->active_parameters.size();

            std::valarray<std::valarray<REAL_T> > system(nops);
            std::valarray<REAL_T> wg(nops);

            REAL_T old_fval;
            REAL_T old_old_fval;
            REAL_T eta;
            this->x.resize(nops);
            this->best.resize(nops);
            this->gradient.resize(nops);
            this->hessian.resize(nops, std::valarray<REAL_T>(nops));
            this->identity.resize(nops, std::valarray<REAL_T>(nops));
            for (int i = 0; i < nops; i++) {
                if (this->active_parameters[i]->IsBounded()) {
                    this->x[i] = this->active_parameters[i]->GetInternalValue();
                } else {
                    this->x[i] = this->active_parameters[i]->GetValue();
                }
                this->gradient[i] = 0;
                for (int j = 0; j < nops; j++) {
                    this->hessian[i][j] = 0;
                }
            }



            std::valarray<REAL_T> b(nops);
            std::valarray<REAL_T> pk(nops);
            std::valarray<REAL_T> update(nops);
            std::valarray<REAL_T> gfk(nops);
            std::valarray<REAL_T> ri(nops);
            std::valarray<REAL_T> psupi(nops);
            std::valarray<REAL_T> xsupi(nops);
            std::valarray<REAL_T> Ap(nops);

            atl::Variable<REAL_T> fx(0.0);
            this->call_objective_function(fx);
            this->get_gradient_and_hessian();
            this->function_value = fx.GetValue();



            //Historical evaluations
            std::valarray<REAL_T> px(nops);
            std::valarray<REAL_T> pg(nops);
            std::valarray<std::valarray<REAL_T> > dxs(std::valarray<REAL_T > (max_history), nops);
            std::valarray<std::valarray<REAL_T> > dgs(std::valarray<REAL_T > (max_history), nops);
            //search direction
            std::valarray<REAL_T> z(nops);

            step = 1.0; //0.0001;

            std::valarray<REAL_T> p(max_history);
            std::valarray<REAL_T>a(max_history);


            for (int j = 0; j < nops; j++) {
                system[j] = std::valarray<REAL_T>(nops + 1);
                for (int k = 0; k < nops; k++) {
                    if (j == k) {
                        identity[j][j] = 1.0;
                    } else {
                        identity[j][j] = 0;
                    }
                }
            }

            int i;
            for (this->iteration = 0; this->iteration < this->max_iterations; this->iteration++) {
                i = this->iteration;


                //
                if (this->verbose) {
                    if ((i % this->iprint) == 0) {
                        //                    this->print();
                        this->Print2(fx, this->gradient, this->active_parameters);
                    }
                }


                for (int j = 0; j < nops; j++) {
                    wg[j] = active_parameters[j]->GetScaledGradient(active_parameters[j]->GetValue()) * gradient[j];
                }

                if (i > 0 && max_history > 0) {

                    size_t h = std::min<size_t > (i, max_history);
                    size_t end = (i - 1) % h;

                    //update histories
                    for (size_t r = 0; r < nops; r++) {
                        // std::cout<<r<<" "<<g[r] <<" - "<< pg[r]<<"\n";
                        dxs[r][end] = active_parameters[r]->GetInternalValue() - px[r];
                        dgs[r][end] = wg[r] - pg[r];
                    }
                }



                //try standard newton update with hessian, else use estimated h.
                if (!this->cholesky_solve(hessian, wg, update)) {

                    z = wg;
                    if (i > 0 && max_history > 0) {
                        size_t h = std::min<size_t > (i, max_history);
                        size_t end = (i - 1) % h;
                        for (size_t j = 0; j < h; ++j) {
                            const size_t k = (end - j + h) % h;
                            p[k] = 1.0 / Dot(Column(dxs, k), Column(dgs, k));

                            a[k] = p[k] * Dot(Column(dxs, k), z);
                            z -= a[k] * Column(dgs, k);
                        }
                        // Scaling of initial Hessian (identity matrix)
                        z *= Dot(Column(dxs, end), Column(dgs, end)) / Dot(Column(dgs, end), Column(dgs, end));

                        for (size_t j = 0; j < h; ++j) {
                            const size_t k = (end + j + 1) % h;
                            const REAL_T b = p[k] * Dot(Column(dgs, k), z);
                            z += Column(dxs, k) * (a[k] - b);
                        }
                        update = z;
                    } else {
                        //                    //                    std::cout<<"not pd!!!\n";
                        for (int j = 0; j < nops; j++) {
                            //                        update[j] = 0;
                            //                        update[j]+=wg[j]*.001*this->step;
                            for (int k = 0; k < nops; k++) {
                                if (j == k) {
                                    identity[j][j] = 1;
                                } else {
                                    identity[j][j] = 0;
                                }
                                system[j][k] = identity[j][k];
                            }
                            system[j][nops] = wg[j];
                        }
                        this->gauss(system, update);
                    }
                } else {
                    std::cout << "was pd!!!\n";
                }

                for (size_t j = 0; j < nops; j++) {

                    px[j] = active_parameters[j]->GetInternalValue();
                    x[j] = px[j];
                    pg[j] = wg[j];

                }

                if (!this->line_search(update, wg, i)) {
                    std::cout << "max line searches...." << std::endl;
                    return false;
                }


                if (maxgc < tolerance) {
                    std::cout << "Successful Convergence.\n";
                    this->Print2(fx, this->gradient, this->active_parameters);
                    std::cout << "Successful Convergence." << std::endl;
                    return true;
                }


            }
            std::cout << "Max iterations.\n";
            return false;
        }

        bool newton_cg(REAL_T tolerance = 1e-4) {
            atl::Variable<REAL_T>::SetRecording(true);
            atl::Variable<REAL_T>::gradient_structure_g.Reset();
            atl::Variable<REAL_T>::gradient_structure_g.derivative_trace_level = atl::GRADIENT_AND_HESSIAN;


            int nops = this->active_parameters.size();

            std::valarray<std::valarray<REAL_T> > system(nops);
            //            std::valarray<REAL_T> update;
            std::valarray<REAL_T> wg(nops);

            REAL_T old_fval;
            REAL_T old_old_fval;
            REAL_T eta;
            this->x.resize(nops);
            this->best.resize(nops);
            this->gradient.resize(nops);
            this->hessian.resize(nops, std::valarray<REAL_T>(nops));

            for (int i = 0; i < nops; i++) {
                system[i].resize(nops);
                if (this->active_parameters[i]->IsBounded()) {
                    this->x[i] = this->active_parameters[i]->GetInternalValue();
                } else {
                    this->x[i] = this->active_parameters[i]->GetValue();
                }
                this->gradient[i] = 0;
                for (int j = 0; j < nops; j++) {
                    this->hessian[i][j] = 0;
                }
            }



            std::valarray<REAL_T> b(nops);
            std::valarray<REAL_T> pk(nops);
            std::valarray<REAL_T> update(nops);
            std::valarray<REAL_T> gfk(nops);
            std::valarray<REAL_T> ri(nops);
            std::valarray<REAL_T> psupi(nops);
            std::valarray<REAL_T> xsupi(nops);
            std::valarray<REAL_T> Ap(nops);

            atl::Variable<REAL_T> fx(0.0);
            this->call_objective_function(fx);
            this->get_gradient_and_hessian();
            this->function_value = fx.GetValue();
            REAL_T maggrad;

            int i;
            for (this->iteration = 0; this->iteration < this->max_iterations; this->iteration++) {
                i = this->iteration;


                //
                if (this->verbose) {
                    if ((i % this->iprint) == 0) {
                        //                    this->print();
                        this->Print2(fx, this->gradient, this->active_parameters);
                    }
                }


                for (int j = 0; j < nops; j++) {
                    wg[j] = active_parameters[j]->GetScaledGradient(active_parameters[j]->GetValue()) * gradient[j];
                }

                maggrad = norm(wg);
                eta = std::min(static_cast<REAL_T> (.5), std::sqrt(maggrad));
                b = REAL_T(-1) * wg;
                REAL_T termcond = eta * maggrad;
                ri = REAL_T(-1) * b;
                xsupi = std::valarray<REAL_T>(nops);
                psupi = REAL_T(-1.0) * ri;
                int i = 0;
                REAL_T dri0 = 0;
                for (int ii = 0; ii < nops; ii++) {
                    dri0 += ri[ii] * ri[ii];
                }


                while (abs_sum(ri) > termcond) {
                    REAL_T alphai;
                    REAL_T betai;
                    REAL_T curv = 0.0;
                    for (int ii = 0; ii < nops; ii++) {
                        Ap[ii] = 0;
                        for (int j = 0; j < nops; j++) {
                            Ap[ii] += psupi[j] * hessian[ii][j];

                        }
                        curv += Ap[ii] * psupi[ii];
                    }

                    if (0 <= curv && curv <= 3.0 * std::numeric_limits<REAL_T>::epsilon()) {
                        break;
                    } else if (curv < 0) {
                        if (i > 0) {
                            break;
                        } else {
                            xsupi = dri0 / (-curv) * b;
                            break;
                        }
                    }

                    alphai = dri0 / curv;
                    xsupi = xsupi + alphai * psupi;
                    ri = ri + alphai * Ap;

                    REAL_T dri1 = 0;
                    for (int i = 0; i < nops; i++) {
                        dri1 += ri[i] * ri[i];
                    }
                    betai = dri1 / dri0;
                    psupi = -ri + betai * psupi;
                    i = i + 1;
                    dri0 = dri1;
                }


                pk = xsupi; //search direction is solution to system.
                //                                gradient = -b; //gradient at xk




                if (!this->line_search(pk, wg, i)) {
                    std::cout << "max line searches...." << std::endl;
                    return false;
                }

                if (maxgc < tolerance) {
                    std::cout << "Successful Convergence.\n";
                    this->Print2(fx, this->gradient, this->active_parameters);
                    std::cout << "Successful Convergence." << std::endl;
                    return true;
                }


            }
            std::cout << "Max iterations.\n";
            return false;
        }

        bool line_search(std::valarray<REAL_T>& z, std::valarray<REAL_T>& wg, int& i) {
            REAL_T descent = 0;

            int nops = this->active_parameters.size();
            std::valarray<REAL_T> nwg(nops);
            std::valarray<REAL_T> ng(nops);

            for (size_t j = 0; j < nops; j++) {
                descent += z[j] * wg[j];
            }//end for

            REAL_T norm_g = this->norm(gradient);
            REAL_T relative_tolerance = this->tolerance * std::max<REAL_T > (REAL_T(1.0), norm_g);

            descent *= REAL_T(-1.0); // * Dot(z, g);
            if (descent > REAL_T(-0.0000000001) * relative_tolerance /* tolerance relative_tolerance*/) {
                z = wg;
                //                this->max_iterations -= i;
                i = 0;
                descent = -1.0 * Dot(z, wg);
            }//end if

            step = 1.0;


            bool down = false;

            int ls;


            //line_search:
            atl::Variable<REAL_T>::SetRecording(false);
            for (int j = 0; j < active_parameters.size(); j++) {
                best[j] = active_parameters[j]->GetValue();
            }

            atl::Variable<REAL_T> fx;
            for (ls = 0; ls < this->max_line_searches; ++ls) {



                // Tentative solution, gradient and loss
                std::valarray<REAL_T> nx = x - step * z;

                for (size_t j = 0; j < nops; j++) {
                    active_parameters[j]->UpdateValue(nx[j]);
                }

                this->call_objective_function(fx);

                if (fx.GetValue() <= this->function_value + tolerance * REAL_T(0.000001) * (1.0 / norm_g) * step * descent) { // First Wolfe condition

                    atl::Variable<REAL_T>::SetRecording(true);
                    this->call_objective_function(fx);
                    this->get_gradient(ng);

                    if (down || (-1.0 * Dot(z, nwg) >= 0.9 * descent)) { // Second Wolfe condition
                        x = nx;
                        gradient = ng;
                        this->function_value = fx.GetValue();
                        return true;
                    } else {
                        atl::Variable<REAL_T>::SetRecording(false);
                        step *= 10.0;
                    }
                } else {
                    step /= 10.0;
                    down = true;
                }
            }
            for (size_t j = 0; j < nops; j++) {
                active_parameters[j]->SetValue(best[j]);
            }

            return false;

        }

        void print22() {
            const char separator = ' ';
            const int idWidth = 8;
            const int nameWidth = 10;
            const int valWidth = 15;

            std::cout << "Phase: " << this->current_phase << " of " << this->max_phase << std::endl;
            std::cout << "Iteration: " << this->iteration << std::endl;
            std::cout << "Convergence Criteria: " << this->tolerance << std::endl;
            std::cout << "Routine: ";
            switch (this->routine) {
                case LBFGS:
                    std::cout << "L-BFGS" << std::endl;
                    break;
                case NEWTON_CG:
                    std::cout << "Newton CG" << std::endl;
                    break;
                case NEWTON:
                    std::cout << "Newton Hybrid" << std::endl;
                    break;
                default:
                    std::cout << "L-BFGS(D)" << std::endl;
                    break;

            }
            std::cout << "Function Value: " << function_value << std::endl;
            std::cout << "Max Gradient:" << maxgc << std::endl;
            if (!min_verbose) {
                std::cout << std::left << std::setw(idWidth) << std::setfill(separator) << "Id";
                std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Name";
                std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Value";
                std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Gradient";
                std::cout << "|";
                std::cout << std::left << std::setw(idWidth) << std::setfill(separator) << "Id";
                std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Name";
                std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Value";
                std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Gradient";
                std::cout << std::endl;


                int stop = this->active_parameters.size() - 2;
                int i;

                for (i = 0; i < stop; i += 2) {

                    if (i > 0 && (i % 20) == 0) {
                        std::cout << "Iteration: " << this->iteration << std::endl;
                        std::cout << "Phase: " << this->current_phase << " of " << this->max_phase << std::endl;
                        std::cout << "Convergence Criteria: " << this->tolerance << std::endl;
                        std::cout << "Routine: ";
                        switch (this->routine) {
                            case LBFGS:
                                std::cout << "L-BFGS" << std::endl;
                                break;
                            case NEWTON_CG:
                                std::cout << "Newton CG" << std::endl;
                                break;
                            default:
                                std::cout << "L-BFGS" << std::endl;
                                break;

                        }
                        std::cout << "Function Value: " << function_value << std::endl;
                        std::cout << "Max Gradient:" << maxgc << std::endl;

                        std::cout << std::left << std::setw(idWidth) << std::setfill(separator) << "Id";
                        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Name";
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Value";
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Gradient";
                        std::cout << "|";
                        std::cout << std::left << std::setw(idWidth) << std::setfill(separator) << "Id";
                        std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << "Name";
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Value";
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << "Gradient";
                        std::cout << std::endl;

                    }

                    std::cout << std::left << std::setw(idWidth) << std::setfill(separator) << this->active_parameters[i]->info->id;
                    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << this->active_parameters[i]->GetName();
                    if (this->print_internal) {
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->active_parameters[i]->GetInternalValue();

                    } else {
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->active_parameters[i]->GetValue();

                    }
                    std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->gradient[i];
                    std::cout << "|";
                    std::cout << std::left << std::setw(idWidth) << std::setfill(separator) << this->active_parameters[i + 1]->info->id;
                    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << this->active_parameters[i + 1]->GetName();
                    if (this->print_internal) {
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->active_parameters[i + 1]->GetInternalValue();

                    } else {
                        std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->active_parameters[i + 1]->GetValue();

                    }
                    std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->gradient[i + 1];
                    std::cout << std::endl;
                }

                for (; i < this->active_parameters.size(); i++) {
                    std::cout << std::left << std::setw(idWidth) << std::setfill(separator) << this->active_parameters[i]->info->id;
                    std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) << this->active_parameters[i]->GetName();
                    std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->active_parameters[i]->GetValue();
                    std::cout << std::left << std::setw(valWidth) << std::setfill(separator) << this->gradient[i];
                    //                if (i < (this->parameters.size() - 1)) {
                    //                    
                    //                                    } 
                }
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }





    private:

        void call_objective_function(atl::Variable<REAL_T>& f) {
            atl::Variable<REAL_T>::gradient_structure_g.Reset();
            this->ObjectiveFunction(f);
        }

        void get_gradient() {
            this->maxgc = 100000;

            atl::Variable<REAL_T>::gradient_structure_g.Accumulate();
            for (int i = 0; i < this->active_parameters.size(); i++) {

                this->gradient[i] = this->active_parameters[i]->info->dvalue;
                if (i == 0) {
                    maxgc = std::fabs(this->gradient[i]);
                } else if (std::fabs(this->gradient[i]) > maxgc) {
                    maxgc = std::fabs(this->gradient[i]);
                }
            }
        }

        void get_gradient(std::valarray<REAL_T>& g) {
            this->maxgc = 100000;
            atl::Variable<REAL_T>::gradient_structure_g.Accumulate();
            for (int i = 0; i < this->active_parameters.size(); i++) {
                g[i] = this->active_parameters[i]->info->dvalue;
                if (i == 0) {
                    maxgc = std::fabs(g[i]);
                } else if (std::fabs(g[i]) > maxgc) {
                    maxgc = std::fabs(g[i]);
                }
            }
        }

        void get_gradient_and_hessian() {

            this->maxgc = 100000;
            atl::Variable<REAL_T>::gradient_structure_g.HessianAndGradientAccumulate();
            for (int i = 0; i < this->active_parameters.size(); i++) {
                gradient[i] = this->active_parameters[i]->info->dvalue;
                if (i == 0) {
                    maxgc = std::fabs(gradient[i]);
                } else if (std::fabs(gradient[i]) > maxgc) {
                    maxgc = std::fabs(gradient[i]);
                }
                for (int j = 0; j < this->active_parameters.size(); j++) {
                    hessian[i][j] = this->active_parameters[i]->info->hessian_row[this->active_parameters[j]->info];
                }
            }

        }

        bool cholesky_decomp(std::valarray< std::valarray<REAL_T> >& el) {
            int n = el.size();

            int i, j, k;

            REAL_T sum;

            for (i = 0; i < n; i++) {
                for (j = i; j < n; j++) {
                    for (sum = el[i][j], k = i - 1; k >= 0; k--)
                        sum -= el[i][k] * el[j][k];
                    if (i == j) {
                        if (sum <= 0.0) { // A, with rounding errors, is not
                            // positive-definite.
                            //                            std::cout << "Cholesky failed: A not positive definite\n";
                            return false;
                        }
                        el[i][i] = ::sqrt(sum);
                    } else
                        el[j][i] = sum / el[i][i];
                }
            }
            for (i = 0; i < n; i++)
                for (j = 0; j < i; j++)
                    el[j][i] = 0.;

            return true;
        }

        /**
         *
         * 
         * 
         * 
         * @param A
         * @param b
         * @param x
         * @return 
         */
        bool cholesky_solve(std::valarray< std::valarray<REAL_T> >& A, std::valarray<REAL_T>& b, std::valarray<REAL_T>& x) {

            int i, j, k;
            int n = A.size();
            REAL_T sum = 0;

            if (!cholesky_decomp(A)) {
                return false;
            }

            for (i = 0; i < n; i++) { // Solve L  y D b, storing y in x.
                for (sum = b[i], k = i - 1; k >= 0; k--)
                    sum -= A[i][k] * x[k];
                x[i] = sum / A[i][i];
            }
            for (i = n - 1; i >= 0; i--) { // Solve LT  x D y.
                for (sum = x[i], k = i + 1; k < n; k++)
                    sum -= A[k][i] * x[k];
                x[i] = sum / A[i][i];
            }
            return true;
        }

        void gauss(std::valarray< std::valarray<REAL_T> >& A, std::valarray<REAL_T>& x) {
            int n = A.size();
            x.resize(n);
            for (int i = 0; i < n; i++) {
                // Search for maximum in this column
                REAL_T maxEl = std::fabs(A[i][i]);
                int maxRow = i;
                for (int k = i + 1; k < n; k++) {
                    if (std::fabs(A[k][i]) > maxEl) {
                        maxEl = std::fabs(A[k][i]);
                        maxRow = k;
                    }
                }

                // Swap maximum row with current row (column by column)
                for (int k = i; k < n + 1; k++) {
                    REAL_T tmp = A[maxRow][k];
                    A[maxRow][k] = A[i][k];
                    A[i][k] = tmp;
                }

                // Make all rows below this one 0 in current column
                for (int k = i + 1; k < n; k++) {
                    REAL_T c = -A[k][i] / A[i][i];
                    for (int j = i; j < n + 1; j++) {
                        if (i == j) {
                            A[k][j] = 0;
                        } else {
                            A[k][j] += c * A[i][j];
                        }
                    }
                }
            }

            // Solve equation Ax=b for an upper triangular matrix A
            for (int i = n - 1; i >= 0; i--) {
                x[i] = A[i][n] / A[i][i];
                for (int k = i - 1; k >= 0; k--) {
                    A[k][n] -= A[k][i] * x[i];
                }
            }
        }

        /**
         * Compute the Norm of the vector v.
         *  
         * @param v
         * @return 
         */
        const REAL_T norm(std::valarray<REAL_T> &v) {

            REAL_T ret = (REAL_T) 0.0;
            unsigned int i;
            for (i = 0; i < v.size(); i++) {
                ret += v[i] * v[i];

            }
            return std::sqrt(ret);
        }

        /**
         * Compute the dot product of two vectors.
         * @param a
         * @param b
         * @return 
         */
        const REAL_T Dot(const std::valarray<REAL_T> &a, const std::valarray<REAL_T> &b) {
            REAL_T ret = 0;
            for (size_t i = 0; i < a.size(); i++) {
                ret += a[i] * b[i];
            }
            return ret;
        }

        /**
         * returns the a column of a matrix as a std::valarray.
         * @param matrix
         * @param column
         * @return 
         */
        const std::valarray<REAL_T> Column(std::valarray<std::valarray<REAL_T> > &matrix, size_t column) {

            std::valarray<REAL_T> ret(this->active_parameters.size());

            for (int i = 0; i < ret.size(); i++) {
                ret[i] = matrix[i][column];
            }
            return ret;
        }

        REAL_T sum(const std::valarray<REAL_T>&array) {
            REAL_T sum = 0.0;
            for (int i = 0; i < array.size(); i++) {
                sum += array[i];
            }
            return sum;
        }

        REAL_T abs_sum(const std::valarray<REAL_T>&array) {
            REAL_T sum = 0.0;
            for (int i = 0; i < array.size(); i++) {
                sum += std::fabs(array[i]);
            }
            return sum;
        }

        void make_identitiy(std::valarray<std::valarray<REAL_T> >& A) {
            for (int i = 0; i < A.size(); i++) {
                for (int j = 0; j < A[0].size(); j++) {
                    if (i == j) {
                        A[i][j] = 1.0;
                    } else {
                        A[i][j] = 0.0;
                    }
                }
            }
        }

        /* Convert double to string with specified number of places after the decimal. */
        const std::string prd(const double &x, const int decDigits) {
            std::stringstream ss;
            ss << std::fixed;
            ss.precision(decDigits); // set # places after decimal
            ss << x;
            return ss.str();
        }

        /* Convert double to string with specified number of places after the decimal
           and left padding. */
        const std::string prd(const double &x, const int decDigits, const int width) {
            std::stringstream ss;
            ss << std::scientific << std::right;
            ss.fill(' '); // fill space around displayed #
            ss.width(width); // set  width around displayed #
            ss.precision(decDigits); // set # places after decimal
            ss << x;
            return ss.str();
        }

        /*! Center-aligns string within a field of width w. Pads with blank spaces
            to enforce alignment. */
        const std::string center(const std::string &s, const int w) {
            std::stringstream ss, spaces;
            int padding = w - s.size(); // count excess room to pad
            for (int i = 0; i < padding / 2; ++i)
                spaces << " ";
            ss << spaces.str() << s << spaces.str(); // format with padding
            if (padding > 0 && padding % 2 != 0) // if odd #, add 1 space
                ss << " ";
            return ss.str();
        }

        /*! Center-aligns string within a field of width w. Pads with blank spaces
           to enforce alignment. */
        const std::string left(const std::string &s, const int w) {
            std::stringstream ss;
            std::string ret;
            if (s.size() >= w) {
                ret.insert(ret.begin(), s.begin(), s.begin()+(w - 4));
                ret += "...";
                ss << ret;
            } else {
                int dif = w - s.size();
                ss << s;
                for (int i = 0; i < dif; i++) {
                    ss << " ";
                }
            }

            return ss.str();
        }

        /**
         * Print current minimizer state to stdout.
         * 
         * @param ret
         * @param gradient
         * @param parameters
         * @param message
         */
        void Print2(const atl::Variable<REAL_T> &ret, const std::valarray<REAL_T> &gradient, const std::vector<atl::Variable<REAL_T>* > &parameters, std::string message = "") {


            std::cout << io::BOLD << message << io::DEFAULT << std::endl;
            std::cout << "Phase: " << io::BOLD << this->current_phase
                    << io::DEFAULT << " of " << io::BOLD << this->max_phase
                    << io::DEFAULT << "\n";
            std::cout << "Iteration: " << io::BOLD
                    << this->iteration << io::DEFAULT << std::endl;
            //            std::cout << "Function Calls: "
            //                    << io::BOLD << this->function_calls_m << io::DEFAULT << " (" << this->unrecorded_calls_m << " unrecorded line searches)" << std::endl;
            //            std::cout << "Average Time in Objective Function: "
            //                    << io::BOLD
            //                    << static_cast<long> (this->average_time_in_user_function_m)
            //                    << " ms\n" << io::DEFAULT;
            //            std::cout << "Average Time Calculating Gradients: " << io::BOLD
            //                    << static_cast<long> (this->average_time_in_grad_calc_m)
            //                    << " ms\n" << io::DEFAULT;
            int prec = std::cout.precision();
            std::cout.precision(50);
            std::cout << "Function Value = " << io::BOLD << ret.GetValue()
                    << std::endl << io::DEFAULT;
            std::cout.precision(prec);
            std::cout << "Active Parameters: " << io::BOLD << parameters.size()
                    << io::DEFAULT << std::endl;
            //            std::cout << std::scientific;
            //            std::cout << "Tolerance: " << io::BOLD << this->GetTolerance() << io::DEFAULT
            //                    << std::endl;
            //            std::cout << "Expression Size: " << BOLD << ret.ExpressionSize() << DEFAULT
            //                    << std::endl;
            std::cout << "Maximum Gradient Component Magnitude: " << io::BOLD
                    << std::scientific << maxgc << io::DEFAULT << std::endl;



            if (!this->min_verbose) {

                double number_good = 0.0;
                std::cout << std::fixed;
                using namespace util;



                for (size_t i = 0; i < this->active_parameters.size(); i += 2) {
                    //
                    //
                    if ((i & 30) == 0) {
                        std::cout << io::BOLD << "------------------------------------------------------------------------------------------------------------------\n|";

                        std::cout << io::BOLD << center("Id", 10) << "| " << center("Parameter", 16) << " | "
                                << center("Value", 10) << " | "
                                << center("Gradient", 10) << " | ";
                        //                    std::cout << center("Variable", 20) << " "
                        //                            << center("Value", 10) << " "
                        //                            << center("Gradient", 10) << " | ";
                        std::cout << io::BOLD << center("Id", 10) << "|  " << center("Parameter", 16) << " | "
                                << center("Value", 10) << " | "
                                << center("Gradient", 10) << " |\n" << io::DEFAULT;
                        std::cout << io::BOLD << "------------------------------------------------------------------------------------------------------------------\n" << io::DEFAULT;
                    }
                    for (int j = i; j < (i + 2); j++) {
                        if (j< this->active_parameters.size()) {


                            std::stringstream sp;
                            sp << this->active_parameters[j]->info->id;
                            std::cout << left(sp.str(), 10) << " | ";

                            if (this->active_parameters[j]->bounded_m) {
                                std::cout << io::BLUE << left(this->active_parameters[j]->GetName(), 16) << io::DEFAULT << " | ";
                            } else {
                                std::cout << left(this->active_parameters[j]->GetName(), 16) << " | ";
                            }
                            if (this->print_internal) {
                                std::cout << prd(this->active_parameters[j]->GetInternalValue(), 3, 10) << " | ";
                            } else {
                                std::cout << prd(this->active_parameters[j]->GetValue(), 3, 10) << " | ";
                            }
                            if (std::fabs(gradient[j])<this->tolerance) {
                                number_good++;
                                std::cout << io::GREEN << prd(gradient[j], 3, 10) << io::DEFAULT << io::BOLD << " | " << io::DEFAULT;
                            } else {
                                std::cout << io::RED << prd(gradient[j], 3, 10) << io::DEFAULT << io::BOLD << " | " << io::DEFAULT;
                            }

                        }
                    }
                    std::cout << std::endl;
                }



                //                        for (size_t i = 0; i < this->active_parameters_m.size(); i++) {
                //                            bool is_max_c = false;
                //            
                //            
                //                            if (std::fabs(gradient[i]) == max_c) {
                //            
                //                                is_max_c = true;
                //                            }
                //            
                //                            //make the output nice and clean by converting to string and 
                //                            //back to double so we can use scientific format for
                //                            //multi precision types.
                //                            std::stringstream ss;
                //                            ss << gradient[i];
                //                            std::string grad_string = ss.str();
                //            
                //                            ss.str("");
                //                            ss << parameters[i]->GetValue();
                //                            std::string param_string = ss.str();
                //            
                //                            ss.str("");
                //                            ss << max_c;
                //                            std::string max_c_string = ss.str();
                //            
                //                            if ((i % 10) == 0) {
                //                                std::cout << BOLD << std::left << std::setw(10)
                //                                        << "Id"
                //                                        << BOLD << std::left << std::setw(20)
                //                                        << "Parameter"
                //                                        << std::left << std::setw(20) << "Value"
                //                                        << "Gradient["
                //                                        << BOLD << RED << "" << toValue<double>(max_c_string) << "" << DEFAULT
                //                                        << BOLD << std::left << std::setw(10) << "]"
                //                                        << "Bounded" << std::endl;
                //                            }
                //                            std::cout << BLUE << BOLD << std::left << std::setw(10)<<this->active_parameters_m[i]->id_m<< std::left << std::setw(20)
                //                                    << this->active_parameters_m[i]->GetName() << DEFAULT << std::left
                //                                    << std::setw(20) << toValue<double>(param_string) << std::left;
                //            
                //            
                //                            if (is_max_c) {
                //                                std::cout << RED << BOLD << std::setw(31) << toValue<double>(grad_string) << DEFAULT;
                //                            } else {
                //                                if (std::fabs(toValue<double>(grad_string)) <= this->tolerance_m) {
                //                                    number_good++;
                //                                    std::cout << BOLD << GREEN << std::setw(31) << toValue<double>(grad_string) << DEFAULT;
                //                                } else {
                //                                    std::cout << std::setw(31) << toValue<double>(grad_string);
                //                                }
                //                            }
                //            
                //                            if (this->active_parameters_m[i]->IsBounded()) {
                //                                std::cout << "True" << "[" << this->active_parameters_m[i]->GetMinBoundary() << ","
                //                                        << this->active_parameters_m[i]->GetMaxBoundary() << "]" << std::endl;
                //                            } else {
                //                                std::cout << "False" << std::endl;
                //                            }
                //            
                //                        }

                std::cout << io::BOLD << "Phase: " << this->current_phase << " ";
                double percent = 100 * (number_good / (double) this->active_parameters.size());
                //                        double p = 100 * (number_good / (double) this->active_parameters_m.size());
                std::cout << io::BOLD << "|";
                for (double i = 0; i < (percent - 1)*.9; i++) {
                    std::cout << io::GREEN << "=";
                }
                if (percent > 1)
                    std::cout << io::GREEN << "=";
                for (int i = 0; i < (90 - percent); i++) {
                    std::cout << io::RED << "=";
                }
                std::cout << io::DEFAULT;
                std::cout << io::BOLD << "|" << io::DEFAULT;
                std::cout.precision(prec);

                std::cout << "[" << (int) (percent) << "%]\n";
                std::cout << std::endl;
                std::cout << std::fixed;
            } else {
                std::cout << "\n\n";
            }
        }



    };


}




#endif