#include "../../../../Optimization/Optimization2/Optimization.hpp"
#include "../../../../Utilities/IO/StreamedDataFile.hpp"

//#define ADMB_VERSION

#ifdef ADMB_VERSION

template<typename T>
class ThetaLog : public atl::ObjectiveFunction<T> {
    std::vector<T> Y;
    std::vector<atl::Variable<T> > X;
    atl::Variable<T> logr0; // = -2.6032947;
    atl::Variable<T> logtheta; // = 0.7625692;
    atl::Variable<T> logK; // = 6.7250075;
    atl::Variable<T> logQ; // = -4.7496015;
    atl::Variable<T> logR; // = -3.1889239;
    
    
    
public:
    
    ThetaLog() {
        
    }
    
    void Initialize() {
        
        
        atl::StreamedDataFile data;
        data.open("ThetaLog.dat");
        
        int size = 0;
        data >> size;

        this->Y.resize(size);
        this->X.resize(size);
        
        
        for (int i = 0; i < X.size(); i++) {
            data >> Y[i];
        }
        
        logr0.SetName("logr0");
        logtheta.SetName("logtheta");
        logK.SetName("logK");
        logQ.SetName("logQ");
        logR.SetName("logR");
        
        /**
         * Register parameters.
         */
        this->RegisterHyperParameter(logr0);
        this->RegisterHyperParameter(logtheta);
        logK.SetBounds(4.6, 7.6);
        this->RegisterHyperParameter(logK);
        this->RegisterHyperParameter(logQ);
        this->RegisterHyperParameter(logR);
        
        /**
         * Register random variables.
         */
        for (int i = 0; i < X.size(); i++) {
            this->RegisterRandomVariable(X[i]);
        }
    }
    
    void step(atl::Variable<T>& jnll, const atl::Variable<T>& x1, const atl::Variable<T>& x2, const atl::Variable<T>& logr0, const atl::Variable<T>& logK, const atl::Variable<T>& logtheta, const atl::Variable<T>& logQ) {
        atl::Variable<T> var = atl::exp(logQ);
        atl::Variable<T> m = (x1 + atl::exp(logr0) * (static_cast<T>(1.0) - atl::pow(atl::exp(x1) / atl::exp(logK), atl::exp(logtheta))));
        jnll += static_cast<T>(0.5) * (atl::log(static_cast<T>(2.0 * M_PI) * var) + atl::pow((x2 - m),static_cast<T>( 2.0)) / var);
    }
    
    void obs(atl::Variable<T>& jnll, const atl::Variable<T>& x, const atl::Variable<T>& logR, int i) {
        atl::Variable<T> var = atl::exp(logR);
        jnll += static_cast<T> (0.5) * (atl::log(static_cast<T> (2.0 * M_PI) * var)+(atl::pow(x - Y[i], static_cast<T> (2.0))) / var);
    }
    
    virtual const atl::Variable<T> Evaluate() {
        
        atl::Variable<T> f = 0.0;
        
        
        int timeSteps = Y.size();
        
        for (int i = 1; i < timeSteps; i++) {
            step(f, X[i - 1], X[i], logr0, logK, logtheta, logQ);
        }
        
        for (int i = 0; i < timeSteps; i++) {
            obs(f, X[i], logR, i);
        }
        return f;
    }
    
    void Report() {
        std::ofstream out;
        out.open("ThetaLog.par");
        out << "logr0 = " << logr0 << "\n";
        out << "logtheta = " << logtheta << "\n";
        out << "logK = " << logK << "\n";
        out << "logQ = " << logQ << "\n";
        out << "logR = " << logR << "\n";
        out << "\nRandom Effects:\n";
        out << "X:\n";
        for (int i = 0; i < X.size(); i++) {
            out << X[i] << " ";
        }
    }
    
};

#else


template<typename T>
class ThetaLog : public atl::ObjectiveFunction<T> {
    std::vector<T> Y;
    std::vector<atl::Variable<T> > X;
    atl::Variable<T> logr0; //  = -2.6032947;
    atl::Variable<T> logtheta; //  = 0.7625692;
    atl::Variable<T> logK = T(6.0); //  = 6.7250075;
    atl::Variable<T> logQ; //  = -4.7496015;
    atl::Variable<T> logR; // = -3.1889239;
    
    
public:
    
    ThetaLog() {
        
    }
    
    void Initialize() {
        
        atl::StreamedDataFile data;
        data.open("ThetaLog.dat");
        
        int size = 0;
        data >> size;
        std::cout << "size = " << size;
   
        this->Y.resize(size);
        this->X.resize(size);
        
        
        for (int i = 0; i < X.size(); i++) {
            double tmp;
            data >> Y[i];
//            Y[i] = tmp;
        }
        
        
        this->RegisterHyperParameter(logr0);
        logr0.SetName("logr0");
        this->RegisterHyperParameter(logtheta);
        logtheta.SetName("logtheta");
        this->RegisterHyperParameter(logK);
        logK.SetName("logK");
        this->RegisterHyperParameter(logQ);
        logQ.SetName("logQ");
        this->RegisterHyperParameter(logR);
        logR.SetName("logR");
        
        for (int i = 0; i < X.size(); i++) {
            this->RegisterRandomVariable(X[i]);
        }
    }
    
    const atl::Variable<T> dnorm(const atl::Variable<T> x,
                                 const atl::Variable<T> mean,
                                 const atl::Variable<T> sd, int give_log = 0) {
        if (sd.GetValue() == 0.0) {
            throw std::overflow_error("Divide by zero exception");
        }
        
        atl::Variable<T> logres = -1.0*atl::log(T(sqrt(2*M_PI))*sd)-T(.5)*atl::pow((x-mean)/sd,2.0);
        if (give_log)return logres;
        else return atl::exp(logres);
    }
    
    virtual const atl::Variable<T> Evaluate() {
        
        atl::Variable<T> r0 = atl::exp(logr0);
        atl::Variable<T> theta = atl::exp(logtheta);
        atl::Variable<T> K = atl::exp(logK);
        atl::Variable<T> Q = atl::exp(logQ);
        atl::Variable<T> R = atl::exp(logR);
        
        
        int timeSteps = Y.size();
        atl::Variable<T> ans = T(0.0);
        
        atl::Variable<T> sqrtq = atl::sqrt(Q);
        for (int i = 1; i < timeSteps; i++) {
            atl::Variable<T> m = X[i - 1] + r0 * (static_cast<T>(1.0) - std::pow(std::exp(X[i - 1]) / K, theta));
            ans -= this->dnorm(X[i], m, sqrtq, true);
        }
        atl::Variable<T> sqrtr = atl::sqrt(R);
        for (int i = 0; i < timeSteps; i++) {
            ans -= this->dnorm(atl::Variable<T>(Y[i]), X[i], sqrtr, true);
        }
        
        return ans;
    }
    
    void Report() {
        std::ofstream out;
        out.open("ThetaLog.par");
        out << "logr0 = " << logr0 << "\n";
        out << "logtheta = " << logtheta << "\n";
        out << "logK = " << logK << "\n";
        out << "logQ = " << logQ << "\n";
        out << "logR = " << logR << "\n";
        out << "\nRandom Effects:\n";
        out << "X:\n";
        for (int i = 0; i < X.size(); i++) {
            out << X[i] << " ";
        }
    }
    
};

#endif

int main(int argc, char** argv) {
    
    typedef double  REAL;

    //create the objective function 
    ThetaLog<REAL> objective_function;
    
    //initialize the objective function
    objective_function.Initialize();
    
    //create an instance of a L-BFGS minimizer
    atl::PortMinimizer<REAL> fm;
    
    //set the objective function
    fm.SetObjectiveFunction(&objective_function);
    
    //run the minimizer
    fm.Run();
//    
//    //dump the results
//    objective_function.Report();
//    
//    std::cout<<objective_function.GetObjectiveFunctionStatistics()<<"\n";
    return 0;
}

